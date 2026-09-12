#load libraries
library(tidyverse)
library(knitr)

#read in common variables
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/R_scripts/common_variables.R")

myTheme <- theme_classic()+
  theme(axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 18),
        legend.title = element_blank())

#create themes----
publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = 20),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(size = 18, color = "black"),
            axis.text.y = element_text(size = 18, color = "black"),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = 18),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.75, "cm"),
            legend.margin = unit(0.5, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(5,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


#read in read counts summaries----
data_list <- list()
for (sample in RPF_sample_names) {
  read_csv(file = file.path(parent_dir, "logs", paste0(sample, "_read_counts.csv"))) %>%
    mutate(sample = rep(sample)) %>%
    inner_join(RPF_sample_info, by = "sample") -> data_list[[sample]]
}

RPF_counts <- do.call("rbind", data_list)
write_csv(file = file.path(parent_dir, "logs/RPF_reads_summary.csv"), RPF_counts)

#input read counts----
RPF_counts %>%
  mutate(cutadapt_in2 = format(cutadapt_in, scientific = T, digits = 3)) %>% 
  ggplot(aes(x = condition, y = cutadapt_in, fill = replicate, label = cutadapt_in2))+
  geom_col(position = position_dodge())+
  geom_text(aes(label = cutadapt_in2), position = position_dodge(.9), angle = 90, hjust = 1.5) +
  ylab("Raw reads")+
  ggtitle("RPFs")+
  scale_fill_brewer(palette = "Pastel1") +
  publication_theme() +
  theme(legend.position = "none",
        legend.direction = "vertical",
        axis.title.x = element_blank())-> RPFs_counts_plot

png(filename = file.path("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/Riboseq_RPFs_raw_reads.png"), res = 300,  width = 1200, height = 1200)
print(RPFs_counts_plot)
dev.off()

# png(filename = file.path(parent_dir, "plots/read_counts_summary/RPF_all_counts.png"), width = 500, height = 300)
# print(RPFs_counts_plot)
# dev.off()

#trimmed percentages----
RPF_counts %>%
  mutate(cutadapt_perc = (cutadapt_out / cutadapt_in) * 100) %>%
  ggplot(aes(x = condition, y = cutadapt_perc, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("trimmed reads %")+
  ylim(c(0,100))+
  ggtitle("RPFs")+
  myTheme -> RPFs_cutadapt_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/RPFs_cutadapt_counts.png"), width = 500, height = 300)
print(RPFs_cutadapt_plot)
dev.off()

#alignment percentages----
alignments_processed <- 
RPF_counts %>%
  mutate(rRNA_perc = rRNA_out / UMI_clipped_out * 100,
         pc_perc = (pc_out / UMI_clipped_out) * 100,
         tRNA_perc = (tRNA_out / UMI_clipped_out) * 100,
         unaligned_perc = ((pc_in - pc_out) / UMI_clipped_out) * 100,
         sample = str_replace(sample, "_RPFs_", "\n")) %>%
  select(sample, rRNA_perc, pc_perc, tRNA_perc, unaligned_perc) %>%
  gather(key = alignment, value = percentage, rRNA_perc, pc_perc, tRNA_perc, unaligned_perc) %>%
  mutate(alignment = factor(alignment, levels = c("unaligned_perc", "tRNA_perc", "pc_perc", "rRNA_perc"), labels = c("Other", "tRNA", "Protein coding", "rRNA"), ordered = T)) %>%
  mutate(sample = str_remove(str_replace_all(sample, "_", " "), "\\sM$"))

alignments_processed$sample <- 
  str_replace(alignments_processed$sample, 
            pattern = "siRPS25", 
            replacement = "sieS25")

alignments_processed$sample %>% unique()

alignments_processed$condition <- str_remove(alignments_processed$sample, "REP\\d\\s")

alignments_processed$replicate <- as.numeric(str_remove(str_extract(alignments_processed$sample, "REP\\d"), "REP"))

alignments_processed %>% 
  ggplot(aes(x = sample, y = percentage, fill = alignment))+
  geom_col()+
  scale_fill_brewer(palette = "Pastel1") +
  xlab("sample")+
  ylab("Read Alignments (%)")+
  ggtitle("RPFs")+
  #scale_x_discrete(labels = alignments_processed$condition) +
  publication_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank(),
        axis.title.x = element_blank())-> RPF_aligments_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/RPF_aligments_plot.png"), width = 900, height = 300)
print(RPF_aligments_plot)
dev.off()

png(filename = file.path("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/Riboseq_RPFs_percentage_alignment.png"), res = 300,  width = 1750, height = 1000)
print(RPF_aligments_plot)
dev.off()

#unique pc percentages----
RPF_counts$condition <- factor(RPF_counts$condition,
                               levels = c("NTC", "siCNOT3", "siRPS25"),
                               labels = c("NTC", "siCNOT3", "sieS25"))

RPF_counts %>%
  mutate(unique_perc = (deduplication_out / deduplication_in) * 100) %>%
  ggplot(aes(x = condition, y = unique_perc, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("Unique pc reads %")+
  scale_fill_brewer(palette = "Pastel2") +
  ylim(c(0,100))+
  ggtitle("RPFs")+
  publication_theme() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        axis.title.x = element_blank())-> RPFs_deduplication_plot

png(filename = file.path("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/RPFs_deduplication_counts.png"), res = 300,  width = 1750, height = 1250)
print(RPFs_deduplication_plot)
dev.off()


#final pc counts----
RPF_counts %>%
  mutate(deduplication_out2 = format(deduplication_out, scientific = T, digits = 3)) %>% 
  ggplot(aes(x = condition, y = deduplication_out, fill = replicate, label = deduplication_out2))+
  geom_col(position = position_dodge())+
  geom_text(aes(label = deduplication_out2), position = position_dodge(.9), angle = 90, hjust = 0.5) +
  scale_fill_brewer(palette = "Pastel2") +
  ylab("Unique protein coding reads")+
  ggtitle("RPFs")+
  scale_y_continuous(limits = c(0,7000000), breaks = c(seq(0, 7000000, 1000000))) +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        legend.position = "right",
        legend.direction = "vertical")-> RPFs_pc_count_plot


png(filename = file.path("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/Riboseq_RPFs_unique_pc_reads.png"), res = 300,  width = 1750, height = 1250)
print(RPFs_pc_count_plot)
dev.off()
