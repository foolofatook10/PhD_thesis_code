#load libraries
library(tidyverse)

#set working directory
setwd("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/R_scripts")

#read in common variables
source("common_variables.R")

myTheme <- theme_classic()+
  theme(axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 18),
        legend.title = element_blank())


#read in read counts summaries----
data_list <- list()
for (sample in RPF_sample_names) {
  read_csv(file = file.path(parent_dir, "logs", paste0(sample, "_read_counts.csv"))) %>%
    mutate(sample = rep(sample)) %>%
    inner_join(RPF_sample_info, by = "sample") -> data_list[[sample]]
}

RPF_counts <- do.call("rbind", data_list)
#write_csv(file = file.path(parent_dir, "logs/RPF_reads_summary.csv"), RPF_counts)

RPF_counts$condition2 <- factor(RPF_counts$condition2, levels = c("M", "D30")) 

# read in publication theme
publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(size = 16, color = "black"),
            axis.text.y = element_text(size = 14, color = "black"),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = 16),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

#input read counts----
RPF_counts %>%
  ggplot(aes(x = condition, y = cutadapt_in, fill = replicate))+
  facet_wrap(~condition2) +
  geom_col(position = position_dodge())+
  scale_fill_brewer(palette = "Set1") +
  ylab("all counts")+
  ggtitle("RPFs")+
  publication_theme() -> RPFs_counts_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/RPF_all_counts.png"), width = 500, height = 300)
print(RPFs_counts_plot)
dev.off()

#trimmed percentages----
RPF_counts %>%
  mutate(cutadapt_perc = (cutadapt_out / cutadapt_in) * 100) %>%
  ggplot(aes(x = condition, y = cutadapt_perc, fill = replicate))+
  facet_wrap(~condition2) +
  scale_fill_brewer(palette = "Set1") +
  geom_col(position = position_dodge())+
  ylab("trimmed reads %")+
  ylim(c(0,100))+
  ggtitle("RPFs")+
  publication_theme() -> RPFs_cutadapt_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/RPFs_cutadapt_counts.png"), width = 500, height = 300)
print(RPFs_cutadapt_plot)
dev.off()

#alignment percentages----
read_alignments <-
RPF_counts %>%
  mutate(rRNA_perc = rRNA_out / UMI_clipped_out * 100,
         pc_perc = (pc_out / UMI_clipped_out) * 100,
         tRNA_perc = (tRNA_out / UMI_clipped_out) * 100,
         unaligned_perc = ((pc_in - pc_out) / UMI_clipped_out) * 100,
         sample = str_replace(sample, "_RPFs_", "\n")) %>%
  select(sample, rRNA_perc, pc_perc, tRNA_perc, unaligned_perc) %>%
  gather(key = alignment, value = percentage, rRNA_perc, pc_perc, tRNA_perc, unaligned_perc) %>%
  mutate(alignment = factor(alignment, levels = c("unaligned_perc", "tRNA_perc", "pc_perc", "rRNA_perc"), 
                            labels = c("other", "tRNA", "pc", "rRNA"), ordered = T)) 

read_alignments$condition2 <-
  str_remove(
  str_remove(
    str_extract(read_alignments$sample, pattern = "_.*_"),
    "_"),"_")

read_alignments$condition2 <- factor(read_alignments$condition2, levels = c("M", "D30"))

read_alignments$new_sample <- str_replace(read_alignments$sample,
                                          pattern = "_.*_",
                                          replacement = "_")

read_alignments$new_sample <- factor(read_alignments$new_sample, levels = c("REP2_Tot",
                                                                    "REP2_CNOT3",
                                                                    "REP3_Tot",
                                                                    "REP3_CNOT3",
                                                                    "REP4_Tot",
                                                                    "REP4_CNOT3"),
                                     labels = c("REP2 Tot",
                                                "REP2 CNOT3",
                                                "REP3 Tot",
                                                "REP3 CNOT3",
                                                "REP4 Tot",
                                                "REP4 CNOT3"))

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
            #axis.title = element_text(face = "bold",size = rel(1)),
            #axis.title.y = element_text(angle=90,vjust =2),
            #axis.title.x = element_text(vjust = -0.2),
            #axis.text.x = element_text(size = 16, color = "black"),
            #axis.text.y = element_text(size = 16, color = "black"),
            #axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = 10),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(1, "cm"),
            legend.margin = unit(0.5, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

read_alignments$alignment <- factor(read_alignments$alignment,
                                    levels = c("other", "tRNA", "pc", "rRNA"),
                                    labels = c("Other", "tRNA", "Protein coding", "rRNA"))

read_alignments$new_sample <- factor(read_alignments$new_sample,
                                     levels = c("REP2 Tot", "REP2 CNOT3",
                                                "REP3 Tot", "REP3 CNOT3",
                                                "REP4 Tot", "REP4 CNOT3"),
                                     labels = c("REP1 Total", "REP1 CNOT3 IP'd",
                                                "REP2 Total", "REP2 CNOT3 IP'd",
                                                "REP3 Total", "REP3 CNOT3 IP'd"))


read_alignments %>% 
  filter(condition2 == "M") %>% 
  ggplot(aes(x = new_sample, y = percentage, fill = alignment))+
  #facet_wrap(~condition2, scales="free_x", nrow = 2) +
  geom_col()+
  scale_fill_brewer(palette = "Pastel1") +
  xlab("sample")+
  ylab("% Aligments")+
  ggtitle("RPF Alignments")+
  publication_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank())-> RPF_aligments_plot

png(filename = file.path("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/SelRiboseq_RPFs_alignments.png"), res = 300,  width = 1500, height = 1250)
print(RPF_aligments_plot)
dev.off()

RPF_counts$condition <- factor(RPF_counts$condition, levels = c("Tot", "CNOT3"))

#unique pc percentages----
RPF_counts %>%
  mutate(unique_perc = (deduplication_out / deduplication_in) * 100) %>%
  ggplot(aes(x = condition, y = unique_perc, fill = replicate))+
  geom_col(position = position_dodge())+
  facet_wrap(~condition2) +
  scale_fill_brewer(palette = "Set1") +
  ylab("Unique pc reads %")+
  ylim(c(0,100))+
  ggtitle("RPFs")+
  publication_theme() -> RPFs_deduplication_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/RPFs_deduplication_counts.png"), width = 500, height = 300)
print(RPFs_deduplication_plot)
dev.off()

#final pc counts----
RPF_counts %>%
  filter(condition2 != "D30")  %>%
  mutate(deduplication_out2 = format(deduplication_out, scientific = T, digits = 3)) %>% 
  mutate(replicate = factor(case_when(replicate == "2" ~ "1",
                               replicate == "3" ~ "2",
                               replicate == "4" ~ "3"),
                            levels = c("1", "2", "3"))) %>% 
  ggplot(aes(x = condition, y = deduplication_out, fill = replicate, label = deduplication_out2)) +
  geom_col(position = position_dodge())+
  geom_text(aes(label = deduplication_out2), position = position_dodge(.9), angle = 90, hjust = 0.5) +
  #facet_wrap(~condition2) +
  scale_fill_brewer(palette = "Pastel2") +
  ylab("Unique protein coding reads")+
  ggtitle("RPFs")+
  scale_y_continuous(limits = c(0,30000000), breaks = c(seq(0, 30000000,5000000))) +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        legend.position = "right",
        legend.direction = "vertical")-> RPFs_pc_count_plot

png(filename = file.path("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/SelRiboseq_RPFs_unique_pc_reads.png"), res = 300,  width = 1750, height = 1250)
print(RPFs_pc_count_plot)
dev.off()
