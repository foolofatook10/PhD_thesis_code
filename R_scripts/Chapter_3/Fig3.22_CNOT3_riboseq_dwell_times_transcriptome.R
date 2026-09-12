library(tidyverse)
library(data.table)

dts <- read_csv("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/dwell_times_CNOT3_riboseq_entire_transcriptome.csv") %>% 
  select(-`...1`)

# deseq <- 
#   read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/Analysis/DESeq2_output/RPFs_CNOT3_DEseq2_apeglm_LFC_shrinkage_M_28nts.csv") %>% 
#   filter(log2FoldChange > 0 & padj < 0.05) %>% 
#   pull(transcript)
# 
# dts <- dts %>% filter(transcript %in% deseq)

dts_summarised <- 
  dts %>% 
  group_by(sample_id, codon) %>% 
  summarise(`-4` = sum(minus4),
            `-3` = sum(minus3),
            `E` = sum(E_site),
            `P` = sum(P_site),
            `A` = sum(A_site),
            `+1` = sum(plus1),
            `+2` = sum(plus2),
            total = sum(total))

dts_summarised_split <- split(dts_summarised, dts_summarised$sample_id)

downstream_processing <- function(df){
  
  df <-
    df %>%
    filter(codon != "tga" & codon != "tag" & codon != "taa") %>%
    mutate(minus_4 = (log2(`-4`)) - log2(total / 7),
           minus_3 = (log2(`-3`)) - log2(total / 7),
           E_site = (log2(E)) - log2(total / 7),
           P_site = (log2(P)) - log2(total / 7),
           A_site = (log2(A)) - log2(total / 7),
           plus_1 = (log2(`+1`)) - log2(total / 7),
           plus_2 = (log2(`+2`)) - log2(total / 7)) %>%
    select(codon, minus_4, minus_3, E_site, P_site, A_site, plus_1, plus_2, sample_id) %>%
    dplyr::rename(sample = sample_id)
  
  return(df)
}



data <- 
  do.call("rbind", lapply(dts_summarised_split, downstream_processing)) %>% 
  ungroup()

#gather the data into tidy format----
data %>%
  pull(codon) %>%
  unique() -> codons

RPF_sample_names <- unique(data$sample)

gathered_list <- list()
for (sample in RPF_sample_names) {
  for (codon in codons) {
    data[data$codon == codon & data$sample == sample,] %>%
      select(-sample) %>%
      gather(key = codon, value = freq) %>%
      rename(position = codon) %>%
      mutate(codon = rep(codon),
             sample = factor(rep(sample))) -> gathered_list[[paste(codon, sample, sep = "_")]]
  }
}
gathered_data <- do.call("rbind", gathered_list)
# 
gathered_data %>%
  mutate(position = as.numeric(case_when(position == "minus_4" ~ -4,
                                         position == "minus_3" ~ -3,
                                         position == "E_site" ~ -2,
                                         position == "P_site" ~ -1,
                                         position == "A_site" ~ 0,
                                         position == "plus_1" ~ 1,
                                         position == "plus_2" ~ 2)),
         wobble = factor(str_sub(codon, 3,3))) -> plot_data

for (sample in RPF_sample_names) {
  plot_title <- str_replace_all(sample, "_", " ")
  plot_title <- str_remove(plot_title, " RPFs")
  
  plot_data[plot_data$sample == sample,] %>%
    ggplot(aes(x = position, y = freq, colour = wobble))+
    geom_point(size = 2)+
    stat_summary(fun=mean, geom="line", size = 1)+
    ylab("normalised codon frequency")+
    scale_x_continuous(limits = c(-4,2), breaks = -4:2)+
    xlab("Codon position (0 = A-site)")+
    ggtitle(plot_title) -> codon_plot
  
  print(codon_plot)
  
}
# 
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "-3") %>% arrange(-freq)})
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "-2") %>% arrange(-freq)})
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "-1") %>% arrange(-freq)})
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "0") %>% arrange(-freq)})
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "1") %>% arrange(-freq)})
# 
# 
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "-1") %>% arrange(-freq)})
# 
#        
mean_data <- 
  plot_data %>% 
  mutate(condition = unlist(lapply(str_split(sample, "_"),function(x){x[[2]]}))) %>% 
  group_by(condition, position, codon, wobble) %>% 
  summarise(mean_freq = mean(freq)) 

# 
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
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

mean_data_final <- 
  mean_data %>% 
  mutate(to_col = case_when(codon == "ggc" & position == -1 ~ "Gly-GGC",
                            codon == "ggt" & position == -1 ~ "Gly-GGT",
                            codon == "gga" & position == -1 ~ "Gly-GGA",
                            codon == "ggg" & position == -1 ~ "Gly-GGG",
                            codon == "gaa" & position == 0 ~ "Glu-GAA",
                            codon == "tgg" & position == 0 ~ "Trp-TGG",
                            codon == "cca" & position == -2 ~ "Pro-CCA",
                            codon == "cct" & position == -2 ~ "Pro-CCT",
                            codon == "ccg" & position == -2 ~ "Pro-CCG",
                            codon == "ccc" & position == -2 ~ "Pro-CCC",.default = "other")) 

mean_data_final$to_col <- 
  factor(mean_data_final$to_col, levels = c("Pro-CCG", "Pro-CCA", "Pro-CCT", "Pro-CCC", "Gly-GGC", "Gly-GGT", "Gly-GGA","Gly-GGG", "Glu-GAA", "Trp-TGG", "other"))


mean_data_final$condition <- factor(mean_data_final$condition,
                                    levels = c("NTC", "siCNOT3"))
set.seed(123456789)


mean_data_final %>% 
  filter(position == 0 & condition == "siCNOT3") %>% 
  arrange(-mean_freq)

mean_data_final

final_gplot <- 
  mean_data_final %>% 
  ggplot(aes(x = position, y = mean_freq, colour = to_col)) +
  geom_point(size = 2, position = position_jitter(width = 0.15)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#984EA3", "#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F","yellow", "#40E0D0", "#E5C494", "#d3d3d3")) +
  facet_wrap(~condition, nrow = 2) +
  scale_x_continuous(limits = c(-4,2), 
                     breaks = -4:2, 
                     labels = c("-2", "-1", "E", "P", "A", "+1", "+2"))+
  #xlab("Codon position (0 = A-site)") +
  ylab("Log2 enrichment score") +
  publication_theme() +
  ylim(-1,1) +
  ggtitle("CNOT3 Riboseq") +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank())


png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/siCNOT3_dwell_times_transcriptome.png",
    res = 300, height = 1750, width = 1750)
print(final_gplot)
dev.off()

test <- 
mean_data_final %>% 
  filter(position == 0) %>% 
  spread(condition,mean_freq) 

cor(test$NTC, test$siCNOT3)
?cor

A_site_correlation <- 
  test %>%   
  ggplot(aes(x = NTC, y = siCNOT3, colour = wobble, label = codon)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = Inf, show.legend = F) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")) +
  ggtitle("A site dwell time") +
  publication_theme()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/siCNOT3_dwell_times_correlation_A.png",
    res = 300, height = 1200, width = 1200)
print(A_site_correlation)
dev.off()

mean_data_final %>% 
  filter(position == -1) %>% 
  spread(condition,mean_freq) %>%   
  ggplot(aes(x = NTC, y = siCNOT3, colour = wobble, label = codon)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = Inf) +
  ggtitle("P site dwell time") +
  publication_theme()

mean_data_final %>% 
  filter(position == -2) %>% 
  spread(condition,mean_freq) %>%   
  ggplot(aes(x = NTC, y = siCNOT3, colour = wobble, label = codon)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = Inf, show.legend = F) +
  ggtitle("E site dwell time") +
  publication_theme()

diffs <- 
  mean_data_final %>% 
  spread(key = condition, value = mean_freq) %>% 
  mutate(differential = siCNOT3-NTC) 

diffs$position %>% unique()

diffs_coloured <- 
  diffs %>% 
  mutate(to_col2 = case_when(codon == "agg" & position == -1 ~ "Arg-AGG",
                             codon == "cga" & position == -1 ~ "Arg-CGA",
                             codon == "cgg" & position == -1 ~ "Arg-CGG",
                             codon == "tcg" & position == -1 ~ "Ser-TCG",
                             codon == "tca" & position == -1 ~ "Ser-TCA",
                             codon == "tgt" & position == -1 ~ "Cys-TGT",
                             codon == "cta" & position == 0 ~ "Leu-CTA",
                             codon =="act" & position == 0 ~ "Thr-ACT",
                             .default = "other"))

diffs_coloured$to_col2 <- factor(diffs_coloured$to_col2,
                                 levels = c("Arg-AGG", "Arg-CGA", "Arg-CGG",
                                            "Ser-TCG", "Ser-TCA", "Cys-TGT",
                                            "Leu-CTA", "Thr-ACT", "other")) 

diffs_coloured %>% filter(position == 0) %>% arrange(-differential)

set.seed(654321)

diff_gplot <-  
  ggplot(diffs_coloured, aes(x = position, y = differential, colour = to_col2)) +
  geom_point(size = 2, position = position_jitter(width = 0.1)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#984EA3", "#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F", "#40E0D0", "#d3d3d3")) +
  #facet_wrap(~condition, ncol = 2) +
  scale_x_continuous(limits = c(-4,2), 
                     breaks = -4:2, 
                     labels = c("-2", "-1", "E", "P", "A", "+1", "+2"))+
  #xlab("Codon position (0 = A-site)") +
  ylab("Differential enrichment") +
  ggtitle("Deseq enriched") +
  publication_theme() +
  ylim(-0.5,0.5) +
  #coord_fixed() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank())


png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/differential_gplot_deseq_enriched.png",
    res = 300, height = 1250, width = 1500)
print(diff_gplot)
dev.off()
