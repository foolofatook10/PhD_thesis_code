library(tidyverse)
library(data.table)

colours_for_ggplot <- c("#1B9E77", "#D95F02", "#7570B3")

# read in publication theme
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

# set working directory
setwd("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL")

# read in transcript half lives
t_half <- read_csv(file = file.path("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/CNOTdata_mRNAhalflives.csv"))

# read in features table
features <- as_tibble(fread("\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/features.csv", drop = "V1")) 

# read in deseq2
deseq2 <- read_csv(file = file.path("Analysis/DESeq2_output/RPFs_CNOT3_DEseq2_apeglm_LFC_shrinkage_M.csv"))

# join deseq2 output with half life data
data <- deseq2 %>% inner_join(t_half, by = c("gene_sym" = "Gene_ID"))

temp <- 
  data %>% mutate(grouping = case_when(padj < 0.05 & log2FoldChange > 0 ~ "enriched",
                                       padj < 0.05 & log2FoldChange < 0 ~ "depleted",
                                       .default = "unchanged")) 

table(temp$grouping)

temp$grouping <- factor(temp$grouping, levels = c("depleted", "unchanged", "enriched"),
                        labels = c("Dep.", "Unch.", "Enr."))

panel_E <- 
  temp %>%   
  ggplot(aes(x = grouping, y = siControl_half_life, fill = grouping)) +
  geom_violin() +
  geom_boxplot(width  = 0.2) +
  scale_y_log10() +
  scale_fill_manual(values = colours_for_ggplot) +
  ylab("mRNA half-life (mins)") +
  stat_summary(
    fun = "mean", 
    geom = "point", 
    shape = 21,       # Circle shape (21 allows fill and color)
    size = 3,         # Size of the circle
    color = "black",  # Border color
    fill = "black"      # Inside color
  ) +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

panel_F <- 
  temp %>%   
  ggplot(aes(x = grouping, y = log2FC_halflife, fill = grouping)) +
  geom_violin() +
  geom_boxplot(width  = 0.2) +
  scale_fill_manual(values = colours_for_ggplot) +
  ylab("Log2FC in mRNA half life\n(siCNOT1-NTC)")  +
  stat_summary(
    fun = "mean", 
    geom = "point", 
    shape = 21,       # Circle shape (21 allows fill and color)
    size = 3,         # Size of the circle
    color = "black",  # Border color
    fill = "black"      # Inside color
  ) +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

panel_C <- 
  temp %>% 
  ggplot(aes(x = log2FoldChange, y = siControl_half_life)) +
  geom_point() +
  scale_y_log10() +
  ggpubr::stat_cor(method="pearson",
                   inherit.aes = F,
                   data = temp,
                   aes(x = log2FoldChange, y = siControl_half_life),
                   colour = "red") +
  ylab("mRNA half-life (mins)") +
  xlab("Log2 Fold Enrichment\n(CNOT3 IP - Tot. RPFs)") +
  publication_theme()

# stats in control conditions
aov_obj <- 
  aov(formula = siControl_half_life ~ grouping, data =  temp)

TukeyHSD(aov_obj)

panel_D <- 
  temp %>% 
  ggplot(aes(x = log2FoldChange, y = log2FC_halflife)) +
  geom_point() +
  #scale_y_log10() +
  ggpubr::stat_cor(method="pearson",
                   inherit.aes = F,
                   data = temp,
                   aes(x = log2FoldChange, y = log2FC_halflife),
                   colour = "red") +
  ylab("Log2FC in mRNA half-life\n(siCNOT1-NTC)") +
  xlab("Log2 Fold Enrichment\n(CNOT3 IP - Tot. RPFs)") +
  publication_theme()

# stats for change in mRNA half life
aov_obj <- 
  aov(formula = log2FC_halflife ~ grouping, data =  temp)

TukeyHSD(aov_obj)


# save files

path_to_thesis_figures <- "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures"

res_pic = 300
height_pic = 1500
width_pic = 1500

png(filename = file.path(path_to_thesis_figures, "CNOT3SelRP_enrichment_versus_mRNA_half_life.png"),
    res = res_pic,
    height = height_pic,
    width = width_pic)
print(panel_C)
dev.off()

png(filename = file.path(path_to_thesis_figures, "CNOT3SelRP_enrichment_versus_change_in_mRNA_half_life_following_siCNOT1.png"),
    res = res_pic,
    height = height_pic,
    width = width_pic)
print(panel_D)
dev.off()

####
res_pic = 300
height_pic = 1250
width_pic = 1250

png(filename = file.path(path_to_thesis_figures, "CNOT3SelRP_enrichment_versus_mRNA_half_life_boxplots.png"),
    res = res_pic,
    height = height_pic,
    width = width_pic)
print(panel_E)
dev.off()

png(filename = file.path(path_to_thesis_figures, "CNOT3SelRP_enrichment_versus_change_in_mRNA_half_life_following_siCNOT1_boxplots.png"),
    res = res_pic,
    height = height_pic,
    width = width_pic)
print(panel_F)
dev.off()


# end of save files
tpms <- 
  read_csv("\\\\data.beatson.gla.ac.uk/data/R11/external_sequencing_data/Gillen_2021_RiboSeq/Analysis/DESeq2_output/tpms.csv")

tpms_long <- 
  tpms %>% 
  gather(key = "sample", value = tpm, Ctrl_1_Totals:CNOT1_3_Totals) 

samples_split <- str_split(tpms_long$sample, "_")
condition <- unlist(lapply(samples_split, function(x){x[[1]]}))

tpms_long$condition <- condition

mean_tpms <- 
  tpms_long %>% 
  group_by(transcript, condition) %>% 
  summarise(mean_tpm = mean(tpm)) %>% 
  filter(mean_tpm > 0)

half_lives <- read_csv("\\\\data.beatson.gla.ac.uk/data/R11/external_sequencing_data/Gillen_2021_halflife/Analysis/DESeq2_output/DEseq2_apeglm_LFC_shrinkage.csv")

half_lives <- 
  half_lives %>% 
  mutate(siCNOT1 = case_when(log2FoldChange > 0 & padj < 0.05 ~ "stabilised",
                             log2FoldChange <0 & padj  < 0.05 ~ "destabilised",
                             .default = "unchanged")) %>% 
  filter(siCNOT1 != "unchanged")

data_Ctrl <- 
  inner_join(deseq2 %>% select(gene_sym, transcript, log2FoldChange, padj),
             mean_tpms,
             by = c("transcript")) %>% 
  mutate(tpm_log10 = log10(mean_tpm)) %>% 
  filter(condition == "Ctrl") %>% 
  inner_join(
    half_lives %>% select(transcript, siCNOT1), by = "transcript")

data_Ctrl <- 
  data_Ctrl %>% 
  mutate(class = case_when(log2FoldChange > 0 & padj < 0.05 ~ "Enriched",
                           log2FoldChange < 0 & padj < 0.05 ~ "Depleted",
                           .default = "Unchanged")) 

library(scales)

panel_A <- 
  data_Ctrl %>%   
  ggplot(aes(x = log2FoldChange, y = mean_tpm, colour = siCNOT1)) +
  geom_point() +
  # ggpubr::stat_cor(method="pearson", 
  #                  inherit.aes = F, 
  #                  data = data_Ctrl, 
  #                  aes(x = log2FoldChange, y = mean_tpm)) +
  scale_y_log10(labels = scales::label_comma()) +
  #ggtitle("Ctrl TPMs") +
  xlab("Log2 Fold Enrichment\n(CNOT3 IP - Tot RPFs)") +
  ylab("Transcript expression\n(mean TPM)") +
  #scale_y_continuous(limits = c(-2, 4.5)) +
  publication_theme() +
  theme(legend.title = element_blank(),
        legend.position = "none")

data_Ctrl$class <- factor(data_Ctrl$class, 
                          levels = c("Depleted", "Unchanged", "Enriched"),
                          labels = c("Dep.", "Unch.", "Enr."))
panel_B <- 
  data_Ctrl %>% 
  ggplot(aes(x = class, y = mean_tpm, fill = class)) +
  scale_fill_brewer(palette = "Dark2") +
  geom_violin() +
  geom_boxplot(width = 0.2)  +
  scale_y_log10(labels = scales::label_comma())+
  xlab("Log2 Fold Enrichment\n(CNOT3 IP - Tot RPFs)") +
  ylab("Transcript expression\n(mean TPM)") +
  stat_summary(
    fun = "mean", 
    geom = "point", 
    shape = 21,       # Circle shape (21 allows fill and color)
    size = 3,         # Size of the circle
    color = "black",  # Border color
    fill = "black"      # Inside color
  )+
  #ylab("Log10 TPMs") +
  #ggtitle("Ctrl") +
  publication_theme() +
  theme(legend.position = "none")

anova_object <- aov(formula = mean_tpm ~ class, data = data_Ctrl)

TukeyHSD(anova_object)

#####
res_pic = 300
height_pic = 1500
width_pic = 1500

png(filename = file.path(path_to_thesis_figures, "CNOT3SelRP_enrichment_versus_mRNA_abundance.png"),
    res = res_pic,
    height = height_pic,
    width = width_pic)
print(panel_A)
dev.off()

res_pic = 300
height_pic = 1500
width_pic = 1500

png(filename = file.path(path_to_thesis_figures, "CNOT3SelRP_enrichment_versus_mRNA_abundance_boxplots.png"),
    res = res_pic,
    height = height_pic,
    width = width_pic)
print(panel_B)
dev.off()
