#load libraries
library(tidyverse)

#read in common variables
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/R_scripts/common_variables.R")

#create a variable for what the treatment is----
treatment <- "siRPS25"

#themes----
mytheme <- theme_classic()+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 16))

#read in DESeq2 output----
totals <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("Totals_", treatment, "_DEseq2_apeglm_LFC_shrinkage.csv")))
RPFs <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("RPFs_", treatment, "_DEseq2_apeglm_LFC_shrinkage.csv")))
TE <- read_csv((file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("TE_", treatment, "_DEseq2.csv"))))

TE %>% 
  filter(padj <0.1) 

#plot volcanos----
RPFs %>%
  filter(!(is.na(padj))) %>%
  mutate(sig = factor(case_when(padj < 0.05 ~ "*",
                                padj >= 0.05 ~ "NS"))) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = sig))+
  geom_point(alpha = 0.5)+
  mytheme+
  xlab("log2FC")+
  ylab("-log10(padj)")+
  ggtitle(paste(treatment, "RPFs")) -> RPFs_volcano

png(filename = file.path(parent_dir, "plots/DE_analysis", paste0("RPFs_", treatment, "_volcano.png")), width = 300, height = 400)
print(RPFs_volcano)
dev.off()

totals %>%
  filter(!(is.na(padj))) %>%
  mutate(sig = factor(case_when(padj < 0.05 ~ "*",
                                padj >= 0.05 ~ "NS"))) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = sig))+
  geom_point(alpha = 0.5)+
  mytheme+
  xlab("log2FC")+
  ylab("-log10(padj)")+
  ggtitle(paste(treatment, "Cytoplasmic RNA")) -> totals_volcano

png(filename = file.path(parent_dir, "plots/DE_analysis", paste0("Totals_", treatment, "_volcano.png")), width = 300, height = 400)
print(totals_volcano)
dev.off()

#plot MAs----
RPFs %>%
  filter(!(is.na(padj))) %>%
  mutate(sig = factor(case_when(padj < 0.1 ~ "*",
                                padj >= 0.1 ~ "NS"))) %>%
  ggplot(aes(x = log10(baseMean), y = log2FoldChange, colour = sig))+
  geom_point(alpha = 0.5)+
  mytheme+
  xlab("log10(mean expression)")+
  ylab("log2FC")+
  ggtitle(paste(treatment, "RPFs")) -> RPFs_MA

png(filename = file.path(parent_dir, "plots/DE_analysis", paste0("RPFs_", treatment, "_MA.png")), width = 400, height = 300)
print(RPFs_MA)
dev.off()

totals %>%
  filter(!(is.na(padj))) %>%
  mutate(sig = factor(case_when(padj < 0.1 ~ "*",
                                padj >= 0.1 ~ "NS"))) %>%
  ggplot(aes(x = log10(baseMean), y = log2FoldChange, colour = sig))+
  geom_point(alpha = 0.5)+
  mytheme+
  xlab("log10(mean expression)")+
  ylab("log2FC")+
  ggtitle(paste(treatment, "Cytoplasmic RNA")) -> totals_MA

png(filename = file.path(parent_dir, "plots/DE_analysis", paste0("Totals_", treatment, "_MA.png")), width = 400, height = 300)
print(totals_MA)
dev.off()

#merge RPF with totals data----
#select apdj thresholds
TE_sig_padj <- 0.1
RPF_sig_padj <- 0.1

TE_non_sig_padj <- 0.9
RPF_non_sig_padj <- 0.5

log2FC_threshold <- 0.2

transcript_info <- 
  RPFs %>% 
  select(gene_sym, gene, transcript) %>% 
  inner_join(totals %>% select(gene_sym, gene, transcript),
             by = c("gene_sym", "gene", "transcript")) %>% 
  inner_join(TE %>% select(gene_sym, gene, transcript),
             by = c("gene_sym", "gene", "transcript"))

#merged data and make groups based on RPF/Total adjusted p-values or TE adjusted p-values
RPFs %>%
  select(gene, gene_sym, log2FoldChange, padj) %>%
  dplyr::rename(RPFs_log2FC = log2FoldChange,
         RPFs_padj = padj) %>%
  inner_join(totals[,c("gene", "log2FoldChange", "padj", "gene_sym")], by = c("gene", "gene_sym")) %>%
  dplyr::rename(totals_log2FC = log2FoldChange,
         totals_padj = padj) %>%
  inner_join(TE[,c("gene","log2FoldChange", "padj")], by = "gene") %>%
  dplyr::rename(TE_log2FC = log2FoldChange,
         TE_padj = padj) %>%
  mutate(TE_group = factor(case_when(TE_padj < TE_sig_padj & TE_log2FC < 0 ~ "TE down",
                                     TE_padj < TE_sig_padj & TE_log2FC > 0 ~ "TE up",
                                     TE_padj >= TE_non_sig_padj ~ "no change",
                                     (TE_padj >= TE_sig_padj & TE_padj <= TE_non_sig_padj) | is.na(TE_padj) ~ "NS"),
                           levels = c("TE down", "TE up", "no change", "NS"), ordered = T),
         RPFs_group = factor(case_when((RPFs_padj < RPF_sig_padj & RPFs_log2FC < -log2FC_threshold) & (totals_padj >= RPF_non_sig_padj | totals_log2FC > log2FC_threshold) ~ "RPFs down",
                                       (RPFs_padj < RPF_sig_padj & RPFs_log2FC > log2FC_threshold) & (totals_padj >= RPF_non_sig_padj | totals_log2FC < -log2FC_threshold)  ~ "RPFs up",
                                       (totals_padj < RPF_sig_padj & totals_log2FC < -log2FC_threshold) & (RPFs_padj >= RPF_non_sig_padj | RPFs_log2FC > log2FC_threshold)  ~ "Totals down",
                                       (totals_padj < RPF_sig_padj & totals_log2FC > log2FC_threshold) & (RPFs_padj >= RPF_non_sig_padj | RPFs_log2FC < -log2FC_threshold)  ~ "Totals up",
                                       RPFs_padj < RPF_sig_padj & totals_padj < RPF_sig_padj & RPFs_log2FC < -log2FC_threshold & totals_log2FC < -log2FC_threshold ~ "both down",
                                       RPFs_padj < RPF_sig_padj & totals_padj < RPF_sig_padj & RPFs_log2FC > log2FC_threshold & totals_log2FC > log2FC_threshold ~ "both up",
                                       totals_padj >= RPF_non_sig_padj & RPFs_padj >= RPF_non_sig_padj ~ "no change")),
         RPFs_group = factor(case_when(is.na(RPFs_group) ~ "NS",
                                       !(is.na(RPFs_group)) ~ RPFs_group), levels = c("RPFs down", "RPFs up", "Totals down", "Totals up", "both down", "both up", "no change", "NS"), ordered = T)) -> merged_data
summary(merged_data)





merged_data <- inner_join(merged_data, transcript_info, by = c("gene", "gene_sym"))

#write out csv
write_csv(merged_data, file.path(parent_dir, "Analysis/DESeq2_output/merged_DESeq2_siRPS25.csv"))
