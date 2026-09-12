#load libraries
library(tidyverse)

#read in common variables
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/R_scripts/common_variables_siRPS25.R")

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
totals <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("FINAL_Totals_", treatment, "_DEseq2_apeglm_LFC_shrinkage_contains_eS25.csv")))
RPFs <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("FINAL_RPFs_", treatment, "_DEseq2_apeglm_LFC_shrinkage_contains_eS25.csv")))
TE <- read_csv((file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("FINAL_TE_", treatment, "_DEseq2_contains_eS25.csv"))))

# #plot volcanos----
# RPFs %>%
#   filter(!(is.na(padj))) %>%
#   mutate(sig = factor(case_when(padj < 0.1 ~ "*",
#                                 padj >= 0.1 ~ "NS"))) %>%
#   ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = sig))+
#   geom_point(alpha = 0.5)+
#   mytheme+
#   xlab("log2FC")+
#   ylab("-log10(padj)")+
#   ggtitle(paste(treatment, "RPFs")) -> RPFs_volcano
# 
# png(filename = file.path(parent_dir, "plots/DE_analysis", paste0("RPFs_", treatment, "_volcano.png")), width = 300, height = 400)
# print(RPFs_volcano)
# dev.off()
# 
# totals %>%
#   filter(!(is.na(padj))) %>%
#   mutate(sig = factor(case_when(padj < 0.1 ~ "*",
#                                 padj >= 0.1 ~ "NS"))) %>%
#   ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = sig))+
#   geom_point(alpha = 0.5)+
#   mytheme+
#   xlab("log2FC")+
#   ylab("-log10(padj)")+
#   ggtitle(paste(treatment, "Cytoplasmic RNA")) -> totals_volcano
# 
# png(filename = file.path(parent_dir, "plots/DE_analysis", paste0("Totals_", treatment, "_volcano.png")), width = 300, height = 400)
# print(totals_volcano)
# dev.off()
# 
# #plot MAs----
# RPFs %>%
#   filter(!(is.na(padj))) %>%
#   mutate(sig = factor(case_when(padj < 0.1 ~ "*",
#                                 padj >= 0.1 ~ "NS"))) %>%
#   ggplot(aes(x = log10(baseMean), y = log2FoldChange, colour = sig))+
#   geom_point(alpha = 0.5)+
#   mytheme+
#   xlab("log10(mean expression)")+
#   ylab("log2FC")+
#   ggtitle(paste(treatment, "RPFs")) -> RPFs_MA
# 
# png(filename = file.path(parent_dir, "plots/DE_analysis", paste0("RPFs_", treatment, "_MA.png")), width = 400, height = 300)
# print(RPFs_MA)
# dev.off()
# 
# totals %>%
#   filter(!(is.na(padj))) %>%
#   mutate(sig = factor(case_when(padj < 0.1 ~ "*",
#                                 padj >= 0.1 ~ "NS"))) %>%
#   ggplot(aes(x = log10(baseMean), y = log2FoldChange, colour = sig))+
#   geom_point(alpha = 0.5)+
#   mytheme+
#   xlab("log10(mean expression)")+
#   ylab("log2FC")+
#   ggtitle(paste(treatment, "Cytoplasmic RNA")) -> totals_MA
# 
# png(filename = file.path(parent_dir, "plots/DE_analysis", paste0("Totals_", treatment, "_MA.png")), width = 400, height = 300)
# print(totals_MA)
# dev.off()

#merge RPF with totals data----
#select apdj thresholds
TE_sig_padj <- 0.1
RPF_sig_padj <- 0.1

TE_non_sig_padj <- 0.9
RPF_non_sig_padj <- 0.5

log2FC_threshold <- 0.1

RPFs %>%
  dplyr::select(gene, gene_sym, log2FoldChange, padj) %>%
  dplyr::rename(RPFs_log2FC = log2FoldChange,
                RPFs_padj = padj) %>%
  dplyr::inner_join(totals[,c("gene", "log2FoldChange", "padj", "gene_sym")], by = c("gene", "gene_sym")) %>%
  dplyr::rename(totals_log2FC = log2FoldChange,
                totals_padj = padj) %>%
  dplyr::inner_join(TE[,c("gene","log2FoldChange", "padj")], by = "gene") %>%
  dplyr::rename(TE_log2FC = log2FoldChange,
                TE_padj = padj) -> test
# 
# test %>% 
#   top_n(341, RPFs_log2FC) %>% 
#   pull(RPFs_padj) %>% min()

#merged data and make groups based on RPF/Total adjusted p-values or TE adjusted p-values
RPFs %>%
  dplyr::select(gene, gene_sym, log2FoldChange, padj) %>%
  dplyr::rename(RPFs_log2FC = log2FoldChange,
                RPFs_padj = padj) %>%
  dplyr::inner_join(totals[,c("gene", "log2FoldChange", "padj", "gene_sym")], by = c("gene", "gene_sym")) %>%
  dplyr::rename(totals_log2FC = log2FoldChange,
                totals_padj = padj) %>%
  dplyr::inner_join(TE[,c("gene","log2FoldChange", "padj")], by = "gene") %>%
  dplyr::rename(TE_log2FC = log2FoldChange,
                TE_padj = padj) %>%
  dplyr::mutate(TE_group = factor(case_when(TE_padj < TE_sig_padj & TE_log2FC < 0 ~ "TE down",
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
                                              totals_padj >= RPF_non_sig_padj & RPFs_padj >= RPF_non_sig_padj ~ "no change",
                                              .default = "NS"), levels = c("RPFs down", "RPFs up", "Totals down", "Totals up", "both down", "both up", "no change", "NS"), ordered = T)) -> merged_data
summary(merged_data$RPFs_group)

nrow(merged_data)


#plot TE scatters----
#based on RPFs/totals logFC
#add "n=" labels
merged_data %>%
  group_by(RPFs_group) %>%
  summarize(num = n()) %>%
  mutate(lab = paste0(RPFs_group, " (n=", num, ") ")) -> RPF_labs

# read in publication theme
publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            #axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2, size = 20),
            axis.title.x = element_text(vjust = -0.2, size = 20),
            axis.text.x = element_text(size = 16, color = "black"),
            axis.text.y = element_text(size = 16, color = "black"),
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
            legend.margin = unit(1, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

#RColorBrewer::display.brewer.all()
#RColorBrewer::brewer.pal(10, "Paired")
#"#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "" "#CAB2D6" "#6A3D9A"

my_colours <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", "#C8C8C8")

RPF_labs %>% filter(RPFs_group == "RPFs down") %>% pull(lab)

data_to_plot <- 
  merged_data %>%
  arrange(desc(RPFs_group)) %>%
  left_join(RPF_labs, by = "RPFs_group") %>% 
  mutate(lab = factor(lab, levels = c(RPF_labs %>% filter(RPFs_group == "both down") %>% pull(lab),
                                      RPF_labs %>% filter(RPFs_group == "both up") %>% pull(lab),
                                      RPF_labs %>% filter(RPFs_group == "Totals down") %>% pull(lab),
                                      RPF_labs %>% filter(RPFs_group == "Totals up") %>% pull(lab),
                                      RPF_labs %>% filter(RPFs_group == "RPFs down") %>% pull(lab),
                                      RPF_labs %>% filter(RPFs_group == "RPFs up") %>% pull(lab),
                                      RPF_labs %>% filter(RPFs_group == "no change") %>% pull(lab),
                                      RPF_labs %>% filter(RPFs_group == "NS") %>% pull(lab)),
                      labels = c(str_remove(RPF_labs %>% filter(RPFs_group == "both down") %>% pull(lab), pattern = "\\s$"),
                                 str_remove(RPF_labs %>% filter(RPFs_group == "both up") %>% pull(lab), pattern = "\\s$"),
                                 str_remove(RPF_labs %>% filter(RPFs_group == "Totals down") %>% pull(lab), pattern = "\\s$"),
                                 str_remove(RPF_labs %>% filter(RPFs_group == "Totals up") %>% pull(lab), pattern = "\\s$"),
                                 str_remove(RPF_labs %>% filter(RPFs_group == "RPFs down") %>% pull(lab), pattern = "\\s$"),
                                 str_remove(RPF_labs %>% filter(RPFs_group == "RPFs up") %>% pull(lab), pattern = "\\s$"),
                                 str_remove(RPF_labs %>% filter(RPFs_group == "no change") %>% pull(lab), pattern = "\\s$"),
                                 str_remove(RPF_labs %>% filter(RPFs_group == "NS") %>% pull(lab), pattern = "\\s$"))))
# mutate(alpha_score = case_when(is.na(RPFs_group) | RPFs_group == "no change" ~ 0.1,
#                                RPFs_group != "NS" & RPFs_group != "no change" ~ 1)) %>%
# #summary()

es25 <- data_to_plot %>% filter(gene_sym == "RPS25") %>% mutate(gene_sym = "eS25")

   
ggplot()+
  geom_point(data = data_to_plot, 
             aes(x = totals_log2FC, y = RPFs_log2FC),
             alpha = 0.75, size = 1.5, colour = "#C8C8C8") +
  geom_point(data = es25, 
             aes(x = totals_log2FC, y = RPFs_log2FC),
               alpha = 0.75, size = 1.5, colour = "black")+
  ggrepel::geom_text_repel(data = es25,
                           aes(x = totals_log2FC, y = RPFs_log2FC, label = gene_sym),
                           vjust = -1,
                           hjust = -0.25,
                           size = 6,
                           colour = "black") +
  xlab("Cytoplasmic RNA log2FC")+
  ylab("RPFs log2FC")+
  ggtitle(paste(treatment, "\nTranslation efficiency"))+
  #xlim(c(-2.5,2.5))+
  #ylim(c(-2.5,2.5))+
  geom_abline(linetype = "dashed", size = 1, colour = "black")+
  geom_hline(yintercept = 0, linetype = "dashed", size = 1, colour = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", size = 1, colour = "black") +
  publication_theme() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank()) -> RPF_groups_scatter_plot

path_to_dir <- "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures"

png(filename = file.path(path_to_dir, paste(treatment, "TE_scatter_with_eS25.png")), width = 1250, height = 1250, res = 300)
print(RPF_groups_scatter_plot)
dev.off()



