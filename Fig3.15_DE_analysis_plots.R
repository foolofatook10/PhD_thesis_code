#load libraries
library(tidyverse)

#read in common variables
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/R_scripts/common_variables_siCNOT3.R")

#create a variable for what the treatment is----
treatment <- "siCNOT3"

#themes----
mytheme <- theme_classic()+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 16))

#read in DESeq2 output----
totals <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("FINAL_Totals_", treatment, "_DEseq2_apeglm_LFC_shrinkage.csv")))
RPFs <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("FINAL_RPFs_", treatment, "_DEseq2_apeglm_LFC_shrinkage.csv")))
TE <- read_csv((file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("FINAL_TE_", treatment, "_DEseq2_siCNOT3.csv"))))

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

merged_data[duplicated(merged_data$gene_sym),]
merged_data %>% filter(gene_sym== "MATR3")

merged_data %>% filter(RPFs_log2FC > 2)

merged_data %>% pull(gene_sym) %>% unique() %>% length()

1715 + 1610 + 219 + 128 + 210 +341 + 1157 + 4268

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


data_to_plot <- 
merged_data %>%
  arrange(desc(RPFs_group)) %>%
  left_join(RPF_labs, by = "RPFs_group") %>% 
  mutate(lab = factor(lab, levels = c("both down (n=1668) ",
                                      "both up (n=1587) ",
                                      "Totals down (n=368) ",
                                      "Totals up (n=331) ",
                                      "RPFs down (n=118) ",
                                      "RPFs up (n=239) ",
                                      "no change (n=1114) ",
                                      "NS (n=4034) "),
                      labels = c("both down (n=1668)",
                                 "both up (n=1587)",
                                 "Totals down (n=368)",
                                 "Totals up (n=331)",
                                 "RPFs down (n=118)",
                                 "RPFs up (n=239)",
                                 "no change (n=1114)",
                                 "non significant (n=4034)")))
  # mutate(alpha_score = case_when(is.na(RPFs_group) | RPFs_group == "no change" ~ 0.1,
  #                                RPFs_group != "NS" & RPFs_group != "no change" ~ 1)) %>%
 # #summary()

data_to_plot %>%    
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = lab, alpha = alpha_score))+
  geom_point(alpha = 0.75, size = 1.5)+
  scale_alpha(guide = "none")+
  scale_color_manual(values= my_colours) +
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

png(filename = file.path(path_to_dir, paste(treatment, "TE_scatter.png")), width = 1250, height = 1250, res = 300)
print(RPF_groups_scatter_plot)
dev.off()

library(cowplot)

leg <- cowplot::get_legend(RPF_groups_scatter_plot)

png(filename = file.path(path_to_dir, paste(treatment, "TE_scatter_legend.png")), width = 1000, height = 750, res = 300)
ggdraw(leg)
dev.off()

#based on TE
#add "n=" labels
merged_data %>%
  group_by(TE_group) %>%
  summarize(num = n()) %>%
  mutate(lab = paste0(TE_group, "\n(n=", num, ")\n")) -> TE_labs

TE_data_to_plot <- 
merged_data %>%
  filter(!(is.na(TE_group))) %>%
  # arrange(desc(TE_group)) %>%
  # mutate(alpha_score = case_when(TE_group =="TE down" | TE_group =="TE up" ~ 1,
  #                                TE_group == "no change" | is.na(TE_group) ~ 0.1)) %>%
  left_join(TE_labs, by = "TE_group") 

unique(TE_data_to_plot$lab)
  
TE_data_to_plot$lab <- factor(TE_data_to_plot$lab, 
                              levels = c("NS\n(n=8914)\n",
                                         "no change\n(n=481)\n",
                                         "TE down\n(n=1)\n",
                                         "TE up\n(n=63)\n"),
                              labels = c("Non significant\n(n=8914)\n",
                                        "No change\n(n=481)\n",
                                        "Decreased TE\n(n=1)\n",
                                        "Increased TE\n(n=63)\n"))

TE_colours <- c("#C8C8C8", "#E6AB02", "#1B9E77", "#7570B3")



TE_data_to_plot %>% 
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = lab))+
  geom_point()+
  scale_colour_manual(values= TE_colours)+
  xlab("Cytoplasmic RNA log2FC")+
  ylab("RPFs log2FC")+
  ggtitle(paste(treatment, "TE scatter"))+
  xlim(c(-2.5,2.5))+
  ylim(c(-2.5,2.5))+
  geom_abline(lty=1)+
  geom_hline(yintercept = 0, lty=1)+
  geom_hline(yintercept = 1, lty=2)+
  geom_hline(yintercept = -1, lty=2)+
  geom_vline(xintercept = 0, lty=1)+
  geom_vline(xintercept = 1, lty=2)+
  geom_vline(xintercept = -1, lty=2) +
  publication_theme() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        legend.direction = "vertical")-> TE_scatter_plot

TE_legend <- ggdraw(get_legend(TE_scatter_plot))

png(filename = file.path(path_to_dir, paste(treatment, "TE_scatter_legend2.png")), width = 750, height = 1000, res = 300)
ggdraw(TE_legend)
dev.off()

non_significant <- TE_data_to_plot %>% filter(lab == "Non significant\n(n=8914)\n")
unchanged <- TE_data_to_plot %>% filter(lab == "No change\n(n=481)\n")
TE_down <- TE_data_to_plot %>% filter(lab == "Decreased TE\n(n=1)\n")
TE_up <- TE_data_to_plot %>% filter(lab == "Increased TE\n(n=63)\n")

line_size = 0.5
line_colour = "black"
point_size = 1.5

CNOT1_label = TE_down

translation_efficiency_ggplot <- 
ggplot() +
  geom_point(data = non_significant, 
             aes(x = totals_log2FC, y = RPFs_log2FC),
             colour = "#C8C8C8",
             size = point_size) +
  geom_point(data = unchanged, 
             aes(x = totals_log2FC, y = RPFs_log2FC),
             colour = "#E6AB02",
             size = point_size) +
  geom_point(data = TE_down, 
             aes(x = totals_log2FC, y = RPFs_log2FC),
             colour = "#1B9E77",
             size = point_size) +
  geom_point(data = TE_up, 
             aes(x = totals_log2FC, y = RPFs_log2FC),
             colour = "#7570B3",
             size = point_size) +
  ggrepel::geom_text_repel(data = CNOT1_label,
                           aes(x = totals_log2FC, y = RPFs_log2FC, label = gene_sym),
                           colour = "#1B9E77",
                           nudge_x = 1,
                           nudge_y = -1,
                           size = 5,
                           fontface = "bold") +
  xlab("Cytoplasmic RNA log2FC")+
  ylab("RPFs log2FC") +
  geom_hline(yintercept = 0, lty=1, colour = line_colour, size = line_size)+
  geom_hline(yintercept = 1, lty=2, colour = line_colour, size = line_size)+
  geom_hline(yintercept = -1, lty=2, colour = line_colour, size = line_size)+
  geom_vline(xintercept = 0, lty=1, colour = line_colour, size = line_size)+
  geom_vline(xintercept = 1, lty=2, colour = line_colour, size = line_size)+
  geom_vline(xintercept = -1, lty=2, colour = line_colour, size = line_size) +
  publication_theme()

png(filename = file.path(path_to_dir, paste(treatment, "TE_scatter_TE.png")), width = 1250, height = 1250, res = 300)
print(translation_efficiency_ggplot)
dev.off()


CCR4_NOT_components <- 
TE_data_to_plot %>%
  mutate(CCR4_NOT = factor(case_when(gene_sym == "CNOT1" ~ "CNOT1",
                              gene_sym == "CNOT2" ~ "CNOT2",
                              gene_sym == "CNOT3" ~ "CNOT3",
                              gene_sym == "CNOT4" ~ "CNOT4",
                              gene_sym == "CNOT6" ~ "CNOT6",
                              gene_sym == "CNOT7" ~ "CNOT7",
                              gene_sym == "CNOT8" ~ "CNOT8",
                              gene_sym == "CNOT9" ~ "CNOT9",
                              gene_sym == "CNOT11" ~"CNOT11",
                              .default = "other"),
                           levels = c(paste0("CNOT", 1:4),
                                      paste0("CNOT", 6:9),
                                      "CNOT11",
                                      "other"))) 

CNOTs <- CCR4_NOT_components %>% filter(CCR4_NOT != "other")

other <- CCR4_NOT_components %>% filter(CCR4_NOT == "other")
CNOT1 <- CCR4_NOT_components %>% filter(CCR4_NOT == "CNOT1")
CNOT2 <- CCR4_NOT_components %>% filter(CCR4_NOT == "CNOT2")
CNOT3 <- CCR4_NOT_components %>% filter(CCR4_NOT == "CNOT3")
CNOT4 <- CCR4_NOT_components %>% filter(CCR4_NOT == "CNOT4")
CNOT6 <- CCR4_NOT_components %>% filter(CCR4_NOT == "CNOT6")
CNOT7 <- CCR4_NOT_components %>% filter(CCR4_NOT == "CNOT7")
CNOT8 <- CCR4_NOT_components %>% filter(CCR4_NOT == "CNOT8")
CNOT9 <- CCR4_NOT_components %>% filter(CCR4_NOT == "CNOT9")
CNOT11 <- CCR4_NOT_components %>% filter(CCR4_NOT == "CNOT11")


# ggplot() +
#   geom_point(data = other, 
#              aes(x = totals_log2FC, y = RPFs_log2FC),
#              colour = "#A3A3A3") +
#   geom_point(data = CNOT1,
#              aes(x = totals_log2FC, y = RPFs_log2FC),
#              colour = "#E41A1C")+
#   geom_point(data = CNOT2,
#              aes(x = totals_log2FC, y = RPFs_log2FC),
#              colour = "#377EB8")+
#   geom_point(data = CNOT3,
#              aes(x = totals_log2FC, y = RPFs_log2FC),
#              colour = "#4DAF4A")+
#   geom_point(data = CNOT4,
#              aes(x = totals_log2FC, y = RPFs_log2FC),
#              colour = "#984EA3")+
#   geom_point(data = CNOT6,
#              aes(x = totals_log2FC, y = RPFs_log2FC),
#              colour = "#FF7F00")+
#   geom_point(data = CNOT7,
#              aes(x = totals_log2FC, y = RPFs_log2FC),
#              colour = "#A65628")+
#   geom_point(data = CNOT8,
#              aes(x = totals_log2FC, y = RPFs_log2FC),
#              colour = "#F781BF") +
#   geom_point(data = CNOT9,
#              aes(x = totals_log2FC, y = RPFs_log2FC),
#              colour = "#A6D854")+
#   geom_point(data = CNOT11,
#              aes(x = totals_log2FC, y = RPFs_log2FC),
#              colour = "#FC8D62") +
#   ggrepel::geom_text_repel(data = CNOT1,
#              aes(x = totals_log2FC, y = RPFs_log2FC, label = gene_sym),
#              colour = "#E41A1C")+
#   ggrepel::geom_text_repel(data = CNOT2,
#              aes(x = totals_log2FC, y = RPFs_log2FC, label = gene_sym),
#              colour = "#377EB8")+
#   ggrepel::geom_text_repel(data = CNOT3,
#              aes(x = totals_log2FC, y = RPFs_log2FC, label = gene_sym),
#              colour = "#4DAF4A")+
#   ggrepel::geom_text_repel(data = CNOT4,
#              aes(x = totals_log2FC, y = RPFs_log2FC, label = gene_sym),
#              colour = "#984EA3")+
#   ggrepel::geom_text_repel(data = CNOT6,
#              aes(x = totals_log2FC, y = RPFs_log2FC, label = gene_sym),
#              colour = "#FF7F00")+
#   ggrepel::geom_text_repel(data = CNOT7,
#              aes(x = totals_log2FC, y = RPFs_log2FC, label = gene_sym),
#              colour = "#A65628")+
#   ggrepel::geom_text_repel(data = CNOT8,
#              aes(x = totals_log2FC, y = RPFs_log2FC, label = gene_sym),
#              colour = "#F781BF") +
#   ggrepel::geom_text_repel(data = CNOT9,
#              aes(x = totals_log2FC, y = RPFs_log2FC, label = gene_sym),
#              colour = "#A6D854")+
#   ggrepel::geom_text_repel(data = CNOT11,
#              aes(x = totals_log2FC, y = RPFs_log2FC, label = gene_sym),
#              colour = "#FC8D62")

CCR4_NOT_components_colours <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#00A087FF", "#A65628", "#3C5488FF", "#A6D854", "#E64B35B2")


CCR4_NOTs <- 
ggplot() +
  geom_point(data = other, 
             aes(x = totals_log2FC, y = RPFs_log2FC),
             colour = "#C8C8C8") +
  geom_point(data = CNOTs,
             aes(x = totals_log2FC, y = RPFs_log2FC, colour = CCR4_NOT),
             show.legend = F, size = 1.5) +
  ggrepel::geom_text_repel(data = CNOTs,
                           aes(x = totals_log2FC, y = RPFs_log2FC, colour = CCR4_NOT, label = gene_sym),
                           show.legend = F,
                           fontface = "bold",
                           size = 5) +
  scale_color_manual(values = CCR4_NOT_components_colours) +
  xlab("Cytoplasmic RNA log2FC")+
  ylab("RPFs log2FC") +
  ggtitle("CCR4-NOT components") +
  publication_theme()
  

png(filename = file.path(path_to_dir, paste(treatment, "TE_scatter_CCR4_NOTs.png")), width = 1250, height = 1250, res = 300)
print(CCR4_NOTs)
dev.off()


#"#66C2A5"  "#8DA0CB" "#E78AC3"  "#FFD92F" "#E5C494" "#B3B3B3"

  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = CCR4_NOT)) +
  geom_point()
  filter(str_detect(gene_sym, "CNOT"))
  

png(filename = file.path(parent_dir, "plots/DE_analysis", paste(treatment, "_TE_scatter.png")), width = 500, height = 400)
print(TE_scatter_plot)
dev.off()

#write out csv
write_csv(merged_data, file.path(parent_dir, "Analysis/DESeq2_output", paste0("FINAL_", treatment, "_merged_DESeq2.csv")))

