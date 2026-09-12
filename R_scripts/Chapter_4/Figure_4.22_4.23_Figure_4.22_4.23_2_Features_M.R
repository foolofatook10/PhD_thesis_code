#### Load in libraries ####
library(data.table)
library(tidyverse)
library(seqinr)
#library(reshape2) 
#library(car)
library(ggpubr)
library(RColorBrewer)
library(purrr)
library(rstatix)
library(gridExtra)
#library(heatmaply)

library(scales)
library(fgsea)

#### Set paths and directories ####

#Set directories
parent_dir = "\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL"

#Make two subdirectories, plots and my_plots
setwd(file.path(parent_dir, "plots/"))

path2deseq = file.path(parent_dir,"Analysis/DESeq2_output/RPFs_CNOT3_DEseq2_apeglm_LFC_shrinkage_M.csv")

#### Publication theme ####

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
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
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
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

#### Functions ####
aov_wrapper <- function(data, formula){
  aov(formula = formula, data = data)
}

anova_fun <- function(data, feature, log = F){
  
  if(log == T){
    
    anova <- 
      data %>%
      mutate(log_feature = log10(eval(parse(text = feature))) + 0.0001) %>%
      aov_wrapper(formula = as.formula(log_feature ~ group))
    
    pval <- ifelse(formatC(summary(anova)[[1]]$'Pr(>F)'[1], format = "f", digits = 3) < 0.001, 
                   yes = "<0.001",
                   no =  formatC(summary(anova)[[1]]$'Pr(>F)'[1], format = "f", digits = 3))
    
    return(pval)
    
  } else
    
    anova <- 
      data %>%
      aov_wrapper(formula = as.formula(eval(parse(text = feature)) ~ group))
  
  pval <- ifelse(formatC(summary(anova)[[1]]$'Pr(>F)'[1], format = "f", digits = 3) < 0.001, 
                 yes = "<0.001",
                 no =  formatC(summary(anova)[[1]]$'Pr(>F)'[1], format = "f", digits = 3))
  
  return(pval)
  
}

TUKEY_to_df <- function(data, feature, log = F){
  if(log == TRUE) {
    
    anova <- 
      data %>%
      mutate(log_feature = log10(eval(parse(text = feature))) + 0.0001) %>%
      aov_wrapper(formula = as.formula(log_feature ~ group))
    
    tukey_dataframe <- data.frame((TukeyHSD(anova))$"group") %>%
      tibble::rownames_to_column("group1") %>%
      mutate(p.adj = case_when(p.adj < 0.001 ~ "<0.001",
                               .default = as.character(formatC(p.adj, format = "f", digits = 3))))
    
    tukey_dataframe[c("group1", "group2")] <- stringr::str_split_fixed(tukey_dataframe$group1, '-', 2)
    
    padjs <- as_tibble(tukey_dataframe %>% select(group1,group2,p.adj, diff:upr))
    padjs <- padjs[c(1,3,2),]
    padjs <- padjs %>% select(group1, group2, p.adj)
    padjs$p.adj <- formatC(padjs$p.adj, format = "f", digits = 3)
    
    return(padjs)
    
  } else {
    
    anova <- 
      data %>%
      aov_wrapper(formula = as.formula(eval(parse(text = feature)) ~ group))
    
    tukey_dataframe <- data.frame((TukeyHSD(anova))$"group") %>%
      tibble::rownames_to_column("group1") %>%
      mutate(p.adj = case_when(p.adj < 0.001 ~ "<0.001",
                               .default = as.character(formatC(p.adj, format = "f", digits = 3))))
    tukey_dataframe[c("group1", "group2")] <- stringr::str_split_fixed(tukey_dataframe$group1, '-', 2)
    
    padjs <- as_tibble(tukey_dataframe %>% select(group1,group2,p.adj, diff:upr))
    padjs <- padjs[c(1,3,2),]
    padjs <- padjs %>% select(group1, group2, p.adj)
    padjs$p.adj <- formatC(padjs$p.adj, format = "f", digits = 3)
    
    return(padjs)
  }
}


kruskal_wrapper <-  function(data, formula){
  kruskal.test(formula = formula, data = data)
}


kruskal_fun <- function(data, feature){
  
  kruskal <- 
    data %>%
    kruskal_wrapper(formula = as.formula(eval(parse(text = feature)) ~ group))
  
  pval <- kruskal$p.value
  
  pval <- as.numeric(formatC(pval, format = "e", digits = 2))
  
  pval <- ifelse(pval < 0.001, 
                 yes = "<0.001",
                 no =  pval)
  
  return(pval)
  
}

kruskal_post_hoc <- function(data, feature){
  
  # Create the wilcox test object
  wilcox_pairwise <- 
    pairwise.wilcox.test(data[[feature]], 
                         data[["group"]],
                         p.adjust.method = "BH", 
                         detailed = T, 
                         conf.level = 0.95)
  
  # Extract the p adjusted values into a tibble
  padjs <- rownames_to_column(as_tibble(wilcox_pairwise$p.value)) %>%
    gather(key = "group2", value = "p.adj") %>%
    filter(group2 == "downregulated" | group2 == "unchanged") %>%
    na.omit() %>%
    mutate(group1 = c("unchanged", "upregulated", "upregulated")) %>%
    select(group1, group2, p.adj)
  
  # Make the p adjusted values numerical not character strings 
  padjs$p.adj<- as.double(padjs$p.adj)
  
  # Reorder the p.adjusted values into the correct order to put on a ggplot
  padjs <- padjs[c(1,3,2),]
  
  # Round the p adjusted values to 2 significant figures
  padjs$p.adj <- as.numeric(formatC(padjs$p.adj, format = "e", digits = 2))
  
  # Replcae values with corresponding labels
  padjs <- padjs %>% 
    mutate(p.adj = case_when(padjs$p.adj < 0.001 ~ "<0.001", 
                             .default = as.character(formatC(padjs$p.adj, format = "f", digits = 3))))
  
  #return p.adjusted values
  return(padjs)
  
}

get_stats <- function(data, feature) {
  stats_for_feature <- data %>%
    group_by(group) %>%
    summarise(mean = mean(eval(parse(text = feature))),
              median = median(eval(parse(text = feature))),
              n = length(eval(parse(text = feature))),
              sd = sd(eval(parse(text = feature))),
              se = sd/sqrt(n),
              t.score = qt(p=0.05/2, df=n-1,lower.tail=F),
              margin_error = se * t.score,
              CI_lower = mean - margin_error,
              CI_upper = mean + margin_error)
  
  return(stats_for_feature)
}

#### Load in tables ####
master <- fread(file = "\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/master.csv", header = T, drop = "V1")
DEseq2 <- read_csv(file = path2deseq)
features <- fread("\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/features.csv", header = T, drop = "V1")


#Import codons table
path2 = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/useful_tables"
codons <- read_csv(paste0(path2, "/codon_box_types.csv"))
codons$codon <- factor(tolower(codons$codon))
codons$codon <- gsub("u","t", codons$codon)
code_ons <- as.character(codons$codon)
code_ons_plus_stop <- append(code_ons, c("taa", "tag", "tga"))

#Define 3 groups
DEseq2 <- DEseq2 %>% 
  mutate(group = factor(case_when(log2FoldChange <= 0 & padj < 0.05 ~ "downregulated",
                                  log2FoldChange >= 0 & padj < 0.05 ~ "upregulated",
                                  TRUE ~ "unchanged")))

#Extract significant genes only
sig_genes_downregulated <- DEseq2 %>% filter(group == "downregulated")
sig_genes_upregulated <- DEseq2 %>% filter(group == "upregulated")
unchanged_genes <- DEseq2 %>% filter(group == "unchanged")

#### volvano plot ####

#Get top 10 significantly DE genes in terms of p.adj
siggenes_down <- rbind(sig_genes_downregulated %>% top_n(-10, padj), 
                       sig_genes_downregulated %>% 
                         filter(log2FoldChange < -2)) %>% 
  unique()
siggenes_up <- right_join(sig_genes_upregulated %>% 
                       top_n(-10, padj), 
                     sig_genes_upregulated %>% 
                       filter(log2FoldChange > 2.75)) %>%
                       filter(gene_sym != "CNOT3") %>% 
                       unique()


siggenes_up <-  sig_genes_upregulated %>% 
  filter(log2FoldChange > 2.6) %>%
  filter(gene_sym != "CNOT3") %>% 
  unique()

#siggenes_up <-siggenes_up %>% filter(gene_sym != c("GABARAP"))

CNOT3 <- DEseq2 %>%
  filter(gene_sym == "CNOT3")

brewer.pal(8, "Set2")
#plot ggplot
gplot1 <- ggplot(DEseq2, aes(x = log2FoldChange,  y = -log10(padj), colour = group)) + 
  geom_point() +
  scale_colour_manual(values = c("#1B9E77", "#D95F02", "#7570B3"),
                      name = "Group",
                      labels = c(paste0("Depleted \n n = ", nrow(sig_genes_downregulated)),
                                 paste0("Unchanged \n n = ", nrow(unchanged_genes)),
                                 paste0("Enriched \n n = ", nrow(sig_genes_upregulated)))) +
  # ggrepel::geom_text_repel(data = siggenes_down, 
  #                          aes(x = log2FoldChange,  y = -log10(padj), label = gene_sym),
  #                          size = 2.5,
  #                          fontface = "bold", 
  #                          show.legend = F) +
  ggrepel::geom_text_repel(data = siggenes_up, 
                           aes(x = log2FoldChange,  y = -log10(padj), label = gene_sym),
                           size = 3,
                           colour = "black",
                           fontface = "bold",
                           show.legend = F) +
  geom_point(data = CNOT3, aes(x = log2FoldChange,  y = -log10(padj)), colour = "black") +
  ggrepel::geom_text_repel(data= CNOT3, 
                           aes(x = log2FoldChange,  y = -log10(padj), label = gene_sym), 
                           min.segment.length = unit(0, 'lines'), 
                           size = 3.25,
                           fontface = "bold",
                           nudge_y = 8, 
                           nudge_x = -1, 
                           color = "black",
                           show.legend = F) +
  #ggtitle("Monosomes") +
  xlab("Log2 Fold Enrichment (CNOT3 IP - Total RPFs)") +
  ylab("-Log10(padj)") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 10, hjust = 0.5),
        axis.title = element_text(size = 2)) +
  publication_theme() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

print(gplot1)

png(filename = "volcano.png",
     height = 1500,
     width = 1500,
    res =300)
print(gplot1)
dev.off()

#### Join features ####
seq_info <- 
  inner_join(features, 
             DEseq2, 
             by = c(
               "ENST" = "transcript"))

seq_info <- na.omit(seq_info)


stats_lyst <- list()
#### Lengths ####

##### 5'UTR #####
stats <- get_stats(seq_info, "fputr_length")

stats_lyst[["5'UTR_length"]] <- stats

p_val <- anova_fun(seq_info, feature = "fputr_length", log = T)
padjs<- TUKEY_to_df(seq_info, feature = "fputr_length", log = T)

gplot_5UTR_length <- 
  ggplot(seq_info, aes(x = group, y = fputr_length)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Dep."), # \nn=", sum(sequence_info_lengths$group == "downregulated")/3),
                              paste0("Unch."), # \nn=", sum(sequence_info_lengths$group == "unchanged")/3),
                              paste0("Enr."))) + # \nn=", sum(sequence_info_lengths$group == "upregulated")/3))) +
  scale_y_continuous(trans = c("log10")) +
  ggtitle("5'UTR Length") +
  ylab("Nucleotides") +
  xlab("Group") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

ylim <- layer_scales(gplot_5UTR_length)$y$range$range[2] * 1.05
gplot_5UTR_length <- gplot_5UTR_length + 
  stat_pvalue_manual(padjs, 
                     y.position = ylim, 
                     step.increase = 0.06, 
                     label = "p.adj", 
                     tip.length = 0.01, 
                     fontface = "bold",
                     bracket.size = 1.2)
#ylim2 <- 10^(layer_scales(gplot_CDS_lengt 
#ylim2 <- 10^(layer_scales(gplot_5UTR_length)$y$range$range[2]) * 2
#gplot_5UTR_length <- gplot_5UTR_length + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold")

png(filename = "5'UTR_lengths.png",
     width = 1200, 
     height = 1500,
     res = 300)
print(gplot_5UTR_length)
dev.off()

##### CDS Length #####
stats <- get_stats(seq_info, "cds_length")

stats_lyst[["CDS_length"]] <- stats

p_val <- anova_fun(seq_info, feature = "cds_length", log = T)
padjs <- TUKEY_to_df(seq_info, feature = "cds_length", log = T)

gplot_CDS_length <- 
  ggplot(seq_info, aes(x = group, y = cds_length)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Dep."), # \nn=", sum(sequence_info_lengths$group == "downregulated")/3),
                              paste0("Unch."), # \nn=", sum(sequence_info_lengths$group == "unchanged")/3),
                              paste0("Enr."))) + # \nn=", sum(sequence_info_lengths$group == "upregulated")/3))) +
  scale_y_continuous(trans = c("log10"), labels = label_comma()) +
  ylab("CDS Length (nucleotides)") +
  xlab("Group") +
  ggtitle("CDS Length") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank())
        ##axis.text.x = element_text(face = "bold", size = 11))

ylim <- layer_scales(gplot_CDS_length)$y$range$range[2] * 1.05
gplot_CDS_length <- gplot_CDS_length + 
  stat_pvalue_manual(padjs, 
                     y.position = ylim, 
                     step.increase = 0.06, 
                     label = "p.adj", 
                     tip.length = 0.01, 
                     fontface = "bold",
                     bracket.size = 1.2)
#ylim2 <- 10^(layer_scales(gplot_CDS_length)$y$range$range[2]) * 1.5
#gplot_CDS_length <- gplot_CDS_length + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold")

png(filename = "CDS_lengths.png",
    width = 1200, 
    height = 1500,
    res = 300)
print(gplot_CDS_length)
dev.off()

##### 3'UTR Length #####
stats <- get_stats(seq_info, "tputr_length")

stats_lyst[["3'UTR_length"]] <- stats

p_val <- anova_fun(seq_info, feature = "tputr_length", log = T)
padjs <- TUKEY_to_df(seq_info, feature = "tputr_length", log = T)

gplot_3UTR_length <- 
  ggplot(seq_info, aes(x = group, y = tputr_length)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Dep."), # \nn=", sum(sequence_info_lengths$group == "downregulated")/3),
                              paste0("Unch."), # \nn=", sum(sequence_info_lengths$group == "unchanged")/3),
                              paste0("Enr."))) + # \nn=", sum(sequence_info_lengths$group == "upregulated")/3))) +
  scale_y_continuous(trans = c("log10"), labels = label_comma()) +
  ylab("3'UTR Length (nucleotides)") +
  xlab("Group") +
  ggtitle("3'UTR Length") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
        #axis.text.x = element_text(face = "bold", size = 10))

ylim <- layer_scales(gplot_3UTR_length)$y$range$range[2] * 1.05
gplot_3UTR_length <- gplot_3UTR_length + 
  stat_pvalue_manual(padjs, 
                     y.position = ylim, 
                     step.increase = 0.06, 
                     label = "p.adj", 
                     tip.length = 0.01, 
                     fontface = "bold",
                     bracket.size = 1.2)
#ylim2 <- 10^(layer_scales(gplot_3UTR_length)$y$range$range[2]) * 2
#gplot_3UTR_length <- gplot_3UTR_length + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold")

png(filename = "3'UTR_lengths.png",
    width = 1200, 
    height = 1500,
    res = 300)
print(gplot_3UTR_length)
dev.off()

combined_figure <- 
  ggarrange(gplot_5UTR_length, gplot_CDS_length, gplot_3UTR_length, nrow = 1)

final_figure <- annotate_figure(combined_figure, 
                left = textGrob("Nucleotides", 
                rot = 90, vjust = 1,
                gp = gpar(cex = 2.5)))

tiff(filename = "lengths.tiff",
     width = 3000, 
     height = 1000,
     res = 200)
print(final_figure)
dev.off()

#### GC content ####

##### 5'UTR #####
stats <- get_stats(seq_info, feature = "GC_cont_fp")

stats_lyst[["5'UTR_GC"]] <- stats

p_val <- anova_fun(seq_info, feature = "GC_cont_fp")
padjs <- TUKEY_to_df(seq_info, feature = "GC_cont_fp")

gplot_5UTR_GC <- 
  ggplot(seq_info, aes(x = group, y = GC_cont_fp * 100)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean* 100), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Dep."), # \nn=", sum(sequence_info_lengths$group == "downregulated")/3),
                              paste0("Unch."), # \nn=", sum(sequence_info_lengths$group == "unchanged")/3),
                              paste0("Enr."))) + # \nn=", sum(sequence_info_lengths$group == "upregulated")/3))) +
  #scale_y_continuous(trans = c("log10"), labels = label_comma()) +
  ylab("GC%") +
  xlab("Group") +
  ggtitle("5'UTR") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.x = element_blank())

ylim <- layer_scales(gplot_5UTR_GC)$y$range$range[2] * 1.05
gplot_5UTR_GC <- gplot_5UTR_GC + 
  stat_pvalue_manual(padjs, 
                     y.position = ylim, 
                     step.increase = 0.08, 
                     label = "p.adj", 
                     tip.length = 0.01, 
                     fontface = "bold", 
                     size = 4,
                     bracket.size = 1.2)
#ylim2 <- (layer_scales(gplot_5UTR_GC)$y$range$range[2]) * 1.055
#gplot_5UTR_GC <- gplot_5UTR_GC + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold", size = 4.5)

png(filename = "5'UTR_GC.png",
     res = 300,
     width = 1200, 
     height = 1500)
print(gplot_5UTR_GC)
dev.off()

#### CDS #####
stats <- get_stats(seq_info, feature = "GC_cont_cds")

stats_lyst[["CDS_GC"]] <- stats

p_val <- anova_fun(seq_info, feature = "GC_cont_cds")
padjs <- TUKEY_to_df(seq_info, feature = "GC_cont_cds")

gplot_CDS_GC <- 
  ggplot(seq_info, aes(x = group, y = GC_cont_cds * 100)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean * 100), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Dep."), # \nn=", sum(sequence_info_lengths$group == "downregulated")/3),
                              paste0("Unch."), # \nn=", sum(sequence_info_lengths$group == "unchanged")/3),
                              paste0("Enr."))) + # \nn=", sum(sequence_info_lengths$group == "upregulated")/3))) +
  #scale_y_continuous(trans = c("log10"), labels = label_comma()) +
  ggtitle("CDS") +
  ylab("GC%") +
  xlab("Group") +
  #scale_y_continuous(limits = c(0,NA)) +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.x = element_blank())

ylim <- layer_scales(gplot_CDS_GC)$y$range$range[2] * 1.02
gplot_CDS_GC <- gplot_CDS_GC + stat_pvalue_manual(padjs, y.position = ylim, 
                                                  step.increase = 0.075, 
                                                  bracket.size = 1.2,
                                                  label = "p.adj", tip.length = 0.01, 
                                                  fontface = "bold", size = 4.5)
ylim2 <- (layer_scales(gplot_CDS_GC)$y$range$range[2]) * 1.04
#gplot_CDS_GC <- gplot_CDS_GC + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold", size = 4.5)

png(filename = "CDS_GC.png",
    res = 300,
    width = 1200, 
    height = 1500)
print(gplot_CDS_GC)
dev.off()

##### 3'UTR #####
stats <- get_stats(seq_info, feature = "GC_cont_tp")

stats_lyst[["3'UTR_GC"]] <- stats

p_val <- anova_fun(seq_info, feature = "GC_cont_tp")
padjs <- TUKEY_to_df(seq_info, feature = "GC_cont_tp")

gplot_3UTR_GC <- 
  ggplot(seq_info, aes(x = group, y = GC_cont_tp * 100)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean * 100), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Dep."), # \nn=", sum(sequence_info_lengths$group == "downregulated")/3),
                              paste0("Unch."), # \nn=", sum(sequence_info_lengths$group == "unchanged")/3),
                              paste0("Enr."))) + # \nn=", sum(sequence_info_lengths$group == "upregulated")/3))) +
  #scale_y_continuous(trans = c("log10"), labels = label_comma()) +
  ylab("GC%") +
  xlab("Group") +
  ggtitle("3'UTR") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.x = element_blank())

ylim <- layer_scales(gplot_3UTR_GC)$y$range$range[2] * 1.025
gplot_3UTR_GC <- 
  gplot_3UTR_GC + 
  stat_pvalue_manual(padjs, 
                     y.position = ylim, 
                     step.increase = 0.075, 
                     label = "p.adj", 
                     tip.length = 0.01, 
                     fontface = "bold", 
                     size = 4.5, 
                     bracket.size = 1.2) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(0,100,20)))
#ylim2 <- (layer_scales(gplot_3UTR_GC)$y$range$range[2]) * 1.05
#gplot_3UTR_GC <- gplot_3UTR_GC + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold", size = 4.5)

png(filename = "3'UTR_GC.png",
    res = 300,
    width = 1200, 
    height = 1500)
print(gplot_3UTR_GC)
dev.off()

##### GC3 #####

stats <- get_stats(seq_info, feature = "GC3_cont")

stats_lyst[["GC3"]] <- stats

p_val <- anova_fun(seq_info, feature = "GC3_cont")
padjs<- TUKEY_to_df(seq_info, feature = "GC3_cont")

gplot_GC3 <- 
  ggplot(seq_info, aes(x = group, y = GC3_cont * 100)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, 
               fill = c("#1B9E77", "#D95F02", "#7570B3"), 
               linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean * 100), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Depl."), # \nn=", sum(sequence_info_lengths$group == "downregulated")/3),
                              paste0("Unch."), # \nn=", sum(sequence_info_lengths$group == "unchanged")/3),
                              paste0("Enrich."))) + # \nn=", sum(sequence_info_lengths$group == "upregulated")/3))) +
  #scale_y_continuous(trans = c("log10"), labels = label_comma()) +
  ylab("GC3%") +
  xlab("Group") +
  ggtitle("GC3") +
  publication_theme() + 
  #ggtitle("GC3") +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

save_path = "\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/plots/transcript_features"

# png(file.path(save_path, "selectiveCNOT3_GC3.png"), res = 300, height = 1000, width = 1000)
# print(gplot_GC3)
# dev.off()

ylim <- layer_scales(gplot_GC3)$y$range$range[2] * 1.05
gplot_GC3 <- gplot_GC3 + 
  stat_pvalue_manual(padjs, 
                     y.position = ylim, 
                     step.increase = 0.075, 
                     label = "p.adj", 
                     tip.length = 0.01, 
                     fontface = "bold",
                     bracket.size = 1.2) +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(0,100,20)))
#ylim2 <- (layer_scales(gplot_GC3)$y$range$range[2]) * 1.05
#gplot_GC3 <- gplot_GC3 + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold", size = 4.5)

png(filename = "GC3.png",
    res = 300,
    width = 1200, 
    height = 1500)
print(gplot_GC3)
dev.off()


p_val <- kruskal_fun(seq_info, feature = "GC3_cont")
padjs<- kruskal_post_hoc(seq_info, feature = "GC3_cont")

gplot_GC3_kruskal <- 
  ggplot(seq_info, aes(x = group, y = GC3_cont * 100)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean * 100), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Downregulated"), # \nn=", sum(sequence_info_lengths$group == "downregulated")/3),
                              paste0("Unchanged"), # \nn=", sum(sequence_info_lengths$group == "unchanged")/3),
                              paste0("Upregulated"))) + # \nn=", sum(sequence_info_lengths$group == "upregulated")/3))) +
  #scale_y_continuous(trans = c("log10"), labels = label_comma()) +
  ylab("GC3%") +
  xlab("Group") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 11))

ylim <- layer_scales(gplot_GC3_kruskal)$y$range$range[2] * 1.02
gplot_GC3_kruskal <- gplot_GC3_kruskal + 
  stat_pvalue_manual(padjs, 
                     y.position = ylim, 
                     step.increase = 0.075, 
                     label = "p.adj", 
                     tip.length = 0.01, 
                     fontface = "bold", 
                     size = 4,
                     bracket.size = 1.2)
#ylim2 <- (layer_scales(gplot_GC3_kruskal)$y$range$range[2]) * 1.05
#gplot_GC3_kruskal <- gplot_GC3_kruskal + annotate("text",x = 1, y = ylim2, label = paste0("Kruskal-wallis, p = ", p_val), fontface = "bold", size = 3.75)

tiff(filename = "GC3_KW.tiff",
     units = "in",
     width = 4.3, 
     height = 4.3,
     res = 750)
print(gplot_GC3_kruskal)
dev.off()


####density plot ####
GC3_content2 <- seq_info %>%
  filter(!group == "unchanged")

mean_GC3 <- GC3_content2 %>%
  group_by(group) %>%
  summarise(meanGC3 = mean(GC3_cont, na.rm = T))

mean_downreg_GC3 <- as.numeric(formatC(mean_GC3$meanGC3[1], format = "f" ,digits = 3))
mean_upreg_GC3 <- as.numeric(formatC(mean_GC3$meanGC3[2], format = "f" ,digits = 3))

median_GC3 <- GC3_content2 %>%
  group_by(group) %>%
  summarise(medianGC3 = median(GC3_cont, na.rm = T))
median_downreg_GC3 <- as.numeric(formatC(median_GC3$medianGC3[1], format = "f" ,digits = 3))
median_upreg_GC3 <- as.numeric(formatC(median_GC3$medianGC3[2], format = "f" ,digits = 3))


density_plot <- ggplot(GC3_content2, aes(x = GC3_cont, colour = group)) + 
  geom_density(show.legend = F, linewidth = 1.5) + 
  stat_density(geom = "line", position="identity", linewidth = 1.5) +
  scale_color_manual(labels = c("Depleted", "Enriched"), values = c("#1B9E77", "#7570B3")) +
  #geom_vline(xintercept = mean_downreg_GC3, colour = "#1B9E77", linetype="dashed", size = 2)+
  geom_vline(xintercept = median_downreg_GC3, colour = "#1B9E77", linetype="dashed", size = 1.5)+
  #geom_vline(xintercept = mean_upreg_GC3, colour = "#7570B3", linetype="dashed", size = 2)+
  geom_vline(xintercept = median_upreg_GC3, colour = "#7570B3", linetype="dashed", size = 1.5)+
  xlab("GC3 Percent") +
  publication_theme() + 
  theme(legend.title = element_blank())

tiff(filename = file.path("GC3_density.tiff"),
     units = "in",
     width = 4, 
     height = 4,
     res = 700)
print(density_plot)
dev.off()

##### AG #####

stats <- get_stats(seq_info, feature = "AG_content")

stats_lyst[["AG"]] <- stats

p_val <- anova_fun(seq_info, feature = "AG_content")
padjs<- TUKEY_to_df(seq_info, feature = "AG_content")

gplot_AG <- 
  ggplot(seq_info, aes(x = group, y = AG_content * 100)) + 
  geom_violin(position = position_dodge(0.9), width=0.70, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean * 100), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(paste0("Depleted"), # \nn=", sum(sequence_info_lengths$group == "downregulated")/3),
                              paste0("Unchanged"), # \nn=", sum(sequence_info_lengths$group == "unchanged")/3),
                              paste0("Enriched"))) + # \nn=", sum(sequence_info_lengths$group == "upregulated")/3))) +
  #scale_y_continuous(trans = c("log10"), labels = label_comma()) +
  ylab("CDS AG%") +
  xlab("Group") +
  ggtitle("CDS AG%") +
  publication_theme() + 
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank())
        # axis.text.x = element_text(face = "bold", size = 11))

ylim <- layer_scales(gplot_AG)$y$range$range[2] * 1.05
gplot_AG <- gplot_AG + 
  stat_pvalue_manual(padjs, 
                     y.position = ylim, 
                     step.increase = 0.06, 
                     label = "p.adj", 
                     tip.length = 0.01, 
                     fontface = "bold",
                     bracket.size = 1.2)
#ylim2 <- (layer_scales(gplot_AG)$y$range$range[2]) * 1.04
#gplot_AG <- gplot_AG + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold", size = 4.5)

tiff(filename = "AG.tiff",
     width = 1000, 
     height = 1100,
     res = 200)
print(gplot_AG)
dev.off()





#### RSCU ####
seq_info_with_sequences <-
  inner_join(seq_info, master %>% select(ENSG, ENST, nucleotide_sequence_CDS), by = c("ENSG", "ENST"))

seq_info_with_sequences$nucleotide_sequence_CDS <- 
  strsplit(seq_info_with_sequences$nucleotide_sequence_CDS, split = "")

downregulated <- seq_info_with_sequences %>%
  filter(group == "downregulated")

upregulated <- seq_info_with_sequences %>%
  filter(group == "upregulated")

downregulated_names <- pull(downregulated, var = ENST)
downregulated_seqs <- pull(downregulated, var = nucleotide_sequence_CDS)
names(downregulated_seqs) <- downregulated_names

upregulated_names <- pull(upregulated, var = ENST)
upregulated_seqs <- pull(upregulated, var = nucleotide_sequence_CDS)
names(upregulated_seqs) <- upregulated_names

RSCU_lyst <- 
  list(lapply(downregulated_seqs, uco, frame = 0, index = c("eff", "freq", "rscu"), as.data.frame = T, NA.rscu = NA),
       lapply(upregulated_seqs,uco, frame = 0, index = c("eff", "freq", "rscu"), as.data.frame = T, NA.rscu = NA))

names(RSCU_lyst) <- c("downregulated", "upregulated")

data_upreg <- bind_rows(purrr::imap(RSCU_lyst[["upregulated"]], ~mutate(.x, transcript = .y))) %>%
  mutate(group = "upregulated")
data_downreg <- bind_rows(purrr::imap(RSCU_lyst[["downregulated"]], ~mutate(.x, transcript = .y))) %>%
  mutate(group = "downregulated")
row.names(data_upreg) <- NULL
row.names(data_downreg) <- NULL

data_combined <- rbind(data_downreg, data_upreg)


RCSUs <- data_combined %>%
  group_by(group, codon) %>%
  summarise(mean_RSCU = mean(RSCU, na.rm = T))

RCSUs <- RCSUs %>% 
  filter(codon != "taa" & codon != "tga" & codon != "tag") %>%
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A", "U", "G","C")))

RCSUs <- spread(RCSUs, group, mean_RSCU)
CTG_codon <- RCSUs %>% filter(codon == "ctg")
TTG_codon <- RCSUs %>% filter(codon == "ttg")

RCSUs %>% filter(downregulated > upregulated & wobble == "G")

gplot_RSCU <- ggplot(RCSUs, aes(x = downregulated, y = upregulated)) +
  geom_point(aes(colour = wobble), shape = 19, alpha = 0.75, size = 3) +
  ggrepel::geom_text_repel(data= CTG_codon, 
                           aes(x = downregulated,  y = upregulated, 
                               label = codon),
                           size = 8) +
  ggrepel::geom_text_repel(data= TTG_codon, 
                           aes(x = downregulated,  y = upregulated, 
                               label = codon),
                           size = 8, nudge_y = -0.3, nudge_x = +0.2) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), name = "Wobble Base Identity" )+
  geom_abline(slope = 1, intercept = 0, linewidth = 1, linetype = "dashed") +
  #ggtitle("Monosomes") +
  #xlim(0,3)+
  #ylim(0,2.5)+
  xlab("Depleted") + 
  ylab("Enriched") +
  ggtitle("RSCU") +
  publication_theme() +
  scale_x_continuous(limits = c(0,3)) +
  scale_y_continuous(limits = c(0,3)) +
  theme(legend.position = "bottom") +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(siz = 20),
    axis.title.x = element_text(siz = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20))

png(filename = "RSCU.png",
    res = 300,
    width = 1750, 
    height = 1750)
print(gplot_RSCU)
dev.off()

AGGC3RSCU <- ggarrange(gplot_AG, 
          gplot_GC3, 
          gplot_RSCU, 
          nrow = 1)

final_AGGC3RSCU <- annotate_figure(AGGC3RSCU, 
                                left = textGrob("Percentage", 
                                                rot = 90, vjust = 1,
                                                gp = gpar(cex = 2.5)))

tiff(filename = "final_AGGC3RSCU.tiff",
width = 3200, 
height = 1000,
res = 200)
print(final_AGGC3RSCU)
dev.off()

#### signal sequences ####
SS_info <- seq_info %>%
  group_by(group) %>%
  summarise(SS = sum(signal_sequence )/ length(signal_sequence ) * 100,
            no_SS = 100 - SS) %>%
  gather(key = "contains_SS", value = "percent", 2:3) %>% 
  filter(contains_SS != "no_SS")

"#FBB4AE" "#B3CDE3" "#CCEBC5" "#DECBE4" "#FED9A6" "#FFFFCC" "#E5D8BD" "#FDDAEC"
"#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"
"#FBB4AE" "" "#CCEBC5" "#DECBE4" "#FED9A6"
"" "#FDCDAC" "#CBD5E8" "#F4CAE4" "#E6F5C9" "#FFF2AE" "#F1E2CC" "#CCCCCC"

"#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF"

SS_plot <- SS_info %>%
  ggplot(aes(x = group, y = percent)) + 
  geom_bar(stat = "identity", fill = "#377EB8") +
  #scale_fill_manual(labels = c("signal sequence"), values = c("#377EB8")) +
  scale_x_discrete(labels = c("Depleted", "Unchanged", "Enriched")) +
  scale_y_continuous(limits = c(0, 11), breaks = c(seq(0,12,2))) +
  ggtitle("Signal sequences") +
  publication_theme() +
  ylab("% of transcripts") +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 18))

tiff(filename = "Signal_seqs.tiff",
     units = "in",
     width = 4, 
     height = 4,
     res = 200)
print(SS_plot)
dev.off()

#### uORFs ####
uORF_info <- seq_info %>%
  group_by(group) %>%
  summarise(uORF_positive = sum(uORF)/ length(uORF) * 100,
            uORF_negative = 100 - uORF_positive) %>%
  gather(key = "contains_uORF", value = "percent", 2:3)%>% 
  filter(contains_uORF != "uORF_negative")


uORF_plot <- uORF_info %>%
  ggplot(aes(x = group, y = percent)) + 
  geom_bar(stat = "identity", fill = "#377EB8") +
  #scale_fill_manual(labels = c("-uORF", "+uORF"), values = c("#377EB8")) +
  scale_x_discrete(labels = c("Downregulated", "Unchanged", "Upregulated")) +
  publication_theme() +
  ylab("% of transcripts") +
  ggtitle("uORFs") +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 18))

tiff(filename = "uORFs.tiff",
     units = "in",
     width = 4, 
     height = 4,
     res = 200)
print(uORF_plot)
dev.off()


#### miRNAs ####
seq_info_miRNAs <- seq_info %>%
  mutate(miRNAs = factor(case_when(miRNA_binding_sites <= 10 ~ "0-10",
                                   miRNA_binding_sites > 10 & miRNA_binding_sites <= 20 ~ "11-20",
                                   miRNA_binding_sites > 20 & miRNA_binding_sites <= 30 ~ "21-30",
                                   miRNA_binding_sites > 30 ~ "30+")))


miRNA_number <- ggplot(seq_info_miRNAs, aes(x = miRNAs, y = log2FoldChange, fill = miRNAs)) + 
  geom_violin() +
  geom_boxplot(position = position_dodge(0.9), width = 0.5) +
  scale_fill_manual(values = c("#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4"),
                    guide = "none") +
  ggtitle("miRNAs") +
  publication_theme()

tiff(filename = "miRNAs_by_number.tiff",
     units = "in",
     width = 4, 
     height = 4,
     res = 200)
print(miRNA_number)
dev.off()


stats <- get_stats(seq_info, feature = "miRNA_binding_sites")

stats_lyst[["miRNAs"]] <- stats

p_val <- anova_fun(seq_info, feature = "miRNA_binding_sites")
padjs<- TUKEY_to_df(seq_info, feature = "miRNA_binding_sites")

miRNA_gplot <- ggplot(seq_info_miRNAs, aes(x = group, y = miRNA_binding_sites)) +
  geom_violin(position = position_dodge(0.9), width=.7, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(trans = c("log10")) +
  scale_x_discrete(labels = c(paste0("Dep."), #; \nn=", sum(sequence_info$group == "downregulated")),
                              paste0("Unch."), #; \nn=", sum(sequence_info$group == "unchanged")),
                              paste0("Enr."))) + #; \nn=", sum(sequence_info$group == "upregulated")))) +
  xlab("Group") +
  ylab("miRNA binding sites") +
  ggtitle("microRNAs") +
  publication_theme() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())

ylim <- (layer_scales(miRNA_gplot)$y$range$range[2]) * 1.025
miRNA_gplot <- miRNA_gplot + stat_pvalue_manual(padjs, y.position = ylim, step.increase = 0.075, label = "p.adj", tip.length = 0.01, fontface = "bold", size = 4.5, bracket.size = 1.2)

#ylim2 <- 10^(layer_scales(miRNA_gplot)$y$range$range[2]) * 1.4
#miRNA_gplot <- miRNA_gplot + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val), fontface = "bold", size = 4.5)


png(filename = "microRNAs.png",
    res = 300,
    width = 1200, 
    height = 1500)
print(miRNA_gplot)
dev.off()

p_val <- kruskal_fun(seq_info, feature = "miRNA_binding_sites")
padjs<- kruskal_post_hoc(seq_info, feature = "miRNA_binding_sites")

miRNA_gplot_kruskal <- ggplot(seq_info_miRNAs, aes(x = group, y = miRNA_binding_sites)) +
  geom_violin(position = position_dodge(0.9), width=.7, aes(fill = group), linewidth = 0.75) +
  geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
  geom_point(data = stats, aes(x = group, y = mean), colour = "black", size = 3, shape = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(trans = c("log10")) +
  scale_x_discrete(labels = c(paste0("Downregulated"), #; \nn=", sum(sequence_info$group == "downregulated")),
                              paste0("Unchanged"), #; \nn=", sum(sequence_info$group == "unchanged")),
                              paste0("Upregulated"))) + #; \nn=", sum(sequence_info$group == "upregulated")))) +
  xlab("Group") +
  ylab("miRNA binding sites") +
  publication_theme() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank())

ylim <- (layer_scales(miRNA_gplot_kruskal)$y$range$range[2]) * 1.025
miRNA_gplot_kruskal <- miRNA_gplot_kruskal + stat_pvalue_manual(padjs, y.position = ylim, step.increase = 0.045, label = "p.adj", tip.length = 0.01, fontface = "bold", size = 4)
ylim2 <- 10^(layer_scales(miRNA_gplot_kruskal)$y$range$range[2]) * 1.5
miRNA_gplot_kruskal <- miRNA_gplot_kruskal + annotate("text",x = 1, y = ylim2, label = paste0("Kruskal Wallis, p = ", p_val), fontface = "bold", size = 3.5)

tiff(filename = "miRNAs_KW.tiff",
     units = "in",
     width = 6.6, 
     height = 6.6,
     res = 200)
print(miRNA_gplot_kruskal)
dev.off()


#### Global_AA_usage ####

aminos <- c(unique(codons$AA))
aminos_plus_Stp <- append(aminos, "Stp")

aa_usage_lyst <- list()

k = 5

for(k in 1:length(aminos_plus_Stp)) {
  
  AA <- aminos_plus_Stp[k] 
  
  print(AA)
  
  int <- seq_info %>%
    select(group, AA) %>%
    gather(key = "amino_acid", value = "transcript_percentage", 2)
  
  stats <- get_stats(int, feature = "transcript_percentage")
  stats_lyst[[AA]] <- stats
  
  p_val <- anova_fun(int, feature = "transcript_percentage")
  padjs<- TUKEY_to_df(int, feature = "transcript_percentage")
  
  gplot <- int %>%
    ggplot(aes(x = group, y = transcript_percentage)) + 
    geom_violin(position = position_dodge(0.9), width=.7, aes(fill = group), linewidth = 0.75) +
    geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
    geom_point(data = stats, aes(x = group, y = mean), colour = "black", size = 3, shape = 15) +
    scale_y_continuous(trans = c("log10")) +
    scale_fill_brewer(palette = "Set2") +
    scale_x_discrete(labels = c(paste0("Downregulated"), #; \nn=", sum(sequence_info$group == "downregulated")),
                                paste0("Unchanged"), #; \nn=", sum(sequence_info$group == "unchanged")),
                                paste0("Upregulated"))) +  #; \nn=", sum(sequence_info$group == "upregulated")))) +
    xlab("Group") +
    ylab(paste0(AA, "%")) +
    
    publication_theme() +
    theme(legend.position = "none") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.y = element_text(siz = 20))
  
  
  
  
  ylim <- (layer_scales(gplot)$y$range$range[2]) * 1.05
  gplot <- gplot + stat_pvalue_manual(padjs, 
                                      y.position = ylim, 
                                      step.increase = 0.06, 
                                      label = "p.adj", 
                                      tip.length = 0.01, 
                                      fontface = "bold",
                                      size = 4)
  ylim2 <- 10^(layer_scales(gplot)$y$range$range[2]) * 1.5
  gplot <- gplot + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val),
                            fontface = "bold", size = 4)
  
  tiff(filename = file.path("features/AA_usage", paste0(AA,".tiff")),
       units = "in",
       width = 6.5, 
       height = 6.5,
       res = 200)
  print(gplot)
  dev.off()
}



#### Global codon usage ####


for(k in 1:length(code_ons_plus_stop)) {
  
  codon <- code_ons_plus_stop[k] 
  
  print(codon)
  
  int <- seq_info %>%
    select(group, codon) %>%
    gather(key = "codon", value = "transcript_percentage", 2)
  
  stats <- get_stats(int, feature = "transcript_percentage")
  stats_lyst[[codon]] <- stats
  
  p_val <- anova_fun(int, feature = "transcript_percentage")
  padjs<- TUKEY_to_df(int, feature = "transcript_percentage")
  
  gplot <- int %>%
    ggplot(aes(x = group, y = transcript_percentage)) + 
    geom_violin(position = position_dodge(0.9), width=.7, aes(fill = group), linewidth = 0.75) +
    geom_boxplot(position = position_dodge(0.9), width = 0.3, fill = c("#1B9E77", "#D95F02", "#7570B3"), linewidth = 0.75) +
    geom_point(data = stats, aes(x = group, y = mean), colour = "black", size = 3, shape = 15) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    scale_fill_brewer(palette = "Set2") +
    scale_x_discrete(labels = c(paste0("Downregulated"), #; \nn=", sum(sequence_info$group == "downregulated")),
                                paste0("Unchanged"), #; \nn=", sum(sequence_info$group == "unchanged")),
                                paste0("Upregulated"))) +  #; \nn=", sum(sequence_info$group == "upregulated")))) +
    xlab("Group") +
    ylab(paste0(codon, "%")) +
    
    publication_theme() +
    theme(legend.position = "none") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.y = element_text(siz = 20))
  
  ylim <- (layer_scales(gplot)$y$range$range[2]) * 1.05
  gplot <- gplot + stat_pvalue_manual(padjs, 
                                      y.position = ylim, 
                                      step.increase = 0.06, 
                                      label = "p.adj", 
                                      tip.length = 0.01, 
                                      fontface = "bold", 
                                      size = 4)
  ylim2 <- 10^(layer_scales(gplot)$y$range$range[2]) * 1.5
  gplot <- gplot + annotate("text",x = 1, y = ylim2, label = paste0("Anova, p = ", p_val),
                            fontface = "bold",
                            size = 4)
  
  tiff(filename = file.path("features/codon_usage", paste0(codon,".tiff")),
       units = "in",
       width = 6.5, 
       height = 6.5,
       res = 200)
  print(gplot)
  dev.off()
}

all_stats <- do.call("rbind", stats_lyst) %>%
  mutate(feature = rep(names(stats_lyst), each = 3))

write.csv(all_stats, file = file.path(parent_dir, "plots/stats/all_stats_M.csv"))


tmhmms <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/TM_domains_FINAL.csv") %>% 
  select(ENST, class) %>% unique()

secretion_volcano <- 
DEseq2 %>% 
  select(transcript, log2FoldChange, padj, group) %>% 
  left_join(tmhmms, by = c("transcript" = "ENST"))%>% 
  mutate(class = case_when(is.na(class) ~ "cytosolic",.default = "ER")) 
  


ggplot() +
  geom_point(data = secretion_volcano %>% filter(class == "cytosolic"), 
             aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.25) +
  geom_point(data = secretion_volcano %>% filter(class == "ER"), 
             aes(x = log2FoldChange, y = -log10(padj)), colour = "blue")
  

secretory_info <- 
  seq_info %>% 
  left_join(tmhmms, by = c( "ENST")) %>% 
  select(group, ENST, class) %>% 
  mutate(class = case_when(is.na(class) ~ "cytosolic",.default = "ER")) 

tmhmm_data <- 
  inner_join(
    secretory_info %>% 
      group_by(group, class) %>% 
      tally(),
    secretory_info %>% 
      group_by(group) %>% 
      tally() %>% 
      rename(sum_n = n),
    by = "group") %>% 
  mutate(per = n/sum_n * 100) %>% 
  filter(class == "ER") %>% 
  mutate(group = factor(case_when(group == "downregulated" ~ "Dep.",
                           group == "upregulated" ~ "Up.",
                           group == "unchanged" ~ "Unch."), 
                        levels = c("Dep.", "Unch.", "Up."))) %>% 
  ggplot(aes(x = group, y = per, fill = group, label = n)) +
  geom_bar(stat = "identity") +
  geom_text(vjust = -0.1, fontface = "bold") +
  scale_fill_brewer(palette = "Dark2") +
  ylab("Percentage") +
  ggtitle("SP/TM presence") +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

tmhmm_data

png("ER_targets.png", res = 300, height = 1200, width = 1200)
print(tmhmm_data)
dev.off()
