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