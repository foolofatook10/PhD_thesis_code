#load libraries----
library(tidyverse)
library(fgsea)


publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = 22, hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = 20),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(size = 16, color = "black"),
            axis.text.y = element_text(size = 14, color = "black"),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = 12),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(1, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(5,5,5,20),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

# setwd()
setwd("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL")

#read in and set common variables----
source("R_scripts/common_variables.R")

#create a variable for what the treatment is
control <- "Tot"
treatment <- "CNOT3"

#set the seed to ensure reproducible results
set.seed(020588)

#themes----
mytheme <- theme_classic()+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.position = "none")

#functions----
run_fgsea <- function(named_vector, pathway) {
  results <- fgsea(pathways = pathway,
                   stats=named_vector,
                   minSize = 20,
                   maxSize = 1000)
  return(results)
}

make_plot <- function(fgsea_result, padj_threshold, title) {
  
  plot <- ggplot(data = fgsea_result[fgsea_result$padj < padj_threshold], aes(reorder(pathway, NES), NES)) +
    geom_col() +
    coord_flip() +
    labs(x="Pathway", y="NES",
         title=title) + 
    # theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none"
          # axis.text.y = element_text(size = 10))
    ) +
    publication_theme()
  return(plot)
}

make_plot2 <- function(fgsea_result, padj_threshold, title) {
  
  df_significant <- 
    fgsea_result %>% 
    filter(padj <0.05) %>% 
    arrange(-NES)
  
  final_df <- 
    rbind(
      df_significant %>% head(10),
      df_significant %>% tail(10))
  
  terms <- 
    unlist(
      lapply(str_split(final_df$pathway, "_"), 
             function(x){paste0(x %>% tail(-1), collapse = " ")}))
  
  final_df$pathway <- terms
  
  # final_df %>% 
  #   mutate(case_when(pathway == "PROTEOGLYCAN BIOSYNTHETIC PROCESS" ~ ))
  
  plot <- ggplot(data = final_df, aes(reorder(pathway, NES), NES)) +
    geom_col() +
    coord_flip() +
    labs(x="Pathway", y="NES",
         title=title) + 
    # theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none"
          # axis.text.y = element_text(size = 10))
    ) +
    scale_y_discrete(labels = \(x) str_wrap(x, width = 20)) +
    publication_theme()
  return(plot)
}

make_plot3 <- function(fgsea_result, padj_threshold, title) {
  
  df_significant <- 
    fgsea_result %>% 
    filter(padj <0.05) %>% 
    arrange(-NES)
  
  final_df <- df_significant
  
  terms <- 
    unlist(
      lapply(str_split(final_df$pathway, "_"), 
             function(x){paste0(x %>% tail(-1), collapse = " ")}))
  
  final_df$pathway <- terms
  
  plot <- ggplot(data = final_df, aes(reorder(pathway, NES), NES)) +
    geom_col() +
    coord_flip() +
    labs(x="Pathway", y="NES",
         title=title) + 
    # theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none"
          # axis.text.y = element_text(size = 10))
    ) +
    publication_theme()
  return(plot)
}

#read in DESeq2 output----
DESeq2_data <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output/RPFs_CNOT3_DEseq2_apeglm_LFC_shrinkage_M.csv"))

#make named vectors----
### This will make a list of named vectors to allow gsea to be carried out separately on the RPFs, Totals and TE log2FC
DESeq2_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(log2FoldChange)) %>%
  deframe() -> RPFs_named_vector


named_vectors <- list(RPFs = RPFs_named_vector)

#read in pathways----
source("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/read_human_GSEA_pathways.R") #This may need to be changed to human

padj <- 0.05


#### biological process

#biological processes----
#carry out fgsea
bio_processes_results <- lapply(named_vectors, run_fgsea, pathway = pathways.bio_processes)


bio_processes_results_df <- 
  bio_processes_results$RPFs %>% 
  filter(NES > 0 & padj < 0.05)

le_bio_process <- bio_processes_results_df$leadingEdge
names(le_bio_process) <- bio_processes_results_df$pathway

fgsea_result = bio_processes_results$RPFs

df_significant <- 
  fgsea_result %>% 
  filter(padj <0.05) %>% 
  arrange(-NES)

final_df <- 
  rbind(
    df_significant %>% head(20))

terms <- 
  unlist(
    lapply(str_split(final_df$pathway, "_"), 
           function(x){paste0(x %>% tail(-1), collapse = " ")}))

final_df$pathway <- terms

# final_df %>% 
#   mutate(case_when(pathway == "PROTEOGLYCAN BIOSYNTHETIC PROCESS" ~ ))

bio_processes <- 
ggplot(final_df, aes(x = NES, reorder(pathway, NES))) + 
  geom_col() + 
  scale_y_discrete( labels = \(x) str_wrap(x, width = 20) )+
  labs(y="Pathway", x="Normalised Enrichment Score") +
  ggtitle("Biological process") +
  publication_theme() +
  theme(axis.title.y = element_blank())



#cellular component----
#carry out fgsea
cell_comp_results <- lapply(named_vectors, run_fgsea, pathway = pathways.cell_comp)




cell_comp_results$RPFs$pathway <- str_remove(cell_comp_results$RPFs$pathway,
                                             "GO_")
cell_comp_results$RPFs$pathway <- str_replace_all(cell_comp_results$RPFs$pathway, "_", " ")



fgsea_result = cell_comp_results$RPFs

df_significant <- 
  fgsea_result %>% 
  filter(padj <0.05) %>% 
  arrange(-NES)

final_df <- 
  rbind(
    df_significant %>% head(20))

cellular_component <- 
ggplot(final_df, aes(x = NES, reorder(pathway, NES))) + 
  geom_col() + 
  scale_y_discrete( labels = \(x) str_wrap(x, width = 20) )+
  labs(y="Pathway", x="Normalised Enrichment Score") +
  ggtitle("Cellular component") +
  publication_theme() +
  theme(axis.title.y = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/GO_terms.png",
    res = 300, height = 4500, width = 4000)
cowplot::plot_grid(cellular_component, bio_processes, nrow = 1, align = "v")
dev.off()
