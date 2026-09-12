#This is written for mouse data, will need to read in human pathways and edit pathway names if to be run on human data

#load libraries----
library(tidyverse)
library(fgsea)


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
            axis.text.y = element_text(size = 16, color = "black"),
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
            plot.margin=unit(c(10,5,5,5),"mm"),
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

#hallmark----
#carry out fgsea
hallmark_results <- lapply(named_vectors, run_fgsea, pathway = pathways.hallmark)

#save results
save(file = file.path("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/GSEA/CNOT3SelRP/hallmarks.RData"), hallmark_results)

padj <- 0.05

le_hallmarks <- 
hallmark_results$RPFs %>% 
  filter(padj < 0.05 & NES > 0) %>% 
  pull(leadingEdge)

le_hallmarks_names <- 
hallmark_results$RPFs %>% 
  filter(padj < 0.05 & NES > 0) %>% 
  pull(pathway)

names(le_hallmarks) <- le_hallmarks_names

hallmark_results$RPFs$pathway <- str_remove(hallmark_results$RPFs$pathway, "HALLMARK_")
hallmark_results$RPFs$pathway <- str_replace_all(hallmark_results$RPFs$pathway, "_", " ")

hallmark_results$RPFs %>% filter(padj < 0.05)

#plot enriched pathways
png(filename = file.path(parent_dir, 
                         "plots/fgsea/", 
                         paste(treatment, "RPFs_hallmark.png", sep = "_")), 
    width = 2500, height = 2500, res = 300)
make_plot(fgsea_result = hallmark_results$RPFs, padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Hallmark genes"))
dev.off()

# transcription factors
hallmark_results <- lapply(named_vectors, run_fgsea, pathway = pathways.transcription_factors)

#save results
save(file = file.path("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/GSEA/CNOT3SelRP/transcription_factors.RData"), hallmark_results)

hallmark_results_NESup <- 
hallmark_results$RPFs %>% 
  filter(padj < 0.05 & NES > 0)

hallmark_results_NESdown <- 
  hallmark_results$RPFs %>% 
  filter(padj < 0.05 & NES < 0)

hallmark_results_NESup$pathway
hallmark_results_NESdown$pathway

make_plot(fgsea_result = hallmark_results_NESup, padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA TF gene sets"))

# microRNAs
# transcription factors
hallmark_results <- lapply(named_vectors, run_fgsea, pathway = pathways.miRNA_targets)

#save results
save(file = file.path("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/GSEA/CNOT3SelRP/microRNAs.RData"), hallmark_results)

hallmark_results_NESup <- 
  hallmark_results$RPFs %>% 
  filter(padj < 0.05 & NES > 0) %>% 
  arrange(-NES)

hallmark_results_NESup$pathway
le_microRNAs <- hallmark_results_NESup$leadingEdge
names(le_microRNAs) <- hallmark_results_NESup$pathway

hallmark_results_NESup %>% slice_head(n = 40)
length(unlist(hallmark_results_NESup[3,]$leadingEdge))

#biological processes----
#carry out fgsea
bio_processes_results <- lapply(named_vectors, run_fgsea, pathway = pathways.bio_processes)

#save results
save(file = file.path("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/GSEA/CNOT3SelRP/biological_processes.RData"), bio_processes_results)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/bio_processes_results.Rdata"), bio_processes_results)

#set adjusted p-value
padj <- 0.05

bio_processes_results_df <- 
bio_processes_results$RPFs %>% 
  filter(NES > 0 & padj < 0.05)

le_bio_process <- bio_processes_results_df$leadingEdge
names(le_bio_process) <- bio_processes_results_df$pathway


#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_bio_processes.png", sep = "_")), res = 300, width = 4200, height = 2500)
make_plot2(fgsea_result = bio_processes_results$RPFs, padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Biological Processes genes"))
dev.off()

bio_processes_results$RPFs %>% 
  filter(padj <0.05) %>% 
  arrange(desc(NES)) %>% 
  filter(NES > 2) %>% 
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  publication_theme()
    #plot.title = element_text(hjust = 0.5),
       # axis.text.y = element_text(size = 5))



#molecular functions----
#carry out fgsea
mol_funs_results <- lapply(named_vectors, run_fgsea, pathway = pathways.mol_funs)

#save results
save(file = file.path("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/GSEA/CNOT3SelRP/molecular_functions.RData"), mol_funs_results)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/mol_funs_results.Rdata"), mol_funs_results)

#set adjusted p-value
padj <- 0.05

le_molfun <- mol_funs_results$RPFs$leadingEdge
names(le_molfun) <- mol_funs_results$RPFs$pathway

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_mol_funs.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = mol_funs_results$RPFs, padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Molecular Functions gene sets"))
dev.off()


#cellular component----
#carry out fgsea
cell_comp_results <- lapply(named_vectors, run_fgsea, pathway = pathways.cell_comp)

#save results
save(file = file.path("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/GSEA/CNOT3SelRP/cellular_component.RData"), cell_comp_results)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/cell_comp_results.Rdata"), cell_comp_results)

#set adjusted p-value
padj <- 0.05

le_cellcomp <- cell_comp_results$RPFs$leadingEdge
names(le_cellcomp) <- cell_comp_results$RPFs$pathway

cell_comp_results$RPFs$pathway <- str_remove(cell_comp_results$RPFs$pathway,
                                             "GO_")
cell_comp_results$RPFs$pathway <- str_replace_all(cell_comp_results$RPFs$pathway, "_", " ")

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_cell_comp.png", sep = "_")), res = 300, width = 4000, height = 3000)
make_plot3(fgsea_result = cell_comp_results$RPFs %>% filter(NES > 0), padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Cellular Component gene sets"))
dev.off()

tidy_data_bio_process <- 
function(x){
  
  genes <- le_bio_process[[x]]
  term <- names(le_bio_process)[x]
  
  df <- 
  data.frame(
    GSEA = "bio_process",
    term = term,
     genes = genes)
  
  return(df) 
}




tidy_data_cell_comp <- 
  function(x){
    
    genes <- le_cellcomp[[x]]
    term <- names(le_cellcomp)[x]
    
    df <- 
      data.frame(
        GSEA = "cell_comp",
        term = term,
                 genes = genes)
    
    return(df) 
  }

tidy_data_hallmarks <- 
  function(x){
    
    genes <- le_hallmarks[[x]]
    term <- names(le_hallmarks)[x]
    
    df <- 
      data.frame(
        GSEA = "hallmarks",
        term = term,
                 genes = genes)
    
    return(df) 
  }

tidy_data_molfun <- 
  function(x){
    
    genes <- le_molfun[[x]]
    term <- names(le_molfun)[x]
    
    df <- 
      data.frame(
        GSEA = "molfun",
        term = term,
                 genes = genes)
    
    return(df) 
  }

tidy_data_microRNAs <- 
  function(x){
    
    genes <- le_microRNAs[[x]]
    term <- names(le_microRNAs)[x]
    
    df <- 
      data.frame(
        GSEA = "microRNA",
        term = term,
        genes = genes)
    
    return(df) 
  }

all_leading_edges <- 
rbind(
do.call("rbind", lapply(1:length(le_bio_process), tidy_data_bio_process)),
do.call("rbind", lapply(1:length(le_cellcomp), tidy_data_cell_comp)),
do.call("rbind", lapply(1:length(le_hallmarks), tidy_data_hallmarks)),
do.call("rbind", lapply(1:length(le_molfun), tidy_data_molfun)),
do.call("rbind", lapply(1:length(le_microRNAs), tidy_data_microRNAs)))

write.csv(all_leading_edges, file = "C:/Users/Jettles/OneDrive - University of Glasgow/Documents/CNOT3SelRP_leading_edge_genes.csv")


#Curated----
### the curated is specific to the mouse GSEA lists
#read in pathways
#carry out fgsea
curated_results <- lapply(named_vectors, run_fgsea, pathway = pathways.curated)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/curated_results.Rdata"), curated_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_curated.png", sep = "_")), width = 1000, height = 2000)
make_plot(fgsea_result = curated_results$RPFs, padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Curated gene sets"))
dev.off()
