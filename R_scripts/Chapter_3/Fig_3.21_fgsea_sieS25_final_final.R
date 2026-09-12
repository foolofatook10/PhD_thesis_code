#This is written for mouse data, will need to read in human pathways and edit pathway names if to be run on human data

#load libraries----
library(tidyverse)
library(fgsea)

setwd("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/R_scripts")

#read in and set common variables----
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/R_scripts/common_variables.R")

#create a variable for what the treatment is
control <- "NTC"
treatment <- "sieS25"

# read in publication theme
publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = 24, face = "italic", hjust = 0.5),
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
            legend.key.size= unit(1, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(2,10,2,25),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

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
  
  
  fgsea_result$pathway <- str_remove(str_replace_all(fgsea_result$pathway, "_", " "), "^\\w+")
  
  plot <- ggplot(data = fgsea_result[fgsea_result$padj < padj_threshold], aes(reorder(pathway, NES), NES)) +
    geom_col(fill = "#984EA3") +
    coord_flip() +
    labs(x="Pathway", y="NES",
         title=title) + 
    publication_theme() +
    theme(axis.title.y = element_blank())
  return(plot)
}

make_plot2 <- function(fgsea_result, padj_threshold, title) {
  
  fgsea_result$pathway <- str_remove(str_replace_all(fgsea_result$pathway, "_", " "), "^\\w+")
  
  sigresults <- 
    fgsea_result %>% filter(padj < 0.05) %>% arrange(NES)
  
  sigresults_up <- sigresults %>% head(n = 5)
  sigresults_down <- sigresults %>% tail(n = 5)
  
  all <- rbind(sigresults_up, sigresults_down)
  
  all2 <- 
    all %>% 
    mutate(
      pathway = case_when(
        pathway == " COTRANSLATIONAL PROTEIN TARGETING TO MEMBRANE" ~ "COTRANSLATIONAL PROTEIN\nTARGETING TO MEMBRANE",
        pathway == " PROTEIN LOCALIZATION TO ENDOPLASMIC RETICULUM" ~ "PROTEIN LOCALIZATION TO\nENDOPLASMIC RETICULUM",
        pathway == " ESTABLISHMENT OF PROTEIN LOCALIZATION TO ENDOPLASMIC RETICULUM" ~ "ESTABLISHMENT OF\nPROTEIN LOCALIZATION TO\nENDOPLASMIC RETICULUM",
        pathway == " DNA REPLICATION DEPENDENT NUCLEOSOME ORGANIZATION" ~ "DNA REPLICATION DEPENDENT\nNUCLEOSOME ORGANIZATION",
        pathway == " NUCLEAR TRANSCRIBED MRNA CATABOLIC PROCESS NONSENSE MEDIATED DECAY" ~ "NUCLEAR TRANSCRIBED MRNA\nCATABOLIC PROCESS NONSENSE\nMEDIATED DECAY",
        pathway == " DNA REPLICATION DEPENDENT NUCLEOSOME ORGANISATION" ~ " DNA REPLICATION DEPENDENT\nNUCLEOSOME ORGANISATION",
        pathway == " VASCULAR PROCESS IN CIRCULATORY SYSTEM" ~ " VASCULAR PROCESS\nIN CIRCULATORY SYSTEM",
        pathway == " MEMBRANE LIPID METABOLIC PROCESS" ~ " MEMBRANE LIPID\nMETABOLIC PROCESS",
        pathway == " GPI ANCHOR METABOLIC PROCESS" ~ " GPI ANCHOR\nMETABOLIC PROCESS",
        pathway == " LIPOSACCHARIDE METABOLIC PROCESS" ~ " LIPOSACCHARIDE\nMETABOLIC PROCESS",
        pathway == " GLYCOLIPID BIOSYNTHETIC PROCESS" ~ " GLYCOLIPID\nBIOSYNTHETIC PROCESS",
        pathway == " CATALYTIC ACTIVITY ACTING ON A TRNA" ~ "CATALYTIC ACTIVITY\nACTING ON A TRNA",
        pathway == " CATALYTIC ACTIVITY ACTING ON RNA" ~ "CATALYTIC ACTIVITY\nACTING ON RNA",
        pathway == " RNA METHYLTRANSFERASE ACTIVITY" ~ " RNA METHYLTRANSFERASE\nACTIVITY",
        pathway == " TRNA METHYLTRANSFERASE ACTIVITY" ~ " TRNA METHYLTRANSFERASE\nACTIVITY",
        pathway == " TRANSMEMBRANE RECEPTOR PROTEIN KINASE ACTIVITY" ~ "TRANSMEMBRANE RECEPTOR\nPROTEIN KINASE ACTIVITY",
        pathway == " UDP GLYCOSYLTRANSFERASE ACTIVITY" ~ " UDP GLYCOSYLTRANSFERASE\nACTIVITY",
        pathway == " TRANSFERASE ACTIVITY TRANSFERRING HEXOSYL GROUPS" ~ "TRANSFERASE ACTIVITY\nTRANSFERRING HEXOSYL\nGROUPS",
        pathway == " TUMOR NECROSIS FACTOR RECEPTOR SUPERFAMILY BINDING" ~ " TUMOR NECROSIS FACTOR\nRECEPTOR SUPERFAMILY BINDING",
        pathway == " PRERIBOSOME LARGE SUBUNIT PRECURSOR" ~ " PRERIBOSOME LARGE\nSUBUNIT PRECURSOR",
        pathway == " CATALYTIC STEP 2 SPLICEOSOME" ~ " CATALYTIC STEP\n2 SPLICEOSOME",
        pathway == " INTRINSIC COMPONENT OF GOLGI MEMBRANE" ~ " INTRINSIC COMPONENT\nOF GOLGI MEMBRANE",
        pathway == " INTRINSIC COMPONENT OF PLASMA MEMBRANE" ~ " INTRINSIC COMPONENT OF\nPLASMA MEMBRANE",
        pathway == " RIBONUCLEOPROTEIN COMPLEX" ~ " RIBONUCLEOPROTEIN\nCOMPLEX",
        pathway == " ENDOPLASMIC RETICULUM LUMEN" ~ " ENDOPLASMIC RETICULUM\nLUMEN",
        pathway == " CYTOSOLIC LARGE RIBOSOMAL SUBUNIT" ~ " CYTOSOLIC LARGE\nRIBOSOMAL SUBUNIT",
        pathway == " SMALL SUBUNIT PROCESSOME" ~ "SMALL SUBUNIT\nPROCESSOME",
        pathway == " INTRINSIC COMPONENT OF SYNAPTIC VESICLE MEMBRANE" ~ " INTRINSIC COMPONENT OF\nSYNAPTIC VESICLE MEMBRANE",
        pathway == " NCRNA PROCESSING" ~ " NCRNA\nPROCESSING",
        pathway == " NCRNA METABOLIC PROCESSING" ~ " NCRNA\nMETABOLIC PROCESSING",
        pathway == " RRNA METABOLIC PROCESSING" ~ " RRNA\nMETABOLIC PROCESSING",
        pathway == " CHROMATIN SILENCING AT RDNA" ~ " CHROMATIN SILENCING\nAT RDNA",
        pathway == " GLYCOPROTEIN BIOSYNTHETIC PROCESS" ~ " GLYCOPROTEIN\nBIOSYNTHETIC PROCESS",
        pathway == " GLYCOPROTEIN METABOLIC PROCESS" ~ " GLYCOPROTEIN\nMETABOLIC PROCESS",
        pathway == " EPITHELIAL MESENCHYMAL TRANSITION" ~ " EPITHELIAL MESENCHYMAL\nTRANSITION",
        pathway == " INTERFERON ALPHA RESPONSE" ~ " INTERFERON ALPHA\nRESPONSE",
        pathway == " INTERFERON GAMMA RESPONSE" ~ " INTERFERON GAMMA\nRESPONSE",
        pathway == " TNFA SIGNALING VIA NFKB" ~ "TNFA SIGNALING\nVIA NFKB",
        pathway == " MITOCHONDRIAL GENE EXPRESSION" ~ " MITOCHONDRIAL GENE\nEXPRESSION",
        pathway == " PROTEIN TARGETING TO MEMBRANE" ~ " PROTEIN TARGETING\nTO MEMBRANE",
        pathway == " PROTEIN HETEROOLIGOMERIZATION" ~ " PROTEIN\nHETEROOLIGOMERIZATION",
        pathway == " REGULATION OF GENE SILENCING" ~ " REGULATION OF\nGENE SILENCING",
        pathway == " NUCLEOSOME ORGANISATION" ~ " NUCLEOSOME\nORGANISATION",
        #pathway == " CHROMATIN ASSEMBLY" ~ " CHROMATIN\nASSEMBLY",
        pathway == " CHROMATIN ASSEMBLY OR DISASSEMBLY" ~ " CHROMATIN ASSEMBLY\nOR DISASSEMBLY",
        pathway == " NUCLEOSOME ORGANIZATION" ~ " NUCLEOSOME\nORGANIZATION",
        pathway == " TRANSLATION INITIATION" ~ " TRANSLATION\nINITIATION",
        pathway == " REGULATION OF EXTRINSIC APOPTOTIC SIGNALING PATHWAY IN ABSENCE OF LIGAND" ~ " REGULATION OF EXTRINSIC\nAPOPTOTIC SIGNALING PATHWAY\nIN ABSENCE OF LIGAND",
        pathway == " NEGATIVE REGULATION OF PROTEIN POLYMERIZATION" ~ " NEGATIVE REGULATION OF\nPROTEIN POLYMERIZATION",
        pathway == " NUCLEAR CHROMOSOME TELOMERIC REGION" ~ " NUCLEAR CHROMOSOME\nTELOMERIC REGION",
        .default = pathway))
  
  plot <- ggplot(data = all2, aes(reorder(pathway, NES), NES)) +
    geom_col(fill = "#984EA3") +
    coord_flip() +
    labs(x="Pathway", y="NES",
         title=title) + 
    publication_theme()+
    theme(axis.title.y = element_blank())
  return(plot)
}







#read in DESeq2 output----
DESeq2_data <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output/FINAL_sieS25_merged_DESeq2.csv"))



#make named vectors----
### This will make a list of named vectors to allow gsea to be carried out separately on the RPFs, Totals and TE log2FC
DESeq2_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(RPFs_log2FC)) %>%
  deframe() -> RPFs_named_vector

DESeq2_data %>%
  group_by(gene_sym) %>%
  #dplyr::rename(totals_log2FC = log2FoldChange) %>% 
  summarise(stat = mean(totals_log2FC)) %>%
  deframe() -> totals_named_vector

DESeq2_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(TE_log2FC)) %>%
  deframe() -> TE_named_vector

named_vectors <- list(RPFs = RPFs_named_vector,
                      totals = totals_named_vector,
                      TE = TE_named_vector)

save_path = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/gsea/sieS25"

#read in pathways----
source("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/read_human_GSEA_pathways.R") #This may need to be changed to human



#hallmark----
#carry out fgsea
hallmark_results <- lapply(named_vectors, run_fgsea, pathway = pathways.hallmark)

#save results
#save(file = file.path(parent_dir, "Analysis/fgsea/hallmark_results.Rdata"), hallmark_results)

#padj <- 0.05

#plot enriched pathways
#png(filename = file.path(save_path, paste(treatment, "RPFs_hallmark.png", sep = "_")), width = 2000, height = 1500, res = 300)
hallmarks_RPFs <- make_plot2(fgsea_result = hallmark_results$RPFs, padj_threshold = padj, title = "Hallmark")
#dev.off()

#png(filename = file.path(save_path, paste(treatment, "totals_hallmark.png", sep = "_")), width = 2000, height = 1500, res = 300)
hallmark_totals <- make_plot2(fgsea_result = hallmark_results$totals, padj_threshold = padj, title = "Hallmark")
#dev.off()

#png(filename = file.path(save_path, paste(treatment, "TE_hallmark.png", sep = "_")), width = 2000, height = 1500, res = 300)
hallmark_TE <- make_plot2(fgsea_result = hallmark_results$TE, padj_threshold = padj, title = "Hallmark")
#dev.off()

#biological processes----
#carry out fgsea
bio_processes_results <- lapply(named_vectors, run_fgsea, pathway = pathways.bio_processes)

#save results
#save(file = file.path(parent_dir, "Analysis/fgsea/bio_processes_results.Rdata"), bio_processes_results)

#set adjusted p-value
#padj <- 0.05

#plot enriched pathways
#png(filename = file.path(save_path, paste(treatment, "RPFs_bio_processes.png", sep = "_")), res = 300, width = 3000, height = 1500)
bio_process_RPFs <- make_plot2(fgsea_result = bio_processes_results$RPFs, padj_threshold = padj, title = " Biological Process")
#dev.off()

#png(filename = file.path(save_path, paste(treatment, "totals_bio_processes.png", sep = "_")), res = 300, width = 3500, height = 1500)
#bio_process_totals <- make_plot2(fgsea_result = bio_processes_results$totals, padj_threshold = padj, title = "Biological Process")
#dev.off()

#png(filename = file.path(save_path, paste(treatment, "TE_bio_processes.png", sep = "_")), res = 300, width = 3750, height = 1500)
bio_process_TE <- make_plot2(fgsea_result = bio_processes_results$TE, padj_threshold = padj, title = "Biological process")
#dev.off()

#molecular functions----
#carry out fgsea
mol_funs_results <- lapply(named_vectors, run_fgsea, pathway = pathways.mol_funs)

#save results
#save(file = file.path(parent_dir, "Analysis/fgsea/mol_funs_results.Rdata"), mol_funs_results)

#set adjusted p-value
#padj <- 0.05

#plot enriched pathways
#png(filename = file.path(save_path, paste(treatment, "RPFs_mol_funs.png", sep = "_")), width = 1000, height = 1000)
make_plot2(fgsea_result = mol_funs_results$RPFs, padj_threshold = padj, title = "Molecular function")
#dev.off()

#png(filename = file.path(save_path, paste(treatment, "totals_mol_funs.png", sep = "_")), width = 1000, height = 1000)
#make_plot2(fgsea_result = mol_funs_results$totals, padj_threshold = padj, title = "Molecular function")
#dev.off()

#png(filename = file.path(save_path, paste(treatment, "TE_mol_funs.png", sep = "_")), width = 2000, height = 800)
make_plot2(fgsea_result = mol_funs_results$TE, padj_threshold = padj, title = "Molecular function")
#dev.off()



#cellular component----
#carry out fgsea
#cell_comp_results <- lapply(named_vectors, run_fgsea, pathway = pathways.cell_comp)

#save results
#save(file = file.path(parent_dir, "Analysis/fgsea/cell_comp_results.Rdata"), cell_comp_results)

#set adjusted p-value
#padj <- 0.05

#plot enriched pathways
#png(filename = file.path(save_path, paste(treatment, "RPFs_cell_comp.png", sep = "_")), res = 300, width = 2500, height = 1500)
cellular_component_RPFs <- make_plot2(fgsea_result = cell_comp_results$RPFs, padj_threshold = padj, title = "Cellular component")
#dev.off()

#png(filename = file.path(save_path, paste(treatment, "totals_cell_comp.png", sep = "_")), res = 300, width = 2500, height = 1500)
#cellular_component_totals <- make_plot2(fgsea_result = cell_comp_results$totals, padj_threshold = padj, title = "Cellular component")
#dev.off()

#png(filename = file.path(save_path, paste(treatment, "TE_cell_comp.png", sep = "_")), res = 300, width = 2750, height = 1500)
cellular_component_TE <- make_plot2(fgsea_result = cell_comp_results$TE, padj_threshold = padj, title = "Cellular component")
#dev.off()

# RPFs
png(file.path(save_path, "sieS25_GSEA.png"), res = 300, height = 5000, width = 4500)
cowplot::plot_grid(cellular_component_RPFs,
                   cellular_component_TE,
                   bio_process_RPFs,
                   bio_process_TE,
                   #hallmarks_RPFs,
                   align = "v", 
                   nrow = 2, 
                   ncol = 2,
                   #rel_heights = c(0.8,1,0.8),
                   axis = "l")
dev.off()




