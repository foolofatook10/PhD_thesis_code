library(tidyverse)
library(data.table)

# read in publication theme
publication_theme <- function(base_size=18, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            #text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            #axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2, size = 20),
            axis.title.x = element_text(vjust = -0.2, size = 20),
            axis.text.x = element_text(size = 18, color = "black"),
            axis.text.y = element_text(size = 18, color = "black"),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = 12),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", size = 14),
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

codons_deseq <- 
  read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/Analysis/DESeq2_output/M_RPFs_28nts.csv") 

master <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/master.csv") %>% 
  select(Gene, ENSG, ENST, nucleotide_sequence_CDS)

# make function split_into_codons
split_into_codons <- function(string) {
  
  # total length of string
  num.chars <- nchar(string)
  
  # the indices where each substr will start
  starts <- seq(1,num.chars, by=3)
  
  # chop it up
  sapply(starts, function(ii) {
    substr(string, ii, ii+2)
  })
}


codons_deseq <- 
  codons_deseq %>% filter(transcript %in% c(master$ENST))

enriched <- codons_deseq %>% filter(log2FoldChange > 0 & padj < 0.05)

CPMs <- 
  fread("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/CDS_CPMs.csv") %>% 
  filter(MD30 == "M") %>% 
  mutate(CPM = (summed_CPM_codon/3)) %>% 
  group_by(IP, transcript, codon) %>% 
  summarise(mean_CPM = mean(CPM))

CPMs_transcript <- 
  CPMs %>% 
  inner_join(master %>% select(Gene, ENST),
             by = c("transcript" = "ENST"))

gene_syms <- c("TP53", "MYC", "CDKN1A", "CDKN1B", "CDKN1C", "CDKN1D",
               "CDKN2A", "CDKN2B", "CDKN2C", "CDKN2D", 
               "CDK1", "CDK2", "CDK3", "CDK4", "CCNB2",  
               "RIPK1", "BAM", "BAX", "AURKB", "MAD1L1",
               "TSC1", "ATG7", "PIK3C3", "SQSTM1", "DDIT3", "EGFR",
               "MAPK1", "JUN", "PRPF31", "TADA2A", "TADA2B")

CPMs_transcript_filtered <- 
CPMs_transcript %>% 
  filter(Gene %in% gene_syms)

CPMs_transcript_filtered_split <- 
split(CPMs_transcript_filtered, CPMs_transcript_filtered$Gene)

CPMs_transcript_filtered_split[[1]]$IP

plot_CPM_gplots <- function(x){
  
  df <- CPMs_transcript_filtered_split[[x]]
  ENST <- unique(df$transcript)
  
  df$IP <- factor(df$IP, levels = c("CNOT3", "Tot"))
  
  final <- 
    df %>% 
    ggplot(aes(x = codon, y = mean_CPM, colour = IP))+
    geom_line(size = 1) +
    scale_color_manual(values = c("#0C7BDC", "#FFC20A")) +
    xlab("Codon position") +
    ylab("CPM") +
    #geom_vline(xintercept = enriched_position, alpha = 0.5) +
    ggtitle(paste0(unique(df$Gene))) +
    publication_theme() +
    theme(legend.position = "none")
  
  return(final)
  
}

genes <- lapply(1:length(CPMs_transcript_filtered_split), plot_CPM_gplots)
gene_names <- unlist(lapply(1:length(CPMs_transcript_filtered_split), 
                     function(x){
                       df <- CPMs_transcript_filtered_split[[x]]
                       gene_name <- df %>% pull(Gene) %>% unique()
                       return(gene_name)
                       }))

names(genes) <- gene_names

genes[["CDKN2A"]]

path_to_thesis_figures <- "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures"

new_gene_syms <- c("CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B", "CDKN2C", "TP53", "CCNB2", "PRPF31", "JUN", "RIPK1")

for(i in 1:length(new_gene_syms)){
  
  name = gene_syms[[i]]
  png(filename = file.path(path_to_thesis_figures, paste0(name, ".png")),
      res = 300, height = 1000, width = 1500)
  print(genes[[name]])
  dev.off()
  
}

png(filename = file.path(path_to_thesis_figures, "genes.png"),
    res = 300, height = 3500, width = 2500)
print(cowplot::plot_grid(genes[["CDKN1A"]], genes[["CDKN1B"]],
                   genes[["CDKN2A"]], genes[["CDKN2B"]],
                   genes[["CDKN2C"]], genes[["TP53"]],
                   genes[["CCNB2"]], genes[["PRPF31"]],
                   genes[["JUN"]], genes[["RIPK1"]],
                   ncol = 2, align = "v"))
dev.off()

png(file.path(path_to_thesis_figures))






CPMs_transcript <- 
  CPMs_transcript %>% filter(transcript %in% unique(enriched$transcript))

splitted <- 
  split(CPMs_transcript, CPMs_transcript$transcript)

names(splitted) <- unique(CPMs_transcript$Gene)

plot_CPM_gplots <- function(x){
  
  df <- splitted[[x]]
  ENST <- unique(df$transcript)
  
  enriched_position <- 
    enriched %>% 
    filter(transcript == {{ENST}}) %>% 
    pull(codon)
  
  final <- 
    df %>% ggplot(aes(x = codon, y = mean_CPM, colour = IP))+
    geom_line() +
    geom_vline(xintercept = enriched_position, alpha = 0.5) +
    ggtitle(paste0(unique(df$Gene))) +
    publication_theme()
  
  return(final)
  
}

gplots <- lapply(1:length(splitted), plot_CPM_gplots)

names(gplots) <- unique(CPMs_transcript$Gene)

unique(CPMs_transcript$Gene)



gplots$SQSTM1

gplots$PRPF40A
gplots$PRPF8
gplots$PRPF6

gplots$TAF10

gplots$CAPRIN1
gplots$EIF2AK1
gplots$TRIM28
gplots$FLRT2
gplots$MLF2
gplots$TRAM2
gplots$AAAS
gplots$NOTCH2
gplots$EEF1B2
gplots$EEF2
gplots$ATF4

gplots$RCN1
gplots$WNK1

gplots$ENO1

gplots$COPE
gplots$IFNGR1
gplots$ITGB1
gplots$G3BP1
gplots$AARS1
gplots$GPI
gplots$SELENOS
gplots$SELENOT

gplots$TMEM30A
gplots$SNAI2
gplots$CHPF2
gplots$GRN
gplots$RCN1
gplots$ADCK2
gplots$CDIPT
gplots$NOMO1
gplots$NDUFB7
gplots$HSP90AA1
gplots$SLC7A5
gplots$EXOSC5

gplots$NCLN
gplots$PSENEN
gplots$CCDC47
gplots$TMCO1
gplots$TRAM1
gplots$TRAM2

gplots$TMEM109
gplots$GORASP2
gplots$SMAD3

gplots$CDKN1A


