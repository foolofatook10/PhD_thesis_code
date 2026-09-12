library(tidyverse)
library(data.table)

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
            axis.text.y = element_text(size = 14, color = "black"),
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

home = "\\\\data.beatson.gla.ac.uk/data"

path_to_geneIDs = file.path("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/features.csv")

gene_IDs <- 
  read_csv(path_to_geneIDs) %>% select(Gene, ENSG, ENST)


path_to_siCNOT3 = file.path("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/riboseq_datasets/RPF_nCPMs_siCNOT3.csv")

siCNOT3_data <- fread(path_to_siCNOT3)

mean_siCNOT3 <- 
siCNOT3_data %>% 
  group_by(KD, transcript, codon) %>% 
  summarise(mean_nCPM = mean(normalised_CPM)) %>% 
  inner_join(gene_IDs, by = c("transcript" = "ENST"))

siCNOT3_deltas <- 
siCNOT3_data %>% 
  select(-mean_tpm, -CPM) %>% 
  spread(key = KD, value = normalised_CPM) %>% 
  na.omit() %>% 
  mutate(delta = siCNOT3 - NTC) %>% 
  select(-NTC, -siCNOT3) %>% 
  group_by(transcript, codon) %>% 
  summarise(mean_nCPM_delta = mean(delta))  %>% 
  inner_join(gene_IDs, by = c("transcript" = "ENST")) 


codons_deseq <- 
  read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/Analysis/DESeq2_output/M_RPFs_28nts.csv") %>% 
  filter(log2FoldChange > 0 & padj <0.05)

codons_deseq %>% filter(gene_sym == "RPS9")

deltas_filt <- 
  siCNOT3_deltas %>% 
  filter(Gene %in% unique(codons_deseq$gene_sym))

deltas_filt_split <- split(deltas_filt, deltas_filt$Gene)

deltas_filt_split[[1]]

plot_CPM_gplots <- function(x){
  
  df <- deltas_filt_split[[x]]
  gene_id <- unique(df$Gene)
  
  enriched_position <- 
    codons_deseq %>% 
    filter(gene_sym == {{gene_id}}) %>% 
    pull(codon)
  
  final <- 
    df %>% ggplot(aes(x = codon, y = mean_nCPM_delta))+
    geom_line() +
    geom_vline(xintercept = enriched_position, alpha = 0.5, lty = "dashed", colour = "#0C7BDC") +
    ylab("Delta (siCNOT3 - Total)") +
    ggtitle(paste0(unique(df$Gene))) +
    publication_theme()
  
  return(final)
  
}

gplots <- lapply(1:length(deltas_filt_split), plot_CPM_gplots)
names(gplots) <- unlist(lapply(deltas_filt_split, function(x){x %>% pull(Gene) %>% unique()}))

ribosomal_proteins <- names(gplots)[str_detect(names(gplots), regex("^RP[LS]\\d+"))]

save_path = "C:/Users/Jettles/OneDrive - University of Glasgow/Documents/siCNOT3_with_stalls_plotted"

for(i in names(gplots)){
  
  png(file = file.path(save_path, "all", paste0(i, ".png")), res = 300, height = 1000, width = 1250)
  print(gplots[[i]])
  dev.off()
}

for(i in ribosomal_proteins){
  
  png(file = file.path(save_path, paste0(i, ".png")), res = 300, height = 1000, width = 1250)
  print(gplots[[i]])
  dev.off()
}


gplots$ATF4
gplots$EEF2

gplots$RPL10
gplots$RPL10A
gplots$RPL12
gplots$RPL7A

gplots$RPS14

enriched_positions_MYC <- codons_deseq %>% filter(gene_sym == "MYC") %>% pull(codon)

siCNOT3_deltas %>% 
  inner_join(gene_IDs, by = c("transcript" = "ENST")) %>% 
  filter(Gene == "MYC") %>% 
  ggplot(aes(x = codon, y = mean_nCPM_delta)) +
  geom_line() +
  geom_vline(xintercept = enriched_positions_MYC, alpha = 0.5, lty = "dashed")


enriched_positions_CDC5L <- codons_deseq %>% filter(gene_sym == "CDC5L") %>% pull(codon)

siCNOT3_deltas %>% 
  inner_join(gene_IDs, by = c("transcript" = "ENST")) %>% 
  filter(Gene == "CDC5L") %>% 
  ggplot(aes(x = codon, y = mean_nCPM_delta)) +
  geom_line() +
  geom_vline(xintercept = enriched_positions_CDC5L, alpha = 0.5, lty = "dashed")





DSOBs <- read_csv ("")

path_to_CNOT1_abundant_transcripts = file.path(home, "R11/James/CNOT3_Riboseq/Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv")
path_to_CNOT3_deseq = file.path(home, "R11/James/CNOT3_Riboseq/Analysis/DESeq2_output/merged_DESeq2.csv")