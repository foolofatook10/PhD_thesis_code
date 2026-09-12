library(tidyverse)
library(data.table)

# read in publication theme
publication_theme <- function(base_size=18, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = 22, hjust = 0.5),
            #text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = 20),
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
  read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/Analysis/DESeq2_output/M_RPFs_28nts_normalised_counts_M_28nts_without_sustained_CNOT3_binding_transcripts.csv") 

master <- 
  read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/master.csv") %>% 
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
  fread("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/CPMs/CDS_CPMs.csv") %>% 
  filter(MD30 == "M") %>% 
  mutate(CPM = (summed_CPM_codon/3)) %>% 
  group_by(IP, transcript, codon) %>% 
  summarise(mean_CPM = mean(CPM))

CPMs_transcript <- 
  CPMs %>% 
  inner_join(master %>% select(Gene, ENST),
             by = c("transcript" = "ENST"))

# Load in TOP mRNAs
TOP_genes <- read_csv("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/Gene_lists/TOPs/pnas.1912864117.sd01.csv", col_names = F) %>% pull(X1)

TOP_CPMs <- 
CPMs_transcript %>% 
  filter(Gene %in% TOP_genes)

TOP_CPMs_split <- split(TOP_CPMs, TOP_CPMs$transcript)

# plot_CPM_gplots <- function(x){
#   
#   df <- CPMs_transcript_filtered_split[[x]]
#   ENST <- unique(df$transcript)
#   
#   df$IP <- factor(df$IP, levels = c("CNOT3", "Tot"))
#   
#   final <- 
#     df %>% 
#     ggplot(aes(x = codon, y = mean_CPM, colour = IP))+
#     geom_line(size = 1) +
#     scale_color_manual(values = c("#0C7BDC", "#FFC20A")) +
#     xlab("Codon position") +
#     ylab("CPM") +
#     #geom_vline(xintercept = enriched_position, alpha = 0.5) +
#     ggtitle(paste0(unique(df$Gene))) +
#     publication_theme() +
#     theme(legend.position = "none")
#   
#   return(final)
#   
# }

plot_CPM_gplots <- function(x){
  
  df <- TOP_CPMs_split[[x]]
  ENST <- unique(df$transcript)
  
  enriched_position <- 
    enriched %>% 
    filter(transcript == {{ENST}}) %>% 
    pull(codon)
  
  df$IP <- factor(df$IP, levels = c("Tot", "CNOT3"))
  
  final <- 
    df %>% ggplot(aes(x = codon, y = mean_CPM, colour = IP))+
    geom_line(size = 1) +
    scale_color_manual(values = c("#FFC20A", "#0C7BDC")) +
    geom_vline(xintercept = enriched_position, alpha = 0.5, lty = "dashed") +
    xlab("Codon position") +
    ylab("Mean CPM") +
    ggtitle(paste0(unique(df$Gene))) +
    publication_theme() +
    theme(legend.position = "none")
  
  return(final)
  
}

gplots <- lapply(1:length(TOP_CPMs_split), plot_CPM_gplots)
names(gplots) <- unlist(lapply(TOP_CPMs_split, function(x){x %>% pull(Gene) %>% unique()}))

save_plots <- function(x){
  
  gene_name <- names(gplots)[[x]]
  
  thesis_figures_dir <- "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/TOP_mRNAs"
  png(filename = file.path(thesis_figures_dir,paste0(gene_name, ".png")),
      res = 300,
      height = 1000,
      width = 1500)
  print(gplots[[x]])
  dev.off()
  
}

lapply(1:length(gplots), save_plots)

TOP_mRNA_selection <- 
cowplot::plot_grid(
gplots$RPS4,
gplots$RPS17,
gplots$RPL3,
gplots$RPL17,
gplots$RPL35,
gplots$RPL18,
gplots$RPL17,
gplots$RPL12,
ncol = 2,
nrow = 4)


png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/TOP_mRNAs/TOP_mRNA_selection.png",
    res = 300, height = 3000, width = 2000)
print(TOP_mRNA_selection)
dev.off()
