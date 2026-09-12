#load libraries
library(tidyverse)
library(readxl)

source('\\\\data.beatson.gla.ac.uk/data/R11/external_sequencing_data/Gillen_2021_halflife/Scripts/R_scripts/common_variables.R')

#makes a label from the output of a correlation test to include r and p values
myR <- function(x) {
  r <- round(as.numeric(x$estimate), digits = 2)
  p <- as.numeric(x$p.value)
  if (p > 0.001) {
    rounded_p <- round(p, digits = 3)
    p_label <- paste('P =', rounded_p)
  }else{
    if (p<2.2E-16) {
      p_label <- 'P < 2.2e-16'
    } else {
      rounded_p <- formatC(p, format = "e", digits = 2)
      p_label <- paste('P =', rounded_p)
    }
  }
  return(paste0('r = ', r, '\n', p_label))
}

#read in rscu data and remove '.' extensions to IDs
load("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_RSCU.Rdata")
rscu_IDs <- names(RSCU_list)
rscu_IDs <- str_remove(rscu_IDs, "\\..+")
names(RSCU_list) <- rscu_IDs

#read in Sarah's HEK293 most abundant transcripts and remove '.' extensions to IDs
most_abundant_transcripts <- read.table("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/most_abundant_transcripts/Gillen_HEK293/most_abundant_transcripts.txt")
most_abundant_transcripts <- str_remove(most_abundant_transcripts$V1, "\\..+")

#read in transcript to gene IDs and filter to most abundant transcripts
gene_to_transcript_ID <- read_csv("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_gene_IDs.csv", col_names = c("ENST", "ENSG", "Gene_ID"))
gene_to_transcript_ID %>%
  mutate(ENST = str_remove(ENST, "\\..+")) %>%
  select(-ENSG) %>%
  filter(ENST %in% most_abundant_transcripts) -> most_abundant_transcript_IDs


#make a data frame of rscu data for just the most abundant transcripts
table(most_abundant_transcript_IDs$ENST %in% names(RSCU_list))
rscu_data <- do.call("rbind", RSCU_list[most_abundant_transcript_IDs$ENST])
rscu_data$transcript <- str_remove(rscu_data$transcript, "\\..+")
n_distinct(rscu_data$transcript)
n_distinct(rscu_data$gene)

#read in Sarah's CNOT1 half-lives data and merge with rscu data
CNOT_data <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/R11/SarahData/CNOT1/mRNA_half_life/CNOTdata_mRNAhalflives.csv")

CNOT_data %>%
  inner_join(most_abundant_transcript_IDs, by = "Gene_ID") %>%
  inner_join(rscu_data, by = c("Gene_ID" = "gene_sym", "ENST" = "transcript")) -> merged_data

summary(merged_data)
n_distinct(merged_data$Gene_ID)
n_distinct(merged_data$ENST)

#calculate CSCs and CNOT-sep CSCs
codons <- as.character(unique(merged_data$codon))
codons <- codons[!codons %in% c("UAA", "UAG", "UGA", "AUG", "UGG")]

features <- 
  read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/features.csv", col_names = T) %>% 
  select(Gene, ENST, signal_sequence) %>% 
  mutate(ENST = str_remove(ENST, "\\..+"))

merged_data <- 
  merged_data %>% 
  inner_join(features, by = c("Gene_ID" = "Gene", "ENST"))

merged_data_split <- split(merged_data, merged_data$signal_sequence)

cor_list <- list()

for (codon in codons) {
  df <- merged_data[merged_data$codon == codon,]
  
  ctrl_r <- cor.test(df$freq, df$siControl_half_life, method = "spearman")
  delta_r <- cor.test(df$freq, df$log2FC_halflife, method = "spearman")
  
  cor_list[[codon]] <- data.frame(codon = codon, ctrl_cor_estimate = ctrl_r$estimate, ctrl_p_value = ctrl_r$p.value,
                                            delta_cor_estimate = delta_r$estimate, delta_p_value = delta_r$p.value)
  
  
}

CSCs <- do.call("rbind", cor_list) %>% 
  mutate(codon = str_replace_all(tolower(codon), pattern = "u", replacement = "t"))

rownames(CSCs) <- NULL

write.csv(CSCs, "C:/Users/Jettles/OneDrive - University of Glasgow/Documents/get_shit_done/chapter5/SGs_CSC_scores.csv")

#plot Ctrl CSC
CSCs %>%
  arrange(ctrl_cor_estimate) %>%
  pull(codon) -> ctrl_ordered_codons

CSCs %>%
  arrange(delta_cor_estimate) %>%
  pull(codon) -> delta_ordered_codons

# read in publication theme
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
            legend.text = element_text(size = 12),
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


CSCs %>%
  mutate(wobble = case_when(str_sub(codon, 3,3) %in% c("A", "U") ~ "A/U ending codons",
                   .default = "G/C ending codons")) %>% 
  ggplot(aes(x = factor(codon, levels = ctrl_ordered_codons, ordered = T), y = ctrl_cor_estimate, fill = wobble))+
  geom_col()+
  scale_fill_manual(values = c("#BD33A4","#33c7d8")) +
  #scale_fill_brewer(palette = "Accent") +
  #xlab("codon")+
  #ylab("CSC")+
  #ggtitle("Cytosolic transccripts") +
  publication_theme()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 11),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.direction = "vertical") -> ctrl_plot


png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/SG_CSCs_ctrl.png", height = 1250, width = 3000, res = 300)
print(ctrl_plot)
dev.off()

"#00FFFF"
"#BD33A4"

CSCs %>%
  mutate(wobble = case_when(str_sub(codon, 3,3) %in% c("A", "U") ~ "A/U ending codons",
                            .default = "G/C ending codons")) %>% 
  ggplot(aes(x = factor(codon, levels = ctrl_ordered_codons, ordered = T), y = ctrl_cor_estimate, fill = wobble))+
  geom_col()+
  scale_fill_manual(values = c("#BD33A4","#33c7d8")) +
  #scale_fill_brewer(palette = "Accent") +
  #xlab("codon")+
  #ylab("CSC")+
  #ggtitle("Cytosolic transccripts") +
  publication_theme()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 11),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.direction = "vertical") -> ctrl_plot


png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/SG_CSCs_ctrl.png", height = 1250, width = 3000, res = 300)
print(ctrl_plot)
dev.off()


CSCs %>%
  mutate(wobble = case_when(str_sub(codon, 3,3) %in% c("A", "U") ~ "A/U ending codons",
                            .default = "G/C ending codons")) %>% 
  ggplot(aes(x = factor(codon, levels = delta_ordered_codons, ordered = T), y = delta_cor_estimate , fill = wobble))+
  geom_col()+
  scale_fill_manual(values = c("#BD33A4","#33c7d8")) +
  #scale_fill_brewer(palette = "Accent") +
  #xlab("codon")+
  #ylab("CSC")+
  #ggtitle("Cytosolic transccripts") +
  publication_theme()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 11),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.direction = "vertical") -> delta_plot

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/SG_CSCs_delta.png", height = 1250, width = 3000, res = 300)
print(delta_plot)
dev.off()
