library(tidyverse)
library(data.table)
library(readxl)



### Part 1: calculating weighted codon usage in HEK293s ###

# read in TPMs of NTC HEK293 cells
CNOT3_riboseq_tpms <- 
read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/Analysis/DESeq2_output/tpms.csv")

# make data into a long format
CNOT3_riboseq_tpms_long <- 
CNOT3_riboseq_tpms %>% 
  gather(key = "sample", value = "tpm", REP1_NTC_Tot:REP5_siRPS25_Tot)

# clean up data
terms <- str_split(CNOT3_riboseq_tpms_long$sample, "_")

CNOT3_riboseq_tpms_long$replicate <- as.numeric(str_remove_all(unlist(lapply(terms, function(x){x[[1]]})), "REP"))
CNOT3_riboseq_tpms_long$condition <-  unlist(lapply(terms, function(x){x[[2]]}))

# only interested in NTC TPMs, calculate mean TPMs across all replicates
NTC_tpms <- 
  CNOT3_riboseq_tpms_long %>% 
  filter(condition == "NTC") %>% 
  select(-sample) %>% 
  group_by(condition, transcript) %>% 
  summarise(mean_tpm = mean(tpm)) %>% 
  ungroup() %>% select(-condition)

# read in master table
master <- as_tibble(fread(file = "\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/master.csv", header = T, drop = "V1"))

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

# make a list with the nucleotide sequences chopped into codons
my_lyst <- sapply(master$nucleotide_sequence_CDS, split_into_codons)

# Name each element of the list by its ENST tag
names(my_lyst) <- pull(master, var = ENST)

# count the codons per transcript
enumerated_codons <- lapply(lapply(my_lyst, table), data.frame)
transcripts <- names(my_lyst)

# make a function to attach 
attach_transcript <- 
  function(x){
    
    df <- enumerated_codons[[x]]
    ENST <- transcripts[[x]]
    
    df_final <- df %>% mutate(transcript = ENST)
    
    return(df_final)
    
  }

# raw table of codons per transcript
codons_counted <- 
  do.call("rbind", lapply(1:length(transcripts), attach_transcript))

colnames(codons_counted) <- c("codon", "raw_number", "transcript")

# calculate weighted codon usage by multiplying the frequency of codon per transcript with transcript's TPM
wcus <- 
  codons_counted %>% 
  group_by(transcript) %>% 
  summarise(total_residues = sum(raw_number)) %>% 
  inner_join(codons_counted, by = "transcript") %>% 
  mutate(freq = raw_number/total_residues * 100) %>% 
  select(transcript, codon, freq) %>% 
  inner_join(NTC_tpms, by = "transcript") %>% 
  mutate(wcu = freq * mean_tpm)

# calculate mean weighted codon usage across transcriptome
mean_wcus <- 
  wcus %>% 
  group_by(codon) %>% 
  summarise(mean_wcu = mean(wcu)) %>% 
  mutate(codon = gsub("T", "U", toupper(codon)))

# import table with codon-anticodon relationships
codon_anticodon <- 
  read_csv("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/useful_tables/codon_anticodon_pairs.csv") %>% 
  filter(Type == "Cognate")

codons_all <- inner_join(codon_anticodon, mean_wcus, by = c("codon")) 

### end of part 1 ###

### Part 2: exploring tRNA anticodon availability ###

data <- 
read_excel(path = "\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/Behrens_et_al.xlsx",
           sheet = 2,
           skip = 1,
           col_names = T)

data <- 
data %>% 
  select(Anticodon, `iPSC rep 1`, `iPSC rep 2`, `HEK293T rep 1`, `HEK293T rep 2`, `K562 rep 1`, `K562 rep 2`) %>% 
  mutate(anticodon = str_remove(str_extract(Anticodon, "-.*"),"-")) %>% 
  rename(AA = Anticodon) %>% 
  mutate(AA = str_remove(AA, "-.*")) %>% 
  select(anticodon, AA, `iPSC rep 1`:`K562 rep 2`) %>% 
  gather(key = sample, value = norm_counts, `iPSC rep 1`: `K562 rep 2`) %>% 
  mutate(rep = as.numeric(str_remove(str_extract(sample, "rep \\d"), "rep "))) %>% 
  mutate(sample = str_remove(sample, " rep \\d")) %>% 
  group_by(sample, AA, anticodon) %>% 
  summarise(norm_counts = mean(norm_counts))

data$anticodon <- gsub("T","U", data$anticodon)



 



data_split <- split(data, f = data$sample)

data_split_freq <- 
lapply(data_split, inner_join, codons_all, by = c("anticodon", "AA"))

attach_leucyl <- function(x){
  
  y <- x %>% 
    mutate(leucyl = case_when(codon == "CUU" ~ "CUU",
                              codon == "CUA" ~"CUA",
                              codon == "CUG" ~ "CUG",
                              codon == "UUA" ~ "UUA",
                              codon == "UUG" ~"UUG",
                              codon == "CUC" ~"CUC",
                              .default = "other"))
  
  return(y)
}

apply_factor <- function(x){
  
  x$leucyl <- factor(x$leucyl, levels = c("other", "CUG", "CUC", "CUU", "UUG", "CUA", "UUA"))
  return(x)
  
}

data_split_freq <- lapply(lapply(data_split_freq, attach_leucyl), apply_factor)

# read in publication theme
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

plot_tRNAs <- function(x){
  
  ggtit <- x$sample %>% unique()
  
  gplot <- 
    x %>% 
      ggplot(aes(x = norm_counts, y = mean_wcu, colour = leucyl)) +
      geom_point(size = 2) +
      geom_smooth(method='lm', formula= y~x)+
      geom_abline(slope = 0, intercept = 0, linewidth = 1, linetype = "dashed") +
      scale_color_manual(values = c("#C8C8C8", "#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F")) +
      ggtitle(paste(ggtit)) +
      publication_theme() +
      xlab("Primary anticodon abundance") +
      ylab("Weighted Codon Usage") +
      theme(legend.direction = "vertical", 
            legend.position = "none")
  
  return(gplot)

}


gplots <- lapply(data_split_freq, plot_tRNAs)

data_split_freq$HEK293T %>% 
  filter(mean_wcu > 100)
data_split_freq$HEK293T %>% 
  filter(codon == "GAG")

label_data <- data_split_freq$HEK293T %>% 
  filter(mean_wcu > 50) %>% 
  filter(codon != "CUG") %>% 
  filter(codon != "CCC") %>% 
  filter(codon != "GCC")

data_split_freq$HEK293T %>% arrange(norm_counts)
data_split_freq$HEK293T %>% filter(codon == "AUC")
data_split_freq$HEK293T %>% filter(codon == "AAG")


wcu_leucyl <- data_split_freq$HEK293T %>% filter(leucyl != "other")
wcu_lysine <- data_split_freq$HEK293T %>% filter(codon == "AAG")
wcu_ile <- data_split_freq$HEK293T %>% filter(codon == "AUC")
wcu_other <-  data_split_freq$HEK293T %>% filter(leucyl == "other")

wcu_leucyl$leucyl <- factor(wcu_leucyl$leucyl, levels = c("CUG", "CUC", "CUU", "UUG", "CUA", "UUA"),
                            labels = c("CTG", "CTC", "CTT", "TTG", "CTA", "TTA")) 


  
HEK_updated <- 
ggplot() +
  geom_smooth(data = data_split_freq$HEK293T, 
              aes(x = norm_counts, y = mean_wcu),
              method='lm', 
              formula= y~x, 
              color = "black", 
              se = F, 
              lty = "longdash") +
  ggpubr::stat_cor(data = data_split_freq$HEK293T,
                   aes(x = norm_counts, y = mean_wcu),
                   method = "pearson", size = 5) +  
  geom_point(data = wcu_other, 
             aes(x = norm_counts, y = mean_wcu), 
             size = 2, colour = "#C8c8c8") +
    geom_point(data = wcu_leucyl, 
               aes(x = norm_counts, y = mean_wcu, colour = leucyl), 
               size = 2) +
    ggrepel::geom_text_repel(data = wcu_leucyl, 
                             aes(x = norm_counts, y = mean_wcu, colour = leucyl, label = leucyl)) +
  #ggrepel::geom_text_repel(data = wcu_lysine, 
   #                        aes(x = norm_counts, y = mean_wcu, label = codon), colour = "black") +
  ggrepel::geom_text_repel(data = wcu_ile, 
                           aes(x = norm_counts, y = mean_wcu, label = codon), colour = "black") +
  ggrepel::geom_text_repel(data = label_data, aes(x = norm_counts, y = mean_wcu, label = codon), 
                           colour = "black") +
  scale_color_manual(values = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F")) +
  
    #geom_abline(slope = 0.001, intercept = 0, linewidth = 1, linetype = "dashed") +
  #scale_color_manual(values = c("#C8C8C8", "#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F")) +
  ggtitle("HEK293") +
  publication_theme() +
  xlab("Primary anticodon abundance") +
  ylab("Weighted Codon Usage") +
  theme(legend.direction = "vertical", 
        legend.position = "none")


png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/HEK293s_tRNAvsWCU_updated.png", res = 300, height = 1500, width = 1500)
print(HEK_updated)
dev.off()
























HEK <- gplots$HEK293T

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/HEK293s_tRNAvsWCU.png", res = 300, height = 1500, width = 1500)
print(HEK)
dev.off()

nTEs <- 
data_split_freq$HEK293T %>% 
  mutate(ratio = mean_wcu/norm_counts * -1)

write.csv(nTEs,
          "C:/Users/Jettles/OneDrive - University of Glasgow/Documents/HEK293_tRNAvsWCU.csv")



codon_usage <- read_csv("\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/bioinformatic_analysis/CU_Hg38.csv")
colnames(codon_usage) <- c("codon", "freq")
codon_usage$codon <- gsub("T","U", codon_usage$codon)

plot_tRNAs2 <- function(x){
  
  ggtit <- x$sample %>% unique()
  
  gplot <- 
    x %>% 
    ggplot(aes(x = norm_counts, y = freq, colour = leucyl)) +
    geom_point(size = 2) +
    geom_smooth(method='lm', formula= y~x)+
    geom_abline(slope = 0, intercept = 0, linewidth = 1, linetype = "dashed") +
    scale_color_manual(values = c("#999999", "#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F")) +
    ggtitle(paste(ggtit)) +
    publication_theme() +
    xlab("Primary anticodon abundance") +
    ylab("Genomic codon frequency") +
    theme(legend.direction = "vertical", 
          legend.position = "none")
  
  return(gplot)
  
}


codons_all <- inner_join(codon_anticodon, codon_usage, by = c("codon")) 
data_split_freq <- 
  lapply(data_split, inner_join, codons_all, by = c("anticodon", "AA"))
data_split_freq <- lapply(lapply(data_split_freq, attach_leucyl), apply_factor)
gplots <- lapply(data_split_freq, plot_tRNAs2)


setwd("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures")

for(i in 1:length(gplots)){

  cell_line <- names(gplots)[i]

  png(filename = paste0("tRNA_abundance_", cell_line, ".png"), res = 300, height = 1100, width = 1100)
  print(gplots[[i]])
  dev.off()
}


