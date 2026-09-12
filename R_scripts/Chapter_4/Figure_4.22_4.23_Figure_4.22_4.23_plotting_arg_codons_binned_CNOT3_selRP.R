library(tidyverse)
library(data.table)


# publication theme
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
            legend.key.size= unit(0.4, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

setwd("C:/Users/Jettles/OneDrive - University of Glasgow/Documents")

master <- as_tibble(fread("master.csv", header = T, drop = "V1"))

CNOT3_transcripts <- fread("CNOT3_monosomes_transcripts.csv", 
                           drop = "V1",
                           header = T) %>% pull(var = x)

deseq <- 
  read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/Analysis/DESeq2_output/RPFs_CNOT3_DEseq2_apeglm_LFC_shrinkage_M_28nts.csv") %>% 
  filter(padj < 0.05)

enriched <- deseq %>% filter(log2FoldChange > 0)
depleted <- deseq %>% filter(log2FoldChange < 0)



master_filt <- 
  master %>% 
  filter(ENST %in% unique(deseq$transcript)) 


master_filt$length <- nchar(master_filt$nucleotide_sequence_CDS)/3

# filter for bare minimum length of CDS which is 100 codons
master_filt <- master_filt %>% filter(length > 100)


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

#The following function creates a new column denoting the bin
calculate_bins <- function(x, bins_n = 50) {
  bins <- bins_n
  cut_size <- 1 / bins
  breaks <- seq(0, 1, cut_size)
  start <- cut_size / 2
  stop <- 1 - start
  bin_centres <- seq(start, stop, cut_size)
  bin <- as.numeric(cut(x, breaks, include.lowest = F, labels = bin_centres))
  return(bin)
}

transcript_processing <- function(x){
  
  
  transscript <- master_filt[x,] %>% pull(ENST)
  seq <- master_filt[x,] %>% pull(nucleotide_sequence_CDS) %>% split_into_codons()
  transcript_length <- master_filt[x,] %>% pull(length)
  position <- 1: transcript_length
  
  
  final <- 
    data.frame(ENST = transscript,
               position = position,
               normalised_position = position/transcript_length,
               codon = seq) %>% 
    as_tibble() %>% 
    mutate(bin = calculate_bins(normalised_position, 50))
  
  return(final)
}


data <- do.call("rbind",
                lapply(1:nrow(master_filt), transcript_processing))

data <- 
  data %>% 
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A", "U", "G","C"))) 

data <- 
data %>% 
  mutate(deseq  = case_when(ENST %in% c(enriched$transcript) ~ "enriched",
                   ENST %in% c(depleted$transcript) ~ "depleted"))

codons_n <- 
  data %>% 
  group_by(deseq, bin, codon) %>% 
  tally()

total_codons <- 
  data %>% 
  group_by(deseq, bin) %>% 
  tally() %>% 
  dplyr::rename(tot_n = n)


final <- inner_join(codons_n, total_codons, by = c("deseq", "bin"))

percentage_data <- 
final %>%
  mutate(bin = bin *2) %>% 
  mutate(percentage = (n/tot_n) * 100)

plot_CU <- function(x){
  
  gplot <- 
    x %>% ggplot(aes(x = bin, y = percentage, colour = deseq)) +
    geom_point() +
    geom_line() +
    geom_smooth() +
    ylab("Codon Percentage") +
    xlab("CDS%") +
    scale_color_manual(values = c("#1B9E77", "#7570B3")) +
    scale_y_continuous(limits = c(0,NA)) +
    ggtitle(paste0(unique(x$codon))) +
    publication_theme() +
    theme(legend.title = element_blank(),
          legend.direction = "vertical")
          
  
  return(gplot)
}

percentage_data$deseq <- factor(percentage_data$deseq,
                                levels = c("depleted", "enriched"),
                                labels = c("Transcripts depleted in CNOT3 IP",
                                           "Transcripts enriched in CNOT3 IP"))

gplots <- lapply(split(percentage_data, percentage_data$codon), plot_CU)

gplots$cga
gplots$cgg
gplots$agg

leg <- cowplot::get_legend(gplots$cga)

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/ARG_binned_legend.png", res = 300, height = 200, width = 1500)
cowplot::ggdraw(leg)
dev.off()


#RColorBrewer::brewer.pal(4, "Dark2")

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/ARG_cga_deseq_binned.png", res = 300, height = 1000, width = 1250)
print(gplots$cga)
dev.off()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/ARG_cgg_deseq_binned.png", res = 300, height = 1000, width = 1250)
print(gplots$cgg)
dev.off()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/ARG_agg_deseq_binned.png", res = 300, height = 1000, width = 1250)
print(gplots$agg)
dev.off()

aa_properties <- 
read_csv("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/useful_tables/codon_amino_acid_pairs_and_properties.csv")

aa_properties <- 
aa_properties %>% 
  mutate(codon = gsub("u", "t", tolower(codon)))

percentage_data_full <- 
percentage_data %>% 
  inner_join(aa_properties, by = "codon") %>% 
  select(-Symbol)

percentage_data_full %>% 
  group_by(bin, Properties) %>% 
  summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = bin, y = percentage, colour = Properties)) +
  geom_line()

percentage_data_full %>% 
  group_by(bin, AA) %>% 
  summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = bin, y = percentage)) +
  geom_line() +
  facet_wrap(~AA)

plot_AA_CU <- function(x){
  
  x %>% 
    mutate(wobble = 
             factor(str_sub(codon, 3,3), 
                    levels = c("a", "t", "g", "c"), 
                    labels = c("A/U", "A/U", "G/C","G/C"))) %>% 
    ggplot(aes(x = bin, y = percentage, group = codon, colour = wobble)) +
    geom_line() +
    scale_y_continuous(limits = c(0,NA)) +
    ggtitle(paste0(unique(x$AA))) +
    publication_theme()
  
}

AA_gplots <- lapply(split(percentage_data_full,percentage_data_full$AA), plot_AA_CU)

AA_gplots[[1]]

plot_AA_usage <- function(x){
  
  x %>% 
    group_by(bin, AA) %>% 
    summarise(percent = sum(percentage)) %>% 
    ggplot(aes(x = bin, y = percent, colour = AA)) +
    geom_line() +
    scale_y_continuous(limits = c(0,NA)) +
    ggtitle(paste0(unique(x$Properties))) +
    publication_theme()
  
}

lapply(split(percentage_data_full, percentage_data_full$Properties), plot_AA_usage)


CU_binned <- 
   %>% 
  ggplot(aes(x = bin, y = percentage, group = codon)) +
  facet_wrap(~codon) +
  geom_point() +
  geom_line() +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, NA)) +
  xlab("CDS%") +
  ylab("Codon%") +
  ggtitle("CNOT3 depleted") +
  publication_theme() 

png("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/thesis/CU_binned_deseq_depleted.png", res = 300, height = 1000, width = 1500)
print(CU_binned)
dev.off()


AG3 <- 
  final %>% 
  mutate(wobble = case_when(wobble == "A" | wobble == "G" ~ "AG3",
                            wobble == "U"| wobble == "C" ~ "UC3")) %>% 
  ungroup() %>% 
  group_by(bin, wobble) %>% 
  summarise(n = sum(n)) %>% 
  inner_join(final %>% select(bin, tot_n) %>% unique(), by = "bin") %>% 
  mutate(percentage = (n/tot_n) * 100)


ggplot(aes(x = bin, y = percentage, colour = wobble)) +
  geom_line() +
  scale_y_continuous(limits = c(0, NA))

GC3 <- 
  final %>% 
  mutate(wobble = case_when(wobble == "G" | wobble == "C" ~ "GC3",
                            wobble == "A"| wobble == "U" ~ "AU3")) %>% 
  ungroup() %>% 
  group_by(bin, wobble) %>% 
  summarise(n = sum(n)) %>% 
  inner_join(final %>% select(bin, tot_n) %>% unique(), by = "bin") %>% 
  mutate(percentage = (n/tot_n) * 100)

GU3 <- 
  final %>% 
  mutate(wobble = case_when(wobble == "G" | wobble == "U" ~ "GU3",
                            wobble == "A"| wobble == "C" ~ "AC3")) %>% 
  ungroup() %>% 
  group_by(bin, wobble) %>% 
  summarise(n = sum(n)) %>% 
  inner_join(final %>% select(bin, tot_n) %>% unique(), by = "bin") %>% 
  mutate(percentage = (n/tot_n) * 100)

AG3GC3GU3 <- 
  rbind(rbind(AG3, GC3), GU3) 

AG3GC3GU3$wobble <- factor(AG3GC3GU3$wobble, levels = c("GC3", "GU3", "AG3",
                                                        "AU3", "AC3", "UC3"))
AG3GC3GU3 %>% 
  ggplot(aes(x = bin, y = percentage)) +
  #geom_point() +
  geom_line(size = 1) +
  facet_wrap(~wobble) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme() +
  ylab("Percentage")


ggplot(aes(x = bin, y = percentage, colour = wobble)) +
  geom_line() +
  scale_y_continuous(limits = c(0, NA))


