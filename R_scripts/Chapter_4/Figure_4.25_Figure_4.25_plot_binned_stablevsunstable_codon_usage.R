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
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

setwd("C:/Users/Jettles/OneDrive - University of Glasgow/Documents")

master <- as_tibble(fread("\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/master.csv", header = T, drop = "V1"))

CNOT3_transcripts <- fread("CNOT3_monosomes_transcripts.csv", 
                           drop = "V1",
                           header = T) %>% pull(var = x)

stabilities <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/stability_extremes.csv") %>% select(-`...1`) %>% 
  dplyr::rename(transcript = "ENST")

master_filt <- 
  master %>% 
  filter(ENST %in% unique(stabilities$transcript)) 

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
                  labels = c("AU3", "AU3", "GC3","GC3"))) 

data <- 
data %>% 
  inner_join(stabilities, by = c("ENST" = "transcript")) %>% 
  select(-Gene)

codons_n <- 
  data %>% 
  group_by(stability, bin, wobble) %>% 
  tally()

total_codons <- 
  data %>% 
  group_by(stability, bin) %>% 
  tally() %>% 
  dplyr::rename(tot_n = n)


final <- inner_join(codons_n, total_codons, by = c("stability", "bin"))

final$stability <- factor(final$stability, 
                          levels = c("unstable", "stable"),
                          labels = c("Unstable transcripts\nn = 1059",
                                     "Stable transcripts\nn = 1207"))

CU_binned <- 
  final %>%
  mutate(bin = bin *2) %>% 
  mutate(percentage = (n/tot_n) * 100) %>% 
  filter(wobble == "GC3") %>% 
  ggplot(aes(x = bin, y = percentage, colour = stability)) +
  geom_point() +
  geom_line() +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, NA)) +
  xlab("CDS%") +
  ylab("GC3%") +
  #ggtitle("Unstable transcripts (n = 1121)") +
  publication_theme() +
  theme(legend.title = element_blank(),
        legend.direction = "vertical",
        legend.position = "right")

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CU_binned_stabilities_unstable.png", res = 300, height = 1500, width = 2000)
print(CU_binned)
dev.off()

