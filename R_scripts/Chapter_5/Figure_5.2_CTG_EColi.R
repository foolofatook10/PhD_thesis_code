library(data.table)
library(tidyverse)
library(ggdendro)
library(RColorBrewer)

setwd("\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/bioinformatic_analysis/plots")

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


# read in zebrafish table
bacteria_transcriptome <- fread("\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/e_coli.csv",
                                 header = T,
                                 drop = "V1")

bacteria_transcriptome <- as_tibble(bacteria_transcriptome)

#Import codons table
path = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/useful_tables"
codons <- read_csv(paste0(path, "/codon_box_types.csv"))

#Make the capitalised codons lower case
codons$codon <- factor(tolower(codons$codon))

#replace letters "u" with "t"
codons$codon <- gsub("u","t", codons$codon)

#Extract a vector of the codons 
code_ons <- as.character(codons$codon)

#Append the stop codons to vector code_ons
code_ons_plus_stop <- append(code_ons, c("taa", "tag", "tga"))

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

# make table with feature information of the longest transcript for each gene
combined <- bacteria_transcriptome

# obtain numbers of transcripts with and without signal sequences
transcripts_with_signal_sequence <-  table(combined$signal_sequence)[["TRUE"]]
transcripts_without_signal_sequence <-  table(combined$signal_sequence)[["FALSE"]]

print(paste("transcripts with signal sequence:", transcripts_with_signal_sequence))
print(paste("transcripts without signal sequence:", transcripts_without_signal_sequence))


# make a list with the nucleotide sequences chopped into codons
my_lyst <- sapply(combined$nucleotide_sequence_CDS, split_into_codons)

# Name each element of the list by its ENST tag
names(my_lyst) <- pull(combined, var = gene)

# obtain logical vector of CDSs that are greater than 100 codons
over_100_codons <- lapply(my_lyst, length) > 100

# filter my_lyst for CDSs that are greater than 100 codons
my_lyst <-  my_lyst[c(over_100_codons)]

# find out how many transcripts contain/ do not contain signal sequences following filtering
transcripts_post_filtering <- 
  combined %>% 
  filter(gene %in% names(my_lyst)) %>% 
  pull(var = signal_sequence) %>% 
  table()

transcripts_with_signal_sequence_post_filtering <- 
  transcripts_post_filtering[["TRUE"]]

transcripts_without_signal_sequence_post_filtering <- 
  transcripts_post_filtering[["FALSE"]]

print(paste("transcripts with signal sequence and over 100 codons:", transcripts_with_signal_sequence_post_filtering))
print(paste("transcripts without signal sequence and over 100 codons:", transcripts_without_signal_sequence_post_filtering))

# only consider the first 50 codons per sequence
codonsN <-lapply(my_lyst, head, n=50)

# condense list into a tbale format to make it easier to work with
codonsN_table <- data.frame(do.call("rbind", codonsN))

# label columns with numbers from 1 to 50
colnames(codonsN_table) <- c(paste0("codon",1:50))

# convert rownames to column with name transcript
codonsN_table <- 
  codonsN_table %>% rownames_to_column("gene")

# convert table to long format
codonsN_table_long <- 
  codonsN_table %>% 
  gather(key = "position", value =  "codon", codon1:codon50)

# arrange by transcript as a sanity check
sanity_check <- codonsN_table_long %>% arrange(gene)

# add amino acid information and signal sequence information to the table
codonsN_table_long <- 
  codonsN_table_long %>%
  inner_join(codons %>% select(codon, AA), by = "codon") %>%
  inner_join(combined %>% select(gene, signal_sequence), by = "gene")

# convert positions to numerics
codonsN_table_long$position <- as.numeric(str_remove(codonsN_table_long$position, pattern = "codon"))

# calculate frequency of codon usage per position
percentages <- cbind(codonsN_table_long %>%
                       group_by(signal_sequence, position) %>%
                       dplyr::count(vars = codon),
                     codonsN_table_long %>%
                       group_by(signal_sequence, position) %>%
                       dplyr::count(vars = codon) %>%
                       group_by(signal_sequence, position) %>%
                       summarise(percent = n/sum(n) * 100) %>%
                       ungroup() %>%
                       select(percent))

# the frequencies of each codon should add up to 100% at each position
sanity_check2 <- 
  percentages %>% 
  group_by(signal_sequence, position) %>% 
  summarise(sum_perc = sum(percent)) %>% 
  pull(var = sum_perc) 

# add column for wobble base position
percentages <- percentages %>%
  dplyr::rename(codon = vars) %>%
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A/U", "A/U", "G/C","G/C")))

# add tags
temp <-
  percentages %>% 
  mutate(tag = case_when(signal_sequence == F ~ "No signal sequence",
                         .default = "Signal sequence"))

# overall GC3 proportion in first 50 codons
GC3 <-
  ggplot(temp, aes(x = position, y = percent, fill = wobble)) +
  geom_bar(stat = "identity") +
  xlab("Codon") +
  ylab("Percent")+
  ggtitle("E. Coli")+
  facet_wrap(~tag, nrow  = 1) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
  publication_theme()+
  theme(legend.position = "right",
        legend.direction = "vertical")
saveRDS(GC3, file = "\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/bioinformatic_analysis/scripts/ggplot_objects/zebrafishGC3.rdata")

tiff(file.path("GC3_zebrafish.tiff"), height = 4, width = 5, units = "in", res = 500)
print(GC3)
dev.off()

# vector of leucyl codons
leucyl_codons <- c(codons %>% filter(AA == "Leu") %>% pull(var = codon))

# tag the ctgs or the leucyl codons in seperate columns
lines_processed <- 
  temp %>% 
  mutate(ctg = case_when(codon == "ctg" ~ "ctg",
                         .default = "other"),
         ctc = case_when(codon == "ctc" ~ "ctc",
                         .default = "other"),
         cta = case_when(codon == "cta" ~ "cta",
                         .default = "other"),
         ctt = case_when(codon == "ctt" ~ "ctt",
                         .default = "other"),
         ttg = case_when(codon == "ttg" ~ "ttg",
                         .default = "other"),
         tta = case_when(codon == "tta" ~ "tta",
                         .default = "other"),
         gcg = case_when(codon == "gcg" ~ "gcg",
                         .default = "other"),
         gcc = case_when(codon == "gcc" ~ "gcc",
                         .default = "other"),
         gca = case_when(codon == "gca" ~ "gca",
                         .default = "other"),
         gct = case_when(codon == "gct" ~ "gct",
                         .default = "other"),
         other = "all_codons",
         leucyl = case_when(codon == "ctt" ~ "ctt",
                            codon == "cta" ~"cta",
                            codon == "ctg" ~ "ctg",
                            codon == "tta" ~ "tta",
                            codon == "ttg" ~"ttg",
                            codon == "ctc" ~"ctc",
                            .default = "other")) %>% 
  ungroup()

# make leucyl codons column a factor
lines_processed$leucyl <- factor(lines_processed$leucyl, 
                                 levels = c("other","ctg", "ctc", "ctt", "ttg", "cta", "tta"))


# plot the leucyl codon usage for transcripts with and without a signal sequence
lines_gplot <- 
  lines_processed %>% 
  ggplot(aes(x = position, y = percent, group = codon, color = leucyl)) +
  geom_line(size = 1, alpha = 0.5) +
  scale_color_manual(values = c("#999999", "#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F")) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  ggtitle("E. Coli transcripts") +
  facet_wrap(~tag) +
  publication_theme() +
  theme(legend.position = "right", 
        legend.direction = "vertical", 
        legend.title = element_blank()) 


lines_gplot <- 
  lines_processed %>% 
  ggplot() +
  geom_line(data = lines_processed %>% 
              filter(leucyl == "other"), 
            aes(x = position, 
                y = percent,  
                group = codon, 
                color = ctt), 
            size = 1, 
            #alpha = 0.5, 
            color = "#999999") +
  geom_line(data = lines_processed %>% 
              filter(ctg == "ctg"), 
            aes(x = position, 
                y = percent,  
                group = codon, 
                color = ctt), 
            size = 1, 
            #alpha = 0.5, 
            color = "#377EB8")+
  geom_line(data = lines_processed %>% 
              filter(ctc == "ctc"), 
            aes(x = position, 
                y = percent,  
                group = codon, 
                color = ctt), 
            size = 1, 
            #alpha = 0.5, 
            color = "#984EA3") +
  geom_line(data = lines_processed %>% 
              filter(ctt == "ctt"), 
            aes(x = position, 
                y = percent,  
                group = codon, 
                color = ctt), 
            size = 1, 
            #alpha = 0.5, 
            color = "#4DAF4A")+
  geom_line(data = lines_processed %>% 
              filter(cta == "cta"), 
            aes(x = position, 
                y = percent,  
                group = codon, 
                color = ctt), 
            size = 1, 
            #alpha = 0.5, 
            color = "#BF5B17")+
  geom_line(data = lines_processed %>% 
              filter(ttg == "ttg"), 
            aes(x = position, 
                y = percent,  
                group = codon, 
                color = ctt), 
            size = 1, 
            #alpha = 0.5, 
            color = "#FF7F00") +
  geom_line(data = lines_processed %>% 
              filter(tta == "tta"), 
            aes(x = position, 
                y = percent,  
                group = codon, 
                color = ctt), 
            size = 1, 
            #alpha = 0.5, 
            color = "#F0027F") +
  #scale_color_manual(values = c("#999999", "#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F")) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  ggtitle("E. Coli transcripts") +
  facet_wrap(~tag) +
  publication_theme() +
  theme(legend.position = "right", 
        legend.direction = "vertical", 
        legend.title = element_blank()) 


lines_gplot <- 
  lines_processed %>% 
  ggplot() +
  geom_line(data = lines_processed %>% 
              filter(other == "all_codons"), 
            aes(x = position, 
                y = percent,  
                group = codon, 
                color = ctt), 
            size = 1, 
            #alpha = 0.5, 
            color = "#999999") +
  geom_line(data = lines_processed %>% 
              filter(gca == "gca"), 
            aes(x = position, 
                y = percent,  
                group = codon, 
                color = gca), 
            size = 1, 
            #alpha = 0.5, 
            color = "#377EB8")+
  geom_line(data = lines_processed %>% 
              filter(gcc == "gcc"), 
            aes(x = position, 
                y = percent,  
                group = codon, 
                color = gcc), 
            size = 1, 
            #alpha = 0.5, 
            color = "#984EA3") +
  geom_line(data = lines_processed %>% 
              filter(gct == "gct"), 
            aes(x = position, 
                y = percent,  
                group = codon, 
                color = gct), 
            size = 1, 
            #alpha = 0.5, 
            color = "#4DAF4A")+
  geom_line(data = lines_processed %>% 
              filter(gcg == "gcg"), 
            aes(x = position, 
                y = percent,  
                group = codon, 
                color = gcg), 
            size = 1, 
            #alpha = 0.5, 
            color = "#BF5B17") +
  facet_wrap(~tag) +
  #scale_color_manual(values = c("#999999", "#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F")) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  ggtitle("E. Coli transcripts") +
  
  publication_theme() +
  theme(legend.position = "right", 
        legend.direction = "vertical", 
        legend.title = element_blank()) 

lines_processed_split <- split(lines_processed, f = lines_processed$position)

lines_processed_split_arranged <- lapply(lines_processed_split, function(x){x %>% filter(signal_sequence == T) %>%  arrange(-percent)})

lines_processed_split_arranged[[24]]

tiff(file.path("CTGenrichment_EColi.tiff"), height = 3, width = 6, units = "in", res = 500)
print(lines_gplot)
dev.off()
