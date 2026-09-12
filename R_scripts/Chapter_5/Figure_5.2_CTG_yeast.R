library(data.table)
library(tidyverse)
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


# read in yeast table
yeast_transcriptome <- fread("\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/saccharomyces_cerevisiae.csv",
                                 header = T,
                                 drop = "V1")

yeast_transcriptome <- as_tibble(yeast_transcriptome)

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
combined <- yeast_transcriptome

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
  mutate(tag = case_when(signal_sequence == F ~ paste0("No signal sequence","\nn = ", 
                                                       transcripts_without_signal_sequence_post_filtering),
                         .default = paste0("Signal sequence","\nn = ", transcripts_with_signal_sequence_post_filtering)))
# overall GC3 proportion in first 50 codons
GC3 <-
  ggplot(temp, aes(x = position, y = percent, fill = wobble)) +
  geom_bar(stat = "identity") +
  xlab("Codon") +
  ylab("Percent")+
  ggtitle("Human")+
  facet_wrap(~tag, nrow  = 1) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  publication_theme()

# tiff(file.path("GC3_content.tiff"), height = 5, width = 5, units = "in", res = 200)
# print(GC3)
# dev.off()

# vector of leucyl codons
leucyl_codons <- c(codons %>% filter(AA == "Leu") %>% pull(var = codon))

# tag the ctgs or the leucyl codons in seperate columns
lines_processed <- 
  temp %>% 
  mutate(ctg = case_when(codon == "ctg" ~ "ctg",
                         .default = "other"),
         leucyl = case_when(codon == "ctt" ~ "ctt",
                            codon == "cta" ~"cta",
                            codon == "ctg" ~ "ctg",
                            codon == "tta" ~ "tta",
                            codon == "ttg" ~"ttg",
                            codon == "ctc" ~"ctc",
                            .default = "other"),
         alanine_codon = case_when(
           codon == "gcc" ~ "Ala-gcc",
           codon == "gca" ~ "Ala-gca",
           codon == "gct" ~ "Ala-gct",
           codon == "gcg" ~ "Ala-gcg",
           codon == "tta" ~ "Leu-tta",
           codon == "ttg" ~ "Leu-ttg",
           codon == "gtt" ~ "Val-gtt",
           codon == "gaa" ~ "Val-gaa",
           .default = "other")) %>% 
  ungroup()


# vector of leucyl codons
leucyl_codons <- c(codons %>% filter(AA == "Leu") %>% pull(var = codon))

# tag the ctgs or the leucyl codons in seperate columns
lines_processed <- 
  temp %>% 
  mutate(ctg = case_when(codon == "ctg" ~ "ctg",
                         .default = "other"),
         leucyl = case_when(codon == "ctt" ~ "ctt",
                            codon == "cta" ~"cta",
                            codon == "ctg" ~ "ctg",
                            codon == "tta" ~ "tta",
                            codon == "ttg" ~"ttg",
                            codon == "ctc" ~"ctc",
                            .default = "other")) %>% 
  ungroup()

# make leucyl codons column a factor
lines_processed$leucyl <- factor(lines_processed$leucyl, levels = c("ctg", "ctc", "ctt", "ttg", "cta", "tta","other"))

#
background <- lines_processed %>% filter(leucyl == "other" & position != 1)
leucyl2plot <- lines_processed %>% filter(leucyl %in% c("ctg", "ctc", "ctt", "ttg", "cta", "tta") & position != 1)

leucyl_colours = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F")

final_gplot <- 
  ggplot() +
  geom_line(data = background,
            aes(x = position, y = percent, group = codon),
            size = 1,
            colour = "#C8C8C8") +
  geom_line(data = leucyl2plot,
            aes(x = position, y = percent, group = codon, colour = leucyl),
            size = 1) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  scale_color_manual(values = leucyl_colours) +
  ggtitle("Saccharomyces cerevisiae") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(0,50)) +
  scale_x_continuous(limits = c(2,50), breaks = c(2,10,20,30,40,50)) +
  publication_theme() +
  theme(legend.position = "none", 
        legend.direction = "vertical", 
        legend.title = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/s_cerevisiae.png",
    res = 300,
    height = 1000,
    width = 1500)
print(final_gplot)
dev.off()






# make leucyl codons column a factor
lines_processed$alanine_codon <- factor(lines_processed$alanine_codon, levels = c("Ala-gcc", "Ala-gca", "Ala-gct", "Ala-gcg", "Leu-tta", "Leu-ttg", "Val-gtt", "Val-gaa", "other"))


# plot the leucyl codon usage for transcripts with and without a signal sequence
lines_gplot <- 
  lines_processed %>% 
  ggplot(aes(x = position, y = percent, group = codon, color = alanine_codon)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F", "black", "red", "#999999")) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  ggtitle("yeast transcripts") +
  facet_wrap(~tag) +
  publication_theme() +
  theme(legend.position = "right", 
        legend.direction = "vertical", 
        legend.title = element_blank()) 


tiff(file.path("CTGenrichment_yeast.tiff"), height = 3, width = 6, units = "in", res = 500)
print(lines_gplot)
dev.off()
