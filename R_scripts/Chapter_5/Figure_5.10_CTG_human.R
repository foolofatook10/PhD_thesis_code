library(data.table)
library(tidyverse)
#library(ggdendro)
library(RColorBrewer)
library(data.table)

setwd("\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/bioinformatic_analysis/plots")


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

# read in publication theme
publication_theme2 <- function(base_size=14, base_family="helvetica") {
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
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


# read in features table
features <- as_tibble(fread("\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/features.csv", drop = "V1")) %>%
  select(ENSG, ENST, signal_sequence)

# tpms <- 
#   read_csv("\\\\data.beatson.gla.ac.uk/data/R11/external_sequencing_data/Gillen_2021_RiboSeq/Analysis/DESeq2_output/tpms.csv")
# 
# expression <- 
# tpms %>% 
#   gather(key = sample, value = tpms, Ctrl_1_Totals:CNOT1_3_Totals) %>% 
#   mutate(replicate = str_remove_all(str_extract(sample, pattern = "_\\d_"), "_"),
#          sample = str_remove(sample, "_\\d_.*")) %>% 
#   inner_join(features %>% select(-ENSG), by = c("transcript" = "ENST"))
# 
# expression %>% 
#   group_by(signal_sequence, transcript,sample) %>% 
#   summarise(mean_tpm = mean(tpms),
#             sd_tpm = sd(tpms)) %>% 
#   ggplot(aes(x = signal_sequence, y = log10(mean_tpm))) +
#   geom_violin()+
#   geom_boxplot()+
#   facet_wrap(~sample) 


# read in master table
master <- as_tibble(fread(file = "\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/master.csv", header = T, drop = "V1"))

# select the longest transcript per gene
longest <- master %>% 
  mutate(polypeptide_length = nchar(master$amino_acid_sequence)) %>%
  group_by(ENSG) %>%
  mutate(the_rank  = rank(-polypeptide_length, ties.method = "first")) %>%
  filter(the_rank == 1) %>% 
  select(-the_rank)

# filter the features tablefor unique transcripts
ENSTs <- unique(features$ENST)

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
combined <- inner_join(features, longest, by = c("ENSG","ENST"))

# obtain numbers of transcripts with and without signal sequences
transcripts_with_signal_sequence <-  table(combined$signal_sequence)[["TRUE"]]
transcripts_without_signal_sequence <-  table(combined$signal_sequence)[["FALSE"]]

print(paste("transcripts with signal sequence:", transcripts_with_signal_sequence))
print(paste("transcripts without signal sequence:", transcripts_without_signal_sequence))


# make a list with the nucleotide sequences chopped into codons
my_lyst <- sapply(combined$nucleotide_sequence_CDS, split_into_codons)

# Name each element of the list by its ENST tag
names(my_lyst) <- pull(combined, var = ENST)

# obtain logical vector of CDSs that are greater than 100 codons
over_100_codons <- lapply(my_lyst, length) > 100

# filter my_lyst for CDSs that are greater than 100 codons
my_lyst <-  my_lyst[c(over_100_codons)]

# find out how many transcripts contain/ do not contain signal sequences following filtering
transcripts_post_filtering <- 
  combined %>% 
  filter(ENST %in% names(my_lyst)) %>% 
  pull(var = signal_sequence) %>% 
  table()

transcripts_with_signal_sequence_post_filtering <- 
  transcripts_post_filtering[["TRUE"]]

transcripts_without_signal_sequence_post_filtering <- 
  transcripts_post_filtering[["FALSE"]]

print(paste("transcripts with signal sequence and over 100 codons:", transcripts_with_signal_sequence_post_filtering))
print(paste("transcripts without signal sequence and over 100 codons:", transcripts_without_signal_sequence_post_filtering))


# only consider the first 50 codons per sequence
codonsN <-lapply(my_lyst, head, n=100)

# rest of the transcript
codons_tail <- lapply(my_lyst, tail, n = 100)

my_lyst[[1]]
codonsN[[1]]
codons_tail[[1]]

# condense list into a tbale format to make it easier to work with
codonsN_table <- data.frame(do.call("rbind", codonsN))
codons_tail_table <- data.frame(do.call("rbind", codons_tail))

max_n <- max(as.numeric(str_remove(colnames(codonsN_table), "X")))
max_n_tail <- -1 * max(as.numeric(str_remove(colnames(codons_tail_table), "X")))

# label columns with numbers from 1 to 50
colnames(codonsN_table) <- c(paste0("codon",1:max_n))
colnames(codons_tail_table) <- c(paste0("codon",max_n_tail:-1))

# convert rownames to column with name transcript
codonsN_table <- 
  codonsN_table %>% rownames_to_column("ENST")
codons_tail_table <- 
  codons_tail_table %>% rownames_to_column("ENST")

# convert table to long format
codonsN_table_long <- 
  codonsN_table %>% 
  gather(key = "position", value =  "codon", codon1:paste0("codon", max_n))
codons_tail_table_long <- 
  codons_tail_table %>% 
  gather(key = "position", value =  "codon", `codon-100`:paste0("codon-1"))

codons_tail_table_long %>% filter(position == "codon-1") %>% pull(var = codon) %>% table()

# arrange by transcript as a sanity check
sanity_check <- codonsN_table_long %>% arrange(ENST)

codons <- 
rbind(codons,
data.frame(codon = c("taa", "tga", "tag"), AA = rep("Stp",3), box = "3-box"))

# add amino acid information and signal sequence information to the table
codonsN_table_long <- 
  codonsN_table_long %>%
  inner_join(codons %>% select(codon, AA), by = "codon") %>%
  inner_join(features %>% select(ENST, signal_sequence), by = "ENST")

codons_tail_table_long <- 
  codons_tail_table_long %>%
  inner_join(codons %>% select(codon, AA), by = "codon") %>%
  inner_join(features %>% select(ENST, signal_sequence), by = "ENST")

# convert positions to numerics
codonsN_table_long$position <- as.numeric(str_remove(codonsN_table_long$position, pattern = "codon"))
codons_tail_table_long$position <- as.numeric(str_remove(codons_tail_table_long$position, pattern = "codon"))

raw_counts <- 
codonsN_table_long %>%
  group_by(signal_sequence, position) %>%
  dplyr::count(vars = codon) %>% 
  ungroup()

raw_counts_tail <- 
codons_tail_table_long %>%
  group_by(signal_sequence, position) %>%
  dplyr::count(vars = codon) %>% 
  ungroup()
  

# calculate frequency of codon usage per position
percentages <- cbind(raw_counts,
                     raw_counts %>%
                       group_by(signal_sequence, position) %>%
                       dplyr::reframe(percent = n/sum(n) * 100) %>%
                       ungroup() %>%
                       select(percent))

percentages_tail <- cbind(raw_counts_tail,
                          raw_counts_tail %>%
                       group_by(signal_sequence, position) %>%
                        dplyr::reframe(percent = n/sum(n) * 100) %>%
                       ungroup() %>%
                       select(percent))

percentages_tail %>% filter(position == -1)

# the frequencies of each codon should add up to 100% at each position
sanity_check2 <- 
  percentages %>% 
  group_by(signal_sequence, position) %>% 
  summarise(sum_perc = sum(percent)) %>% 
  pull(var = sum_perc) 

percentages_tail %>% 
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

percentages_tail <- percentages_tail %>%
  dplyr::rename(codon = vars) %>%
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A/U", "A/U", "G/C","G/C")))

# add tags
temp <-
  percentages %>% 
  mutate(tag = case_when(signal_sequence == F ~ paste0("No signal sequence", 
                                                      "\nn = ", 
                                                      transcripts_without_signal_sequence_post_filtering),
                         .default = paste0("Signal sequence", "\nn = ", transcripts_with_signal_sequence_post_filtering)))

temp_tail <-
  percentages_tail %>% 
  mutate(tag = case_when(signal_sequence == F ~ paste0("No signal sequence", 
                                                       "\nn = ", 
                                                       transcripts_without_signal_sequence_post_filtering),
                         .default = paste0("Signal sequence", "\nn = ", transcripts_with_signal_sequence_post_filtering)))

# # overall GC3 proportion in first 50 codons
# # GC3 <-
# #   ggplot(temp, aes(x = position, y = percent, fill = wobble)) +
# #   geom_bar(stat = "identity") +
# #   xlab("Codon") +
# #   ylab("Percent")+
# #   ggtitle("Human")+
# #   facet_wrap(~tag, nrow  = 1) +
# #   scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
# #   #coord_cartesian(xlim = c(100,200)) +
# #   publication_theme() +
# #   theme(legend.position = "right",
# #         legend.direction = "vertical")
# 
# GC3 <-
#   temp %>% 
#   ungroup() %>% 
#   group_by(signal_sequence, position, wobble) %>% 
#   summarise(wobble_percentage = sum(percent)) %>% 
#   filter(wobble == "G/C") %>% 
#   ggplot(aes(x = position, y = wobble_percentage, colour = signal_sequence)) +
#   geom_line(size = 1.5) +
#   geom_point(size = 3) +
#   xlab("Codon") +
#   ylab("Percent")+
#   ggtitle("Human")+
#   #facet_wrap(~tag, nrow  = 1) +
#   scale_colour_manual(values = c("#66C2A5", "#FC8D62")) +
#   ylim(0,100) +
#   coord_cartesian(xlim = c(0,150)) +
#   publication_theme() +
#   theme(legend.position = "right",
#         legend.direction = "vertical")
# 
# get_only_legend <- function(plot) { 
#   
#   # get tabular interpretation of plot 
#   plot_table <- ggplot_gtable(ggplot_build(plot))  
#   
#   #  Mark only legend in plot 
#   legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")  
#   
#   # extract legend 
#   legend <- plot_table$grobs[[legend_plot]] 
#   
#   # return legend 
#   return(legend)  
# }
# 
# leg <- get_only_legend(GC3)

# GC3 <-
#   ggplot(temp, aes(x = position, y = percent, fill = wobble)) +
#   geom_bar(stat = "identity") +
#   xlab("Codon") +
#   ylab("Percent")+
#   ggtitle("Human")+
#   facet_wrap(~tag, nrow  = 1) +
#   scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
#   publication_theme() +
#   theme(legend.position = "right",
#         legend.direction = "vertical")

# saveRDS(GC3, file = "\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/bioinformatic_analysis/scripts/ggplot_objects/humanGC3.rdata")
# saveRDS(leg, file = "\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/bioinformatic_analysis/scripts/ggplot_objects/legend.rdata")

 # tiff(file.path("GC3_human.tiff"), height = 300, width = 400)
 # print(GC3)
 # dev.off()

# vector of leucyl codons
leucyl_codons <- c(codons %>% filter(AA == "Leu") %>% pull(var = codon))

# vector of arginine codons
arginine_codons <- c(codons %>% filter(AA == "Arg") %>% pull(var = codon))

ala_codons <- codons %>% filter(AA == "Ala")
val_codons <- codons %>% filter(AA == "Val")
ile_codons <- codons %>% filter(AA == "Ile")
pro_codons <- codons %>% filter(AA == "Pro")
phe_codons <- codons %>% filter(AA == "Phe")


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
         arginine = case_when(codon == "cgt" ~ "cgt",
                              codon == "cga" ~ "cga",
                              codon == "cgg" ~ "cgg",
                              codon == "aga" ~ "aga",
                              codon == "agg" ~ "agg",
                              codon == "cgc" ~ "cgc",
                              .default = "other"),
         alanine = case_when(codon == "gct" ~ "gct",
                             codon == "gca" ~ "gca",
                             codon == "gcg" ~ "gcg",
                             codon == "gcc" ~ "gcc",
                             .default = "other"),
         valine = case_when(codon == "gtt" ~ "gtt",
                            codon == "gta" ~ "gta",
                            codon == "gtg" ~ "gtg",
                            codon == "gtc" ~ "gtc",
                            .default = "other"),
         isoleucine = case_when(codon == "att" ~ "att",
                            codon == "atc" ~ "atc",
                            codon == "ata" ~ "ata",
                            .default = "other"),
         proline = case_when(codon == "cct" ~ "cct",
                             codon == "cca" ~ "cca",
                             codon == "ccg" ~ "ccg",
                             codon == "ccc" ~ "ccc",
                             .default = "other"),
         phenylalanine = case_when(codon == "ttc" ~ "ttc",
                                   codon == "ttt" ~ "ttt",
                                   .default = "other"),
         CGNCG = case_when(
           codon == "ctg" ~ "Leu-CTG",
           codon == "ctc" ~ "Leu-CTC",
           codon == "gcc" ~ "Ala-GCC",
           codon == "gcg" ~ "Ala-GCG",
           codon == "gct" ~ "Ala-GCT",
           codon == "gca" ~ "Ala-GCA",
           codon == "gtc" ~ "Val-GTC",
           codon == "gtg" ~ "Val-GTG",
           codon == "tgg" ~ "Trp-TGG",
           .default = "other"),
         hydrophobic_AU = case_when(
           codon == "gtt" ~ "Val-GTT",
           codon == "gta" ~ "Val-GTA",
           codon == "att" ~ "Ile-ATT",
           codon == "atc" ~ "Ile-ATC",
           codon == "ata" ~ "Ile-ATA",
           codon == "ctt" ~ "Leu-CTT",
           codon == "cta" ~ "Leu-CTA",
           codon == "ttg" ~ "Leu-TTG",
           codon == "tta" ~ "Leu-TTA",
           codon == "ttt" ~ "Phe-TTT",
           codon == "ttc" ~ "Phe-TTC",
           codon == "tat" ~ "Tyr-TAT",
           codon == "tac" ~ "Tyr-TAC",
           codon == "atg" ~ "Met-ATG",
           .default = "other")) %>% 
  ungroup()

lines_processed_tail <- 
  temp_tail %>% 
  mutate(ctg = case_when(codon == "ctg" ~ "ctg",
                         .default = "other"),
         leucyl = case_when(codon == "ctt" ~ "ctt",
                            codon == "cta" ~"cta",
                            codon == "ctg" ~ "ctg",
                            codon == "tta" ~ "tta",
                            codon == "ttg" ~"ttg",
                            codon == "ctc" ~"ctc",
                            .default = "other"),
         arginine = case_when(codon == "cgt" ~ "cgt",
                              codon == "cga" ~ "cga",
                              codon == "cgg" ~ "cgg",
                              codon == "aga" ~ "aga",
                              codon == "agg" ~ "agg",
                              codon == "cgc" ~ "cgc",
                              .default = "other")) %>% 
  ungroup()

# make leucyl codons column a factor
lines_processed$leucyl <- factor(lines_processed$leucyl, levels = c("ctg", "ctc", "ctt", "ttg", "cta", "tta","other"))
lines_processed_tail$leucyl <- factor(lines_processed_tail$leucyl, levels = c("ctg", "ctc", "ctt", "ttg", "cta", "tta","other"))

lines_processed$arginine <- factor(lines_processed$arginine, levels = c("other", "cgg", "cgc", "cga", "cgt", "aga", "agg"))
lines_processed$alanine <- factor(lines_processed$alanine, levels = c("other", "gcg", "gcc", "gca", "gct"))
lines_processed$valine <- factor(lines_processed$valine, levels = c("other", "gtg", "gtc", "gta", "gtt"))
lines_processed$isoleucine <- factor(lines_processed$isoleucine, levels = c("other", "atc", "att", "ata"))
lines_processed$proline <- factor(lines_processed$proline, levels = c("other", "ccg", "ccc", "cca", "cct"))
lines_processed$phenylalanine <- factor(lines_processed$phenylalanine, levels = c("other", "ttc", "ttt"))

lines_processed$CGNCG <- factor(lines_processed$CGNCG, levels = c("other", "Leu-CTG", "Leu-CTC", "Ala-GCG", "Ala-GCC", "Val-GTC", "Val-GTG", "Ala-GCT", "Ala-GCA", "Trp-TGG"))
lines_processed$hydrophobic_AU <- factor(lines_processed$hydrophobic_AU, levels = c("other", 
                                                                                    "Ile-ATT",
                                                                                    "Ile-ATC",
                                                                                    "Ile-ATA",
                                                                                    "Leu-CTT",
                                                                                    "Leu-CTA",
                                                                                    "Leu-TTG",
                                                                                    "Leu-TTA",
                                                                                    "Val-GTT",
                                                                                    "Val-GTA",
                                                                                    "Phe-TTT",
                                                                                    "Phe-TTC",
                                                                                    "Tyr-TAT",
                                                                                    "Tyr-TAC",
                                                                                    "Met-ATG"))





CSCs <-
  read_csv("\\\\data.beatson.gla.ac.uk/data/R11/external_sequencing_data/Gillen_2021_halflife/Analysis/CSCs/CSCs.csv")%>%
  mutate(codon = tolower(codon)) %>%
  select(codon, ctrl_cor_estimate, delta_cor_estimate)

CSCs$codon <- gsub("u","t", CSCs$codon)

CSC_info <-
codonsN_table_long %>%
  inner_join(CSCs, by = "codon") %>%
  ungroup() %>%
  group_by(signal_sequence, position) %>%
  summarise(ctrl_opt = mean(ctrl_cor_estimate),
            delta_opt = mean(delta_cor_estimate))
# 
# CSC_info_C3 <-
#   codonsN_table_long %>%
#   filter(ENST %in% CNOT3_transcripts) %>%
#   inner_join(CSCs, by = "codon") %>%
#   ungroup() %>%
#   group_by(signal_sequence, position) %>%
#   summarise(ctrl_opt = mean(ctrl_cor_estimate),
#             delta_opt = mean(delta_cor_estimate))
# 
CSC_ctrl <-
CSC_info %>%
  filter(position != 1 & position <= 100) %>%
  ggplot(aes(x = position, y = ctrl_opt, colour = signal_sequence)) +
  scale_color_manual(values = c("#FF7F0EFF", "#1F77B4FF"),
                   labels = c("No signal sequene", "Signal sequence")) +
  ylab("Mean CSC") +
  xlab("Codon position") +
  geom_line(linewidth = 1) +
  geom_smooth() +
  geom_point() +
  scale_x_continuous(limits = c(2,50), breaks = c(2,10,20,30,40,50)) +
  scale_y_continuous(limits = c(-0.055, 0.00), breaks = c(0.00,-0.01,-0.02,-0.03,-0.04,-0.05)) +
  #ggtitle("") +
  publication_theme() +
  theme(legend.position = "none")
# 
png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/TranscriptomeCSC_ctrl.png",
    res = 300, height = 1000, width = 1000)
print(CSC_ctrl)
dev.off()
# 
# CSC_ctrl_C3 <-
#   CSC_info_C3 %>%
#   filter(position != 1) %>% 
#   ggplot(aes(x = position, y = ctrl_opt, colour = signal_sequence)) +
#   scale_color_manual(name = NULL, values = c("#1B9E77", "#D95F02"),
#                      labels = c("No signal sequene", "Signal sequence")) +
#   ylab("Mean CSC") +
#   xlab("Codon position") +
#   geom_line(linewidth = 1) +
#   geom_point() +
#   ggtitle("CNOT3: NTC") +
#   scale_y_continuous(limits = c(-0.075, 0.025)) +
#   publication_theme() +
#   theme(legend.position = "none")
# 
# png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3CSCs_NTC.png",
#     res = 300, height = 1000, width = 2000)
# print(CSC_ctrl_C3)
# dev.off()
# 
CSC_CNOT1KD <-
CSC_info %>%
  filter(position != 1 & position <= 50) %>%
  ggplot(aes(x = position, y = delta_opt, colour = signal_sequence)) +
  scale_color_manual(values = c("#FF7F0EFF", "#1F77B4FF"),
                     labels = c("No signal sequene", "Signal sequence")) +
  ylab(paste0("\u0394", "CSC")) +
  xlab("Codon position") +
  geom_line(linewidth = 1) +
  geom_smooth() +
  geom_point() +
  #scale_x_continuous(limits = c(2,50), breaks = c(2,10,20,30,40,50)) +
  #scale_y_continuous(limits = c(-0.075, 0.025)) +
  #ggtitle("") +
  publication_theme() +
  theme(legend.position = "none")
# 


png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/TranscriptomeCSCs_siCNOT1.png",
    res = 300, height = 1000, width = 1000)
print(CSC_CNOT1KD)
dev.off()

# CNOT3CSC_CNOT1KD <- 
# CSC_info_C3 %>%
#   filter(position != 1) %>% 
#   ggplot(aes(x = position, y = delta_opt, colour = signal_sequence)) +
#   scale_color_manual(name = NULL, values = c("#1B9E77", "#D95F02"),
#                      labels = c("No signal sequene", "Signal sequence")) +
#   ylab(paste0("\u0394", "CSC")) +
#   xlab("Codon position") +
#   geom_line(linewidth = 1) +
#   geom_point() +
#   ggtitle("CNOT3 transcripts: siCNOT1") +
#   scale_y_continuous(limits = c(-0.01, 0.15)) +
#   publication_theme() +
#   theme(legend.position = "none")
# 
# png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3CSCs_siCNOT1.png",
#     res = 300, height = 1000, width = 2000)
# print(CNOT3CSC_CNOT1KD)
# dev.off()
# 
# 
# png("\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/bioinformatic_analysis/plots/CSC_signal_peptides_Ctrl.png", res = 300, height = 1200, width = 2000)
# print(CSC_ctrl)
# dev.off()
# 
# png("\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/bioinformatic_analysis/plots/CSC_signal_peptides_CNOT1KD.png", res = 300, height = 1200, width = 2000)
# print(CSC_CNOT1KD)
# dev.off()
# 
# CSCs %>% ggplot(aes())
# 
# CSCs <-
#   read_csv("\\\\data.beatson.gla.ac.uk/data/R11/external_sequencing_data/Gillen_2021_halflife/Analysis/CSCs/CSCs.csv") %>%
#   arrange(-delta_cor_estimate) %>%
#   select(codon, delta_cor_estimate ) %>%
#   mutate(codon = tolower(codon)) %>%
#   mutate(optimality_tag = case_when(delta_cor_estimate < 0 ~ "optimal",
#                                  .default = "non-optimal"))
# 
# CSCs$codon <- gsub("u","t", CSCs$codon)
# 
# 
# nTEs <- 
# read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/HEK293_tRNAvsWCU.csv") %>% 
#   mutate(nTE = ratio * -1) %>% 
#   select(codon,nTE) %>% 
#   mutate(codon = gsub("u", "t", tolower(codon)))
# 
# nTEs_transcriptome <- 
# codonsN_table_long %>% 
#   inner_join(nTEs, by = "codon") %>%
#   ungroup() %>%
#   group_by(signal_sequence, position) %>%
#   summarise(ctrl_nTE = mean(nTE)) %>% 
#   filter(position != 1) %>% 
#   ggplot(aes(x = position, y = ctrl_nTE, color = signal_sequence)) +
#   geom_point() +
#   geom_line() +
#   geom_smooth() +
#   scale_y_continuous(limits = c(0,NA)) +
#   scale_color_manual(name = NULL, values = c("#1B9E77", "#D95F02")) +
#   xlab("Codon Position") +
#   ylab("nTE") +
#   ggtitle("Transcriptome") +
#   publication_theme() +
#   theme(legend.position = "none")
# 
# png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/nTEs_SPs_transcriptome.png",
#     res = 300, height = 1000, width = 1500)
# print(nTEs_transcriptome)
# dev.off()
# 
# nTE_CNOT3 <- 
# codonsN_table_long %>% 
#   inner_join(nTEs, by = "codon") %>%
#   ungroup() %>%
#   filter(ENST %in% CNOT3_transcripts) %>% 
#   group_by(signal_sequence, position) %>%
#   summarise(ctrl_nTE = mean(nTE)) %>% 
#   filter(position != 1) %>% 
#   ggplot(aes(x = position, y = ctrl_nTE, color = signal_sequence)) +
#   geom_point() +
#   geom_line() +
#   geom_smooth() +
#   scale_y_continuous(limits = c(0,NA)) +
#   scale_color_manual(name = NULL, values = c("#1B9E77", "#D95F02")) +
#   xlab("Codon Position") +
#   ylab("nTE") +
#   ggtitle("CNOT3 transcripts") +
#   publication_theme() +
#   theme(legend.position = "none")
# 
# png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/nTEs_SPs_CNOT3.png",
#     res = 300, height = 1000, width = 1500)
# print(nTE_CNOT3)
# dev.off()

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
  ggtitle("Homo sapiens") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(0,50)) +
  scale_x_continuous(limits = c(2,50), breaks = c(2,10,20,30,40,50)) +
  publication_theme() +
  theme(legend.position = "bottom", 
        legend.direction = "vertical", 
        legend.title = element_blank()) + guides(colour = guide_legend(nrow = 1))

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/hsapiens.png",
    res = 300,
    height = 1000,
    width = 1500)
print(final_gplot)
dev.off()

get_only_legend <- function(plot) { 
  
  # get tabular interpretation of plot 
  plot_table <- ggplot_gtable(ggplot_build(plot))  
  
  #  Mark only legend in plot 
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")  
  
  # extract legend 
  legend <- plot_table$grobs[[legend_plot]] 
  
  # return legend 
  return(legend)  
}

leg <- get_only_legend(final_gplot)

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/leucines_legend.png",
    res = 300,
    height = 750,
    width = 1000)
grid.newpage()
grid.draw(leg)
dev.off()

HS5CDS <- 
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
  ggtitle("5'CDS") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(0,50)) +
  scale_x_continuous(limits = c(2,50), breaks = c(2,10,20,30,40,50)) +
  publication_theme() +
  theme(legend.position = "none", 
        legend.direction = "vertical", 
        legend.title = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/hs5cds.png",
    res = 300,
    height = 1000,
    width = 1500)
print(HS5CDS)
dev.off()


lines_processed$arginine <- factor(lines_processed$arginine, levels = c("cgg", "cga", "cgt", "cgc", "agg", "aga", "other"))
background_arginine <- lines_processed %>% filter(arginine == "other" & position != 1)
argininel2plot <- lines_processed %>% filter(arginine %in% c("cgg", "cga", "cgt", "cgc", "agg", "aga") & position != 1)
arginine_colours <- c("#D51317FF", "#F39200FF", "#164194FF", "#6F286AFF", "#007B3DFF", "#31B7BCFF")

arginines <- 
ggplot() +
  # geom_line(data = background_arginine,
  #           aes(x = position, y = percent, group = codon),
  #           size = 1,
  #           colour = "#C8C8C8") +
  geom_line(data = argininel2plot,
            aes(x = position, y = percent, group = codon, colour = arginine),
            size = 1) +
  #ylim(0,5) +
  xlab("Codon Position") +
  ylab("Percent") +
  scale_color_manual(values = arginine_colours) +
  ggtitle("5'CDS (Arg codons)") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(0,50)) +
  scale_x_continuous(limits = c(2,50), breaks = c(2,10,20,30,40,50)) +
  publication_theme() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.title = element_blank())+ guides(colour = guide_legend(nrow = 1))

get_only_legend <- function(plot) { 
  
  # get tabular interpretation of plot 
  plot_table <- ggplot_gtable(ggplot_build(plot))  
  
  #  Mark only legend in plot 
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")  
  
  # extract legend 
  legend <- plot_table$grobs[[legend_plot]] 
  
  # return legend 
  return(legend)  
}

leg <- get_only_legend(arginines)

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/arginine_legend.png",
    res = 300,
    height = 750,
    width = 1000)
grid.newpage()
grid.draw(leg)
dev.off()

arginines <- 
  ggplot() +
  # geom_line(data = background_arginine,
  #           aes(x = position, y = percent, group = codon),
  #           size = 1,
  #           colour = "#C8C8C8") +
  geom_line(data = argininel2plot,
            aes(x = position, y = percent, group = codon, colour = arginine),
            size = 1) +
  #ylim(0,5) +
  xlab("Codon Position") +
  ylab("Percent") +
  scale_color_manual(values = arginine_colours) +
  ggtitle("5'CDS (Arg codons)") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(0,50)) +
  scale_x_continuous(limits = c(2,50), breaks = c(2,10,20,30,40,50)) +
  publication_theme() +
  theme(legend.position = "none", 
        legend.direction = "horizontal", 
        legend.title = element_blank())+ guides(colour = guide_legend(nrow = 1))


png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/hs5cds_arg.png",
    res = 300,
    height = 1000,
    width = 1500)
print(arginines)
dev.off()



lines_processed <- 
lines_processed %>% 
  mutate(serine = factor(case_when(codon == "agc" ~ "agc",
                            codon == "agt" ~ "agt",
                            codon == "tcg" ~ "tcg",
                            codon == "tcc" ~ "tcc",
                            codon == "tct" ~ "tct",
                            codon == "tca" ~ "tca",
                            .default = "other"),
                         levels = c("tcg", "tcc", "tca", "tct", "agc", "agt", "other")),
         threonine = factor(case_when(codon == "acg" ~ "acg",
                                     codon == "acc" ~ "acc",
                                     codon == "act" ~ "act",
                                     codon == "aca" ~ "aca",
                                     .default = "other"),
                           levels = c("acg", "acc", "act", "aca", "other")),
         asparagine = factor(case_when(codon == "aac" ~ "aac",
                                       codon == "aat" ~ "aat",
                                       .default = "other"),
                             levels = c("aac", "aat","other")),
         glutamine = factor(case_when(codon == "caa" ~ "caa",
                                       codon == "cag" ~ "cag",
                                       .default = "other"),
                             levels = c("caa", "cag","other")),
         histidine = factor(case_when(codon == "cat" ~ "cat",
                                      codon == "cac" ~ "cac",
                                      .default = "other"),
                            levels = c("cat", "cac","other")),
         lysine = factor(case_when(codon == "aag" ~ "aag",
                                   codon == "aaa" ~ "aaa",
                                   .default = "other"),
                         levels = c("aag", "aaa","other")))




background_lysine = lines_processed %>% filter(lysine == "other")
lysine = lines_processed %>% filter(lysine != "other")


ggplot() +
  geom_line(data = background_lysine,
            aes(x = position, y = percent, group = codon),
            size = 1,
            colour = "#C8C8C8") +
  geom_line(data = lysine,
            aes(x = position, y = percent, group = codon, colour = lysine),
            size = 1) +
  ylim(0,5) +
  xlab("Codon Position") +
  ylab("Percent") +
  scale_color_manual(values = arginine_colours) +
  ggtitle("5'CDS (Lys codons)") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(0,50)) +
  scale_x_continuous(limits = c(2,50), breaks = c(2,10,20,30,40,50)) +
  publication_theme() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.title = element_blank())+ guides(colour = guide_legend(nrow = 1))




background_histidine = lines_processed %>% filter(histidine == "other")
histidine = lines_processed %>% filter(histidine != "other")


ggplot() +
  geom_line(data = background_histidine,
            aes(x = position, y = percent, group = codon),
            size = 1,
            colour = "#C8C8C8") +
  geom_line(data = histidine,
            aes(x = position, y = percent, group = codon, colour = histidine),
            size = 1) +
  ylim(0,5) +
  xlab("Codon Position") +
  ylab("Percent") +
  scale_color_manual(values = arginine_colours) +
  ggtitle("5'CDS (His codons)") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(0,50)) +
  scale_x_continuous(limits = c(2,50), breaks = c(2,10,20,30,40,50)) +
  publication_theme() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.title = element_blank())+ guides(colour = guide_legend(nrow = 1))





background_serine = lines_processed %>% filter(serine == "other")
serine = lines_processed %>% filter(serine != "other")


ggplot() +
  geom_line(data = background_histidine,
            aes(x = position, y = percent, group = codon),
            size = 1,
            colour = "#C8C8C8") +
  geom_line(data = serine,
            aes(x = position, y = percent, group = codon, colour = serine),
            size = 1) +
  ylim(0,5) +
  xlab("Codon Position") +
  ylab("Percent") +
  scale_color_manual(values = arginine_colours) +
  ggtitle("5'CDS (His codons)") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(0,50)) +
  scale_x_continuous(limits = c(2,50), breaks = c(2,10,20,30,40,50)) +
  publication_theme() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.title = element_blank())+ guides(colour = guide_legend(nrow = 1))













# plot the leucyl codon usage for transcripts with and without a signal sequence
lines_gplot_CUG <- 
  lines_processed %>%
  filter(codon %in% leucyl_codons) %>% 
  ggplot(aes(x = position, y = percent, group = codon, color = leucyl)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F", "#999999")) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  ggtitle("Transcriptome: 5'ORF") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(0,50)) +
  publication_theme2() +
  theme(legend.position = "none", 
        legend.direction = "vertical", 
        legend.title = element_blank())

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CUG_human_start.png",
    res = 300, height  = 1000, width = 1850)
print(lines_gplot_CUG)
dev.off()

background_tail <- lines_processed_tail %>% filter(leucyl == "other" & position != 1)
leucyl2plot_tail <- lines_processed_tail %>% filter(leucyl %in% c("ctg", "ctc", "ctt", "ttg", "cta", "tta") & position != 1)

final_gplot_tail <- 
  ggplot() +
  geom_line(data = background_tail,
            aes(x = position, y = percent, group = codon),
            size = 1,
            colour = "#C8C8C8") +
  geom_line(data = leucyl2plot_tail,
            aes(x = position, y = percent, group = codon, colour = leucyl),
            size = 1) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  scale_color_manual(values = leucyl_colours) +
  ggtitle("3'CDS") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(-50,-2)) +
  scale_x_continuous(limits = c(-50,0), breaks = c(-50,-40,-30,-20,-10,-2)) +
  publication_theme() +
  theme(legend.position = "none", 
        legend.direction = "vertical", 
        legend.title = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/HS3CDS.png",
    res = 300, height = 1000, width = 1500)
print(final_gplot_tail)
dev.off()

lines_gplot_tail <- 
  lines_processed_tail %>%
  filter(position!= -1) %>% 
  ggplot(aes(x = position, y = percent, group = codon, color = leucyl)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F", "#999999")) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  ggtitle("Transcriptome: 3'ORF") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(-100,-2)) +
  publication_theme2() +
  theme(legend.position = "none", 
        legend.direction = "vertical", 
        legend.title = element_blank())

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CUG_human_end.png",
    res = 300, height  = 1000, width = 1850)
print(lines_gplot_tail)
dev.off()


ggsci::pal_d3()(5)
"#1F77B4FF"  "#2CA02CFF" "#D62728FF" 

hydrophobic_GC <- 
lines_processed %>%
  filter(CGNCG != "other") %>% 
  ggplot(aes(x = position, y = percent, color = tag)) +
  geom_line(size = 1) +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = c("#FF7F0EFF", "#1F77B4FF"), labels = c("No signal sequence", "Signal sequence")) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  ggtitle("Hydrophobic GC rich codons") +
  facet_wrap(~CGNCG, ncol = 2) +
  coord_cartesian(xlim = c(0,50)) +
  publication_theme() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.title = element_blank())

leg <- get_only_legend(hydrophobic_GC)

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/Signal_sequence_legend.png", res = 300, height = 200, width = 1000)
grid.newpage()
grid.draw(leg)
dev.off()

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/SP_hydrophobic_GC_codons.png",
    res = 300, height  = 2500, width = 1250)
print(hydrophobic_GC)
dev.off()

hydrophobic_AU <- 
lines_processed %>%
  filter(hydrophobic_AU != "other") %>% 
  ggplot(aes(x = position, y = percent, color = tag)) +
  geom_line(size = 1) +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = c("#FF7F0EFF", "#1F77B4FF")) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  ggtitle("Hydrophobic AU rich codons") +
  facet_wrap(~hydrophobic_AU, ncol = 2) +
  coord_cartesian(xlim = c(0,50)) +
  publication_theme() +
  theme(legend.position = "none", 
        legend.direction = "horizontal",
        legend.title = element_blank())

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/SP_hydrophobic_AU_codons.png",
    res = 300, height  = 2750, width = 1250)
print(hydrophobic_AU)
dev.off()



lines_gplot <- 
  lines_processed %>%
  filter(codon %in% arginine_codons) %>% 
  ggplot(aes(x = position, y = percent, group = codon, color = arginine)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F")) +
  ylim(0,5) +
  xlab("Codon Position") +
  ylab("Percent") +
  ggtitle("Human transcripts (Signal P)") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(0,50)) +
  publication_theme() +
  theme(legend.position = "right", 
        legend.direction = "vertical")

lines_processed %>%
  filter(codon %in% ala_codons$codon) %>% 
  ggplot(aes(x = position, y = percent, group = codon, color = alanine)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F")) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  ggtitle("Human transcripts (Signal P)") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(0,50)) +
  publication_theme() +
  theme(legend.position = "right", 
        legend.direction = "vertical")

lines_processed %>%
  filter(codon %in% val_codons$codon) %>% 
  ggplot(aes(x = position, y = percent, group = tag, color = tag)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F")) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  ggtitle("Human transcripts (Signal P)") +
  facet_wrap(~codon) +
  coord_cartesian(xlim = c(0,50)) +
  publication_theme() +
  theme(legend.position = "right", 
        legend.direction = "vertical")

lines_processed %>%
  filter(codon %in% ile_codons$codon) %>% 
  ggplot(aes(x = position, y = percent, group = codon, color = isoleucine)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F")) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  ggtitle("Human transcripts (Signal P)") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(0,50)) +
  publication_theme() +
  theme(legend.position = "right", 
        legend.direction = "vertical")

lines_processed %>%
  filter(codon %in% pro_codons$codon) %>% 
  ggplot(aes(x = position, y = percent, group = codon, color = proline)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F")) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  ggtitle("Human transcripts (Signal P)") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(0,50)) +
  publication_theme() +
  theme(legend.position = "right", 
        legend.direction = "vertical")
  
  
lines_gplot <- 
    lines_processed %>%
    inner_join(CSCs, by = "codon") %>% 
  mutate(optimality = percent * delta_cor_estimate) %>% 
    ungroup() %>% 
    group_by(signal_sequence, position) %>% 
    summarise(total_optimality = sum(optimality)) %>% 
  ggplot(aes(x = position, y = total_optimality, colour = signal_sequence)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F", "#999999")) +
    #ylim(0,20) +
    xlab("Codon Position") +
    #ylab("Percent") +
    ggtitle("Human transcripts") +
    #facet_wrap(~signal_sequence) +
  
    coord_cartesian(xlim = c(0,500)) +
    publication_theme() +
    theme(legend.position = "right", 
          legend.direction = "vertical", 
          legend.title = element_blank()) 


tiff(file.path("CTGenrichment_humans.tiff"), height = 3, width = 6, units = "in", res = 500)
print(lines_gplot_CUG)
dev.off()

only_leu_codons <- 
  lines_processed %>% 
  mutate(Leucine_AA = case_when(codon %in% leucyl_codons ~ T,
                                .default = F)) %>% 
  filter(Leucine_AA == T) %>% 
  ggplot(aes(x = position, y = percent, group = codon, color = leucyl)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F", "#999999")) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  ggtitle("Human transcripts") +
  facet_wrap(~tag) +
  publication_theme() +
  theme(legend.position = "right", 
        legend.direction = "vertical", 
        legend.title = element_blank())

tiff(file.path("CTGenrichment_humans_onlyLeu_codons.tiff"), height = 3, width = 6, units = "in", res = 500)
print(only_leu_codons)
dev.off()

lines_processed <- 
  lines_processed %>% 
  ungroup() %>% 
  inner_join(codons, by= "codon") %>% 
  mutate(box = as.numeric(str_remove(box, "-box"))) 

lines_processed <- 
  lines_processed %>% 
  group_by(signal_sequence, position, AA) %>% 
    summarise(total_codons_per_AA = sum(n)) %>% 
  inner_join(lines_processed, by = c("signal_sequence", "position", "AA")) %>% 
  select(signal_sequence, position, AA, codon, n, total_codons_per_AA, box, wobble, tag, ctg, leucyl)

RSCU <- lines_processed %>% 
  mutate(RSCU = n/total_codons_per_AA * box) %>%
  ungroup() %>% 
  select(position, codon, RSCU, wobble, tag) %>% 
  ungroup() %>% 
  spread(key = tag, value = RSCU) 

RSCU[is.na(RSCU)] <- 0

leucine_codons <- c("ctg", "ctc", "ctt", "cta", "ttg", "tta")

RSCU <- RSCU %>% 
  mutate(leucyl_codons = case_when(codon == "ctg" ~"ctg",
                                   codon == "ctc" ~ "ctc",
                                   codon == "cta" ~ "cta",
                                   codon == "ctt" ~ "ctt",
                                   codon == "ttg" ~ "ttg",
                                   codon == "tta" ~ "tta",
                                   .default = "other"))

RSCU$leucyl_codons <- factor(RSCU$leucyl_codons, levels = c("ctg", "ctc", "ctt", "ttg", "cta", "tta","other"))

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
            #plot.margin=unit(c(,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

RSCU %>% 
  filter(position <= 30) %>% 
  ggplot()

RSCU_30 <- 
RSCU %>% 
  #filter(position != 1) %>% 
  filter(position <= 30) %>% 
  dplyr::rename(`No signal sequence` = `No signal sequence\nn = 15648`,
                `Signal sequence` = `Signal sequence\nn = 2828`) %>% 
  ggplot(aes(x = `No signal sequence`, 
             y = `Signal sequence`, 
             colour = leucyl_codons)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F", "#999999")) +
  #geom_text(show.legend = F) +
  #ggtitle(paste0("Position ", x$position[1])) +
  geom_abline(lty = "dashed") +
  ggtitle("Transcriptome") +
  facet_wrap(~position, ncol = 6, nrow = 5) +
  publication_theme() +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/RSCU_transcriptome_30.png",
    res = 300, height = 2000, width = 2200)
print(RSCU_30)
dev.off()

RSCU_positions_15_30 <- 
RSCU %>% 
  filter(position != 1) %>% 
  filter(position > 15 & position <= 30) %>% 
  ggplot(aes(x = `No signal sequence\nn = 15648`, 
             y = `Signal sequence\nn = 2828`, 
             colour = leucyl_codons)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F", "#999999")) +
  #geom_text(show.legend = F) +
  #ggtitle(paste0("Position ", x$position[1])) +
  geom_abline(lty = "dashed") +
  scale_x_continuous(limits = c(0,5)) +
  scale_y_continuous(limits = c(0,5)) +
  facet_wrap(~position, ncol = 5, nrow = 6) +
  publication_theme() +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/RSCU_positional_15_30.png",
    res = 300, height = 1750, width = 1750)
print(RSCU_positions_15_30)
dev.off()

RSCU <- 
  RSCU %>% 
  dplyr::rename(`No signal sequence` = `No signal sequence\nn = 15648`,
                `Signal sequence` = `Signal sequence\nn = 2828`) %>% 
  mutate(diff = `Signal sequence` - `No signal sequence`)
   
RSCU_diff_list <- split(RSCU, f = RSCU$position)

plot_differentials <- function(x){
  
  y <- 
    x %>% ggplot(aes(reorder(codon, -diff), diff, fill = leucyl_codons)) +
    scale_fill_manual(values = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F", "#999999"))+ 
    geom_col() +
    ggtitle(paste0("Position: ", x$position[1]))+
    publication_theme() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
  
  return(y)
}

RSCU_diff_plots <- lapply(RSCU_diff_list, plot_differentials)  

RSCU_diff_plots[[19]]

RSCU_lyst <- split(RSCU, f = RSCU$position)

plot_RSCU <- function(x){
  
  gplot <- 
    x %>% 
    ggplot(aes(x = `No signal sequence`, y = `Signal sequence`, colour = wobble, label = codon)) +
    geom_point()+
    geom_text(show.legend = F) +
    ggtitle(paste0("Position ", x$position[1])) +
    geom_abline() +
    publication_theme()
  
  return(gplot)
}

my_gplots <- lapply(RSCU_lyst, plot_RSCU)

my_gplots[[19]]

RSCU

test_fun <- function(obj){
  
  sequence <- 1:length(obj)
  lapply(sequence, function(y){paste0(obj[y], obj[(y + 1)])})
}

lapply(my_lyst[[1]], test_fun)

obj <- my_lyst[[1]]
sequence <- 1:length(obj)
unlist(lapply(sequence, function(y){paste0(obj[y], obj[(y + 1)])}))

extract_dicodons <- function(split_cods){
  
  sequence <- 1:length(split_cods)
  fin <- unlist(lapply(sequence, function(y){paste0(split_cods[y], split_cods[(y + 1)])}))
  fin <- head(fin, -1)
  return(fin)
}

dicodons_lyst <- lapply(codonsN, extract_dicodons)

dicodons_lyst[[1]]

leu_arg_dicodons <- merge(leucyl_codons, arginine_codons) %>% 
  mutate(dicodon = paste0(x,y)) %>% 
  pull(dicodon)



count_leuarg_dicodons <- function(x){
  
  final <- unlist(lapply(lapply(dicodons_lyst, str_count, pattern = x), sum))
  names(final) <- NULL
  return(final)
}

test <- lapply(leu_arg_dicodons, count_leuarg_dicodons)
names(test) <- leu_arg_dicodons
leuarg_dicodons <- data.frame(t(do.call("rbind", test)))
rownames(leuarg_dicodons) <- names(dicodons_lyst)
leuarg_dicodons <- leuarg_dicodons %>% 
  rownames_to_column("ENST") %>% 
  gather(key = "dicodon_motif", value = raw_count, cttcgt:ctccgc) %>% 
  inner_join(features %>% select(ENST, signal_sequence)) 

total_motif_counts <- 
  leuarg_dicodons %>%
  group_by(signal_sequence) %>% 
  count()  

leuarg_dicodons %>% 
  group_by(signal_sequence, dicodon_motif) %>% 
  summarise(total = sum(raw_count)) %>% 
  inner_join(total_motif_counts, by = c("signal_sequence"), relationship = "many-to-many") %>% 
  mutate(perc = total/n * 100) %>% 
  ggplot(aes(x = reorder(dicodon_motif, +perc), y = perc)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~signal_sequence) +
  xlab("Dicodon")

test$ctgcgc


