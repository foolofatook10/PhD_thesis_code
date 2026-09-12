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


# read in features table
features <- as_tibble(fread("\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/features.csv", drop = "V1")) %>%
  select(ENSG, ENST, signal_sequence)

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

#combined

#combined %>% filter(signal_sequence == F) %>% pull(ENST)
#cytoplasmic_transcripts_tmhmm

# make a list with the nucleotide sequences chopped into codons
my_lyst <- sapply(combined$nucleotide_sequence_CDS, str_split, pattern = "")

# Name each element of the list by its ENST tag
names(my_lyst) <- pull(combined, var = ENST)

# obtain logical vector of CDSs that are greater than 100 codons (or 300 nts)
over_100_codons <- lapply(my_lyst, length) > 300

# filter my_lyst for CDSs that are greater than 100 codons
my_lyst <-  my_lyst[c(over_100_codons)]

# only consider the first 50 codons per sequence
codons_head <-lapply(my_lyst, head, n=90)

# rest of the transcript
codons_tail <- lapply(my_lyst, tail, n = -90)



RSCUs <- lapply(codons_head, seqinr::uco, as.data.frame = T, NA.rscu = 0)

signal_sequence_table <- do.call("rbind", RSCUs)

signal_sequence_table <- 
signal_sequence_table %>% 
  rownames_to_column("transcript") %>% 
  mutate(transcript = str_remove(transcript, pattern = ".[:lower:]{3}$")) %>% 
  inner_join(features, by = c("transcript" = "ENST"))



# ### quick check ####
# 
# transcripts_to_consider <-   
# unique(signal_sequence_table$transcript)
# 
# sequences <- master %>% filter(ENST %in% transcripts_to_consider) %>% pull(nucleotide_sequence_CDS)
# names <- master %>% filter(ENST %in% transcripts_to_consider) %>% pull(ENST)
#   
# sequences_split <- lapply(sequences, split_into_codons)
# names(sequences_split) <- names
# 
# temp <- lapply(lapply(sequences_split, function(x){x[1:30]}), function(x){
#   
#   y <- as.data.frame(table(x))
#   return(y)
#   
# })
# 
# temp2 <- lapply(1:length(names), function(x){
#   
#   transcript <- names[[x]]
#   data <- temp[[x]]
#   
#   return(
#     data %>% 
#     mutate(transcript = transcript) %>%
#     dplyr::rename(codon = x, freequency = Freq)
#       )
#   
# })
# 
# codons <- 
# codons %>% 
#   mutate(box = case_when(AA == "Ile" ~ "3-box",
#                          .default = box))
# 
# data_temp <- 
# do.call("rbind", temp2) %>% 
#   inner_join(codons, by = "codon") %>% 
#   mutate(box = as.numeric(str_remove(box, "-box")))
# 
# AA_count <- 
# data_temp %>% 
#   group_by(transcript, AA) %>% 
#   summarise(sum_AA = sum(freequency))
# 
# data_temp <-  
# data_temp %>% 
#   inner_join(AA_count, by = c("transcript", "AA")) %>% 
#   mutate(RSCU2 = freequency/sum_AA * box) %>% 
#   select(transcript, codon, RSCU2)
#  
# comparison <- 
# inner_join(
# signal_sequence_table,
# data_temp,
# by = c("transcript", "codon"))
# 
# comparison[!c(round(comparison$RSCU, 5) == round(comparison$RSCU2, 5)),]

#RSCUs$ENST00000216962.9

#### 

signal_sequence_table_filt <- 
signal_sequence_table %>% 
  filter(!codon %in% c("taa", "tga", "tag"))

signal_sequence_table_filt$codon <- droplevels(signal_sequence_table_filt$codon)

signal_sequence_table_filt <- 
  signal_sequence_table_filt %>% 
  mutate(signal_sequence = factor(case_when(signal_sequence == F ~ "no",
                                     .default = "yes"),
                                  levels = c("yes", "no")))

signal_sequence_table_split <- split(signal_sequence_table_filt, signal_sequence_table_filt$codon)

conf_int_lyst = list()


for(i in 1:length(signal_sequence_table_split)){
  
  dataframee <- signal_sequence_table_split[[i]]
  
  AA <- unique(dataframee$AA)
  codon <- unique(dataframee$codon)
  
  print(codon)
  
  
  stat_test <- wilcox.test(RSCU ~ signal_sequence, data = dataframee, conf.int = T)
  conf_int_upper <- stat_test$conf.int[2]
  conf_int_lower <- stat_test$conf.int[1]
  p_val <- stat_test$p.value
  #print(conf_int_upper)
  #print(conf_int_lower)
  
  
  
  df <- 
  data.frame(AA = AA,
             codon = codon,
             conf_int_upper = conf_int_upper,
             conf_int_lower = conf_int_lower,
             p_val = p_val)
  
  conf_int_lyst[[paste0(codon)]] <- df
  
}

codon_data <- do.call("rbind", conf_int_lyst)

colours <- c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F", "#999999")

colour_vec <- 
case_when(codon_data$codon == "ctg" ~ "#377EB8",
          codon_data$codon == "ctc" ~ "#984EA3",
          codon_data$codon == "ctt" ~ "#4DAF4A",
          codon_data$codon == "cta" ~ "#FF7F00",
          codon_data$codon == "ttg" ~ "#BF5B17",
          codon_data$codon == "tta" ~ "#F0027F",
          .default = "#999999")

codon_data <- 
codon_data %>% 
  mutate(signif = case_when(p_val < 0.05 & p_val > 0.01 ~ "*",
                            p_val < 0.01 & p_val > 0.001 ~ "**",
                            p_val < 0.001 ~ "***",
                            .default = ""))
 
CI_plot <- 
ggplot(data = codon_data,
       aes(x = codon, label = signif)) +
  geom_linerange(aes(ymin = conf_int_lower,
                     ymax = conf_int_upper, x = codon ), size = 1) +
  geom_point(aes(x = codon, y = conf_int_upper), colour = "red", alpha = 0.5)+
  geom_point(aes(x = codon, y = conf_int_lower), alpha = 0.5) +
  geom_text(aes(y = 0.8), size = 6, colour = "black", vjust = 0.75) +
  xlab("Codon") +
  coord_flip() +
  publication_theme() +
  theme(axis.text.y = element_text(colour = colour_vec, size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CI_plot.png",
    height = 2500, width = 1250, res = 300)
print(CI_plot)
dev.off()


CTG_RSCU_per_transcript <- 
ggplot(signal_sequence_table %>% filter(codon == "ctg"), 
       aes(x = signal_sequence, y = RSCU, fill = signal_sequence)) +
  geom_violin() +
  geom_boxplot(width = 0.05)+
  ggtitle("CTG") +
  scale_fill_manual(values = c("#FF7F0EFF", "#1F77B4FF"), 
                    labels = c("No signal sequence", "Signal sequence")) +
  scale_x_discrete(labels = c("Transcripts\nwithout SP", "Transcripts\nwith SP")) +
  stat_summary(fun.y=mean, geom="point", shape=15, size=1, color="black", fill="black") +
  ylab("RSCU (Codons 0-30)") +
  publication_theme() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20))

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CTG_RSCU30.png",
    height = 1000, width = 1000, res = 300)
print(CTG_RSCU_per_transcript)
dev.off()

CTC_RSCU_per_transcript <- 
ggplot(signal_sequence_table %>% filter(codon == "ctc"), 
       aes(x = signal_sequence, y = RSCU, fill = signal_sequence)) +
  geom_violin() +
  geom_boxplot(width = 0.05)+
  ggtitle("CTC") +
  scale_fill_manual(values = c("#FF7F0EFF", "#1F77B4FF"), 
                    labels = c("No signal sequence", "Signal sequence")) +
  scale_x_discrete(labels = c("Transcripts\nwithout SP", "Transcripts\nwith SP")) +
  stat_summary(fun.y=mean, geom="point", shape=15, size=1, color="black", fill="black") +
  ylab("RSCU (Codons 0-30)") +
  publication_theme() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20))

ggplot(signal_sequence_table %>% filter(codon == "tca"), 
       aes(x = signal_sequence, y = RSCU, fill = signal_sequence)) +
  geom_violin() +
  geom_boxplot(width = 0.05)+
  ggtitle("CTC") +
  scale_fill_manual(values = c("#FF7F0EFF", "#1F77B4FF"), 
                    labels = c("No signal sequence", "Signal sequence")) +
  scale_x_discrete(labels = c("Transcripts\nwithout SP", "Transcripts\nwith SP")) +
  stat_summary(fun.y=mean, geom="point", shape=15, size=1, color="black", fill="black") +
  ylab("RSCU (Codons 0-30)") +
  publication_theme() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20))

cta_data <- signal_sequence_table %>% filter(codon == "cta")


wilcox.test(RSCU ~ signal_sequence, data = signal_sequence_table %>% filter(codon == "ctg"), conf.int = T)
wilcox.test(RSCU ~ signal_sequence, data = signal_sequence_table %>% filter(codon == "ctc"), conf.int = T)

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CTC_RSCU30.png",
    height = 1000, width = 1000, res = 300)
print(CTC_RSCU_per_transcript)
dev.off()

RSCU_info <- 
signal_sequence_table %>% 
  group_by(signal_sequence, codon) %>% 
  summarise(mean_RSCU = mean(RSCU)) %>% 
  spread(key = "signal_sequence", value = mean_RSCU) %>% 
  dplyr::rename(No_SP = "FALSE",
         SP = "TRUE") %>%
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A/U", "A/U", "G/C","G/C"))) %>% 
  mutate(leucyl_codons = case_when(codon == "ctg" ~"ctg",
                                   codon == "ctc" ~ "ctc",
                                   codon == "cta" ~ "cta",
                                   codon == "ctt" ~ "ctt",
                                   codon == "ttg" ~ "ttg",
                                   codon == "tta" ~ "tta",
                                   .default = "other"),
         serine_codons = case_when(codon == "tct" ~"ctg",
                                   codon == "ctc" ~ "ctc",
                                   codon == "cta" ~ "cta",
                                   codon == "ctt" ~ "ctt",
                                   codon == "ttg" ~ "ttg",
                                   codon == "tta" ~ "tta",
                                   .default = "other"))


RSCU_info <- RSCU_info %>% filter(!codon %in% c("tga", "taa", "tag"))

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

# make box column a numeric
codons <- codons %>% mutate(box = as.numeric(str_extract(box, pattern = "\\d")))

AAs <- unique(codons$AA)


RSCU_info$leucyl_codons <- factor(RSCU_info$leucyl_codons, levels = c("other", "ctg", "ctc", "ctt", "ttg", "cta", "tta"))

RSCU_gplot <- 
RSCU_info %>% 
  ggplot(aes(x = No_SP, y = SP, colour = leucyl_codons)) + 
  geom_point()+
  scale_color_manual(name = NULL,
                     values = c("#999999", 
                                "#377EB8", 
                                                          "#984EA3",
                                                          "#4DAF4A", 
                                                          "#BF5B17", 
                                                          "#FF7F00", 
                                                          "#F0027F")) +
  geom_point(shape = 19, alpha = 0.5, size = 3) +
  #geom_text(show.legend = F) +
  geom_abline(slope = 1, intercept = 0, linewidth = 1, linetype = "dashed") +
  xlab("No signal peptide") +
  ylab("Signal peptide") +
  ggtitle("Codons 0-30 RSCU") +
  publication_theme()  +
  theme(legend.position = "right",
        legend.direction = "vertical")

#
png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/bioinformatic_analysis/plots/RSCU_LEU_Codons.png",
     height = 1200, width = 2000, res = 300)
print(RSCU_gplot)
dev.off()

aminos <- split(codons, f= codons$AA)

my_gplot_lyst <- list()

for (amino_acid in AAs) {
  
  y <- aminos[[amino_acid]]
  
  box <- unique(y$box)
  
  if(box == 6){
  
   z <- RSCU_info %>% 
     mutate(amino = case_when(codon == y$codon[1] ~ y$codon[1],
                              codon == y$codon[2] ~ y$codon[2],
                              codon == y$codon[3] ~ y$codon[3],
                              codon == y$codon[4] ~ y$codon[4],
                              codon == y$codon[5] ~ y$codon[5],
                              codon == y$codon[6] ~ y$codon[6],
                              .default = "other"))
   z$amino <- factor(z$amino, levels = c("other", y$codon[1],y$codon[2],y$codon[3],
                                         y$codon[4],y$codon[5],y$codon[6]))
   
  } else if(box == 4){
    
    z <- RSCU_info %>% 
      mutate(amino = case_when(codon == y$codon[1] ~ y$codon[1],
                               codon == y$codon[2] ~ y$codon[2],
                               codon == y$codon[3] ~ y$codon[3],
                               codon == y$codon[4] ~ y$codon[4],
                               .default = "other"))
    
    z$amino <- factor(z$amino, levels = c("other", y$codon[1],y$codon[2],y$codon[3],
                                          y$codon[4]))
    
  } else if(box == 2){
    
    z <- RSCU_info %>% 
      mutate(amino = case_when(codon == y$codon[1] ~ y$codon[1],
                               codon == y$codon[2] ~ y$codon[2],
                               .default = "other"))
    
    z$amino <- factor(z$amino, levels = c("other", y$codon[1],y$codon[2]))
    
  } else{
    
    z <- RSCU_info %>% 
      mutate(amino = case_when(codon == y$codon[1] ~ y$codon[1],
                               .default = "other"))
    
    z$amino <- factor(z$amino, levels = c("other", y$codon[1]))
  }
  
  gplot <- ggplot(data = z, aes(x = No_SP, 
                               y = SP, 
                               group = codon, 
                               colour = amino)) +
    scale_color_manual(name = paste(amino_acid), values = c("#999999",
                                                     "#377EB8", 
                                                     "#984EA3",
                                                     "#4DAF4A", 
                                                     "#BF5B17", 
                                                     "#FF7F00", 
                                                     "#F0027F")) +
    geom_point(shape = 19, alpha = 0.5, size = 3) +
    #geom_text(show.legend = F) +
    geom_abline(slope = 1, intercept = 0, linewidth = 1, linetype = "dashed") +
    xlab("No signal peptide") +
    ylab("Signal peptide") +
    ggtitle(paste(amino_acid)) +
    publication_theme() +
    theme(legend.position = "right",
          legend.direction = "vertical") 
  
  my_gplot_lyst[[amino_acid]] <- gplot 
}

AAs1 <- 
cowplot::plot_grid(my_gplot_lyst[[1]], 
                   my_gplot_lyst[[2]], 
                   my_gplot_lyst[[3]],
                   my_gplot_lyst[[4]],
                   ncol = 2,
                   nrow = 2)

AAs2 <- 
cowplot::plot_grid( 
                   my_gplot_lyst[[5]], 
                   my_gplot_lyst[[6]],
                   my_gplot_lyst[[7]], 
                   my_gplot_lyst[[8]],
                   ncol = 2,
                   nrow = 2)
AAs3 <- 
cowplot::plot_grid(my_gplot_lyst[[9]],
                   my_gplot_lyst[[10]], 
                   my_gplot_lyst[[11]], 
                   my_gplot_lyst[[12]],
                   ncol = 2,
                   nrow = 2)
AAs4 <- 
cowplot::plot_grid(
                   my_gplot_lyst[[13]], 
                   my_gplot_lyst[[14]], 
                   my_gplot_lyst[[15]],
                   my_gplot_lyst[[16]],
                   ncol = 2,
                   nrow = 2)
AAs5 <- 
cowplot::plot_grid( 
                   my_gplot_lyst[[17]], 
                   my_gplot_lyst[[18]],
                   my_gplot_lyst[[19]],
                   my_gplot_lyst[[20]],
                   ncol = 2,
                   nrow = 2)

png("RSCU1.png", height = 2000, width = 2500, res = 300)
print(AAs1)
dev.off()

png("RSCU2.png", height = 2000, width = 2500, res = 300)
print(AAs2)
dev.off()

png("RSCU3.png", height = 2000, width = 2500, res = 300)
print(AAs3)
dev.off()

png("RSCU4.png", height = 2000, width = 2500, res = 300)
print(AAs4)
dev.off()

png("RSCU5.png", height = 2000, width = 2500, res = 300)
print(AAs5)
dev.off()


# RSCU_info %>% 
#   mutate(delta = SP/No_SP) %>% 
#   arrange(-delta) %>% 
#   ggplot(aes(x = reorder(codon, -delta), y = delta, fill = leucyl_codons)) +
#   geom_col()
