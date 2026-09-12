library(tidyverse)
library(seqinr)

# home
home <- "\\\\data.beatson.gla.ac.uk/data"

# read in master
master <- read_csv(file.path(home, "R11/James/sequences/master.csv"))

# select the longest transcript per gene
longest <- master %>% 
  mutate(polypeptide_length = nchar(master$amino_acid_sequence)) %>%
  group_by(ENSG) %>%
  mutate(the_rank  = rank(-polypeptide_length, ties.method = "first")) %>%
  filter(the_rank == 1) %>% 
  select(-the_rank)

over_100_codons <- longest %>% filter(nchar(nucleotide_sequence_CDS) > 300)

features <- read_csv(file.path(home, "R11/James/sequences/features.csv")) %>% select(ENSG,ENST, signal_sequence)

signalP <- 
  inner_join(over_100_codons %>% select(ENSG,ENST,nucleotide_sequence_CDS),
             features, by = c("ENSG", "ENST")) %>% 
  mutate(nucleotide_sequence_CDS = substr(nucleotide_sequence_CDS, 1, 90))

signalP_SP <- 
  signalP %>% 
  filter(signal_sequence == T)

signalP_noSP <- 
  signalP %>% 
  filter(signal_sequence == F)

#lapply(1:length(signalP_SP$ENST)
RSCU_signalP_SP <- 
  lapply(
    signalP_SP %>% 
      pull(nucleotide_sequence_CDS) %>% 
      str_split(pattern = ""),
    uco, 
    as.data.frame = T, 
    NA.rscu = 0)

RSCU_signalP_SP <- 
  do.call("rbind", lapply(1:length(signalP_SP$ENST), function(x){
    
    data <- RSCU_signalP_SP[[x]]
    transcript <- signalP_SP$ENST[[x]]
    data$transcript <- transcript
    
    return(data)
    
  }))

#lapply(1:length(signalP_SP$ENST)
RSCU_signalP_noSP <- 
  lapply(
    signalP_noSP %>% 
      pull(nucleotide_sequence_CDS) %>% 
      str_split(pattern = ""),
    uco, 
    as.data.frame = T, 
    NA.rscu = 0)

RSCU_signalP_noSP <- 
  do.call("rbind", lapply(1:length(signalP_noSP$ENST), function(x){
    
    data <- RSCU_signalP_noSP[[x]]
    transcript <- signalP_noSP$ENST[[x]]
    data$transcript <- transcript
    
    return(data)
    
  }))

SigPFinal <- 
  rbind(
    RSCU_signalP_SP %>% mutate(transcript_class = "signal_peptide"),
    RSCU_signalP_noSP %>% mutate(transcript_class = "no_signal_peptide")) 

SigPFinal %>% 
  filter(codon == "ctc") %>% 
  ggplot(aes(x = transcript_class, y = RSCU)) +
  geom_violin() +
  geom_boxplot(width = 0.1)
#####  

tmhmm <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/TM_domains_FINAL.csv") %>% 
  select(ENST:codon_position)

tmhmm$class %>% unique()

tmhmm_filt <- tmhmm %>% filter(ENST %in% longest$ENST)

T1_transcripts <- tmhmm_filt %>% filter(class %in% c("type1")) %>% pull(ENST) %>% unique()
SPonly_transcripts <- tmhmm_filt %>% filter(class %in% c("sponly")) %>% pull(ENST) %>% unique()
multipassT1 <- tmhmm_filt %>% filter(class %in% c("multipassT1")) %>% pull(ENST) %>% unique()
multipassT2 <- tmhmm_filt %>% filter(class %in% c("multipassT2")) %>% pull(ENST) %>% unique()
multipassT3 <- tmhmm_filt %>% filter(class %in% c("multipassT3")) %>% pull(ENST) %>% unique()
type2 <- tmhmm_filt %>% filter(class %in% c("type2")) %>% pull(ENST) %>% unique()
type3 <- tmhmm_filt %>% filter(class %in% c("type3")) %>% pull(ENST) %>% unique()
TA <- tmhmm_filt %>% filter(class %in% c("TA")) %>% pull(ENST) %>% unique()

# unique(SigPFinal$transcript) %>% length()
# over_100_codons$ENST %>% unique() %>% length()
#  
# SigPFinal %>% filter(transcript_class == "no_signal_peptide") %>% pull(transcript) %>% unique() %>% length()

SigPvstmhmm <- 
  SigPFinal %>% mutate(tmhmm = case_when(transcript %in% T1_transcripts~ "type1",
                                         transcript %in% SPonly_transcripts~ "sponly",
                                         transcript %in% multipassT1 ~"multipassT1",
                                         transcript %in% multipassT2 ~"multipassT2",
                                         transcript %in% multipassT3 ~"multipassT3",
                                         transcript %in% type2 ~ "type2",
                                         transcript %in% type3 ~ "type3",
                                         transcript %in% TA ~ "TA",
                                         .default = "cytoplasmic")) 

SigPvstmhmm %>% 
  select(transcript_class, tmhmm, transcript) %>% unique() %>% 
  group_by(transcript_class, tmhmm) %>% 
  summarise(n_count = n())

SigPvstmhmm %>% 
  select(transcript_class, tmhmm, transcript) %>% unique() %>% 
  group_by( tmhmm) %>% 
  summarise(n_count = n())

transcripts_to_consider <- 
  SigPvstmhmm %>% 
  filter(transcript_class == "no_signal_peptide" & tmhmm == "TA" |
           transcript_class == "no_signal_peptide" & tmhmm == "cytoplasmic"|
           transcript_class == "no_signal_peptide" & tmhmm == "multipassT2"|
           transcript_class == "no_signal_peptide" & tmhmm == "multipassT3"|
           transcript_class == "no_signal_peptide" & tmhmm == "type2"|
           transcript_class == "no_signal_peptide" & tmhmm == "type3" | 
           transcript_class == "signal_peptide" & tmhmm == "multipassT1" |
           transcript_class == "signal_peptide" & tmhmm == "sponly" |
           transcript_class == "signal_peptide" & tmhmm == "type1") %>% 
  pull(transcript) %>% unique() 

tmhmm %>% 
  filter(ENST %in% transcripts_to_consider) %>% 
  filter(region == "TS" & ID != 0) %>% 
  pull(ENST) %>% unique() %>% length()

tmhmm_refiltered <- 
  tmhmm %>% 
  filter(ENST %in% transcripts_to_consider) %>% 
  filter(region == "TS") %>% 
  mutate(tag = paste0(ENST, "_", ID)) %>% 
  inner_join(codons, by = "codon")

tmhmm_refiltered_split <- split(tmhmm_refiltered, tmhmm_refiltered$tag)

calculate_RSCU <- function(x){
  
  codon_freq <- as.data.frame(table(x$codon))
  colnames(codon_freq) <- c("codon", "codon_freq")
  
  int <- inner_join(codon_freq, codons, by = "codon")
  
  final_df <- 
    int %>% 
    group_by(AA) %>% 
    summarise(total_residues = sum(codon_freq)) %>% 
    inner_join(int, by = "AA") %>% 
    mutate(RSCU = codon_freq/total_residues * box) %>% 
    mutate(transcript = unique(x$ENST),
           ID = unique(x$ID),
           class = unique(x$class))
  
  return(final_df)
  
}

test <- lapply(tmhmm_refiltered_split, calculate_RSCU)

transcripts_already_considered <- 
  c(do.call("rbind", test) %>% pull(transcript) %>% unique())

transcripts_to_consider[!transcripts_to_consider %in% transcripts_already_considered]

do.call("rbind", test) %>% 
  mutate(id = paste0(transcript, "_", ID)) %>% 
  mutate(new_classification = case_when(class %in% c("multipassT1", "multipassT2", "multipassT3", "type1", "type2", "type3") & ID != 0 ~ "TM domains",
                                        ID == 0 ~  "SP",
                                        class == "TA" ~ "TA")) %>% 
  filter(codon == "ctc")-> TM_info #%>% 

# ggplot(aes(x = new_classification, y = RSCU)) +
#   geom_violin() +
#   geom_boxplot(width = 0.1)

# controls
cytoplasmic_transcripts <- 
  SigPvstmhmm %>% 
  filter(transcript_class == "no_signal_peptide" & tmhmm == "cytoplasmic")%>% 
  pull(transcript) %>% unique() 

cytoplasmics_filtered <- master %>% 
  filter(ENST %in% cytoplasmic_transcripts)

sequences <- str_split(cytoplasmics_filtered %>% 
                         pull(nucleotide_sequence_CDS), "")

cytoplasmic_transcripts_RSCU <- 
  lapply(1:length(sequences), function(x){
    
    a_sequence <- sequences[[x]][1:90]
    transcript <- cytoplasmics_filtered$ENST[[x]]
    df1 <- 
      uco(a_sequence, as.data.frame = T, NA.rscu = 0) %>% 
      mutate(transcript = transcript,
             part_of_transcript = "30codons")
    
    b_sequence <- sequences[[x]][91:length(sequences[[x]])]
    df2 <- 
      uco(b_sequence, as.data.frame = T, NA.rscu = 0) %>% 
      mutate(transcript = transcript,
             part_of_transcript = "RoT")
    return(rbind(df1,df2))
  })

cytoplasmic_information <- do.call("rbind", cytoplasmic_transcripts_RSCU) %>% filter(codon == "ctc")
#TM_info




####
tmhmm_non_coding <- 
  tmhmm %>% 
  filter(ENST %in% transcripts_to_consider) %>% 
  filter(region != "TS") %>% 
  mutate(tag = paste0(ENST, "_", ID)) %>%
  filter(rel_position > -40 & rel_position < 40) %>% 
  inner_join(codons, by = "codon")

tmhmm_non_coding_split <- split(tmhmm_non_coding, tmhmm_non_coding$tag)

test2 <- lapply(tmhmm_non_coding_split, calculate_RSCU)

do.call("rbind", test2) %>% 
  mutate(id = paste0(transcript, "_", ID)) %>% 
  mutate(new_classification = "Non_TM") %>% 
  filter(codon == "ctc")-> non_TM_info



final_data <- 
  rbind(rbind(cytoplasmic_information %>% select(transcript, 
                                                 new_classification = part_of_transcript,
                                                 RSCU),
              TM_info %>% select(transcript, new_classification, RSCU)),
        non_TM_info %>% select(transcript, new_classification, RSCU))  

final_data$new_classification <- factor(final_data$new_classification,
                                        levels = c("30codons",
                                                   "RoT",
                                                   "SP",
                                                   "TM domains",
                                                   "TA",
                                                   "Non_TM"),
                                        labels = c("Codons\n0-30",
                                                   "Rest\nof mRNA",
                                                   "SP",
                                                   "TM\ndomains",
                                                   "TM\ndomains", 
                                                   "Non\nSP/TM"))
final_data <- 
  final_data %>% 
  mutate(`mRNA type` = case_when(new_classification %in% c("Codons\n0-30", "Rest\nof mRNA") ~ "Cytoplasmic",
                                 .default = "ER targeted"))

final_data$`mRNA type` <- factor(final_data$`mRNA type`, levels = c("Cytoplasmic", "ER targeted"))

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
            axis.text.y = element_text(size = 16, color = "black"),
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

gplot <- 
  ggplot(final_data, aes(x = new_classification, y = RSCU, fill = `mRNA type`)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  xlab("Region of mRNA") +
  ylab("CTC RSCU") +
  stat_summary(fun.y=mean, geom="point", shape=15, size=1.5, color="black", fill="black") +
  scale_fill_manual(values = c("#FF7F0EFF", "#1F77B4FF")) +
  publication_theme() + theme(legend.title = element_blank())

TukeyHSD(aov(RSCU ~ new_classification, final_data))

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CTC_RSCU_all_regions.png",
    res = 300, height = 1000, width = 1750)
print(gplot)
dev.off()

####
SigPvstmhmm %>% 
  filter(codon == "ctg") %>% 
  mutate(tmhmm_SP = case_when(tmhmm %in% c("type1", "sponly") ~ "SP",
                              .default = "cytoplasmic")) %>% 
  ggplot(aes(x = tmhmm_SP, y = RSCU)) +
  geom_violin() +
  geom_boxplot(width = 0.1)

SigPvstmhmm %>% 
  filter(tmhmm == "type1" & transcript_class == "signal_peptide" | 
           tmhmm == "sponly" & transcript_class == "signal_peptide" | 
           tmhmm == "cytoplasmic" & transcript_class== "no_signal_peptide") %>% 
  filter(codon == "ctg") %>% 
  ggplot(aes(x = tmhmm, y = RSCU)) +
  geom_violin() +
  geom_boxplot(width = 0.1)

ambigious_transcripts <- 
  SigPvstmhmm %>% 
  filter(transcript_class == "signal_peptide" & tmhmm == "cytoplasmic") %>% pull(transcript) %>% unique()


SigPvstmhmm %>% filter(!transcript %in% ambigious_transcripts) %>% 
  filter(tmhmm == "sponly") %>% pull(transcript_class) %>% table()
