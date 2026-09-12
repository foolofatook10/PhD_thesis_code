library(data.table)
library(tidyverse)

#### Functions ####
# define function
extract_counts <- function(x){
  
  y <- str_replace_all(x$counts, "[\\[\\]\\s+\\'\\']", "")
  
  z <- lapply(str_split(y, ","), as.numeric)
  
  names(z) <- x$ID
  
  return(z)
  
}

extract_sequences <- function(x){
  
  seqs_lyst <- str_split(x$sequence, "")
  names(seqs_lyst) <- x$ID
  
  return(seqs_lyst)
  
}


splice_cds <- function(x, alyst){
  
  ID = transcripts[x]
  df_subsetted = counts_table[x,]
  
  fputr_length = df_subsetted$fputr_length
  cds_length = df_subsetted$cds_length
  tputr_length = df_subsetted$tputr_length
  
  cds_start = fputr_length + 1
  cds_end = fputr_length + cds_length
  cds_length = cds_end - (cds_start -1)
  
  spliced_lyst <- alyst[[{{ID}}]][cds_start:cds_end]
  
  return(spliced_lyst)
  
  
}


count_codon_occupancy <- function(int){
  
  #print(paste("% complete:", round(int/length(transcripts) * 100, digits = 2)))
  
  ID = transcripts[[int]]
  print(int)
  
  start = 10
  end =  length(spliced_counts[[{{ID}}]]) - 11
  
  series <- seq(start, end, 3)
  
  
  minus4_lyst = list("tct" = 0, "tca" = 0, "tcg" = 0, "agc" = 0, "agt" = 0, "tcc" = 0, "cgt" = 0, "cga" = 0, "cgg" = 0, "aga" = 0, 
                     "agg" = 0, "cgc" = 0, "ctt" = 0, "cta" = 0, "ctg" = 0, "tta" = 0, "ttg" = 0, "ctc" = 0, "gct" = 0, "gca" = 0, 
                     "gcg" = 0, "gcc" = 0, "cct" = 0, "cca" = 0, "ccg" = 0,"ccc" = 0, "gtt" = 0, "gta" = 0, "gtg" = 0,"gtc" = 0, 
                     "act" = 0, "aca" = 0, "acg" = 0, "acc" = 0, "ggc" = 0, "gga" = 0, "ggg" = 0, "ggt" = 0, "att" = 0, "atc" = 0, 
                     "ata" = 0, "ttc" = 0, "ttt" = 0, "aac" = 0, "aat" = 0, "cac" = 0, "cat" = 0, "gac" = 0, "gat" = 0, "tac" = 0, 
                     "tat" = 0, "tgc" = 0, "tgt" = 0, "gaa" = 0, "gag" = 0, "caa" = 0, "cag" = 0, "aaa" = 0, "aag" = 0, "tgg" = 0, 
                     "atg" = 0, "tag" = 0, "taa" = 0, "tga" = 0)
  
  minus3_lyst = list("tct" = 0, "tca" = 0, "tcg" = 0, "agc" = 0, "agt" = 0, "tcc" = 0, "cgt" = 0, "cga" = 0, "cgg" = 0, "aga" = 0, 
                     "agg" = 0, "cgc" = 0, "ctt" = 0, "cta" = 0, "ctg" = 0, "tta" = 0, "ttg" = 0, "ctc" = 0, "gct" = 0, "gca" = 0, 
                     "gcg" = 0, "gcc" = 0, "cct" = 0, "cca" = 0, "ccg" = 0,"ccc" = 0, "gtt" = 0, "gta" = 0, "gtg" = 0,"gtc" = 0, 
                     "act" = 0, "aca" = 0, "acg" = 0, "acc" = 0, "ggc" = 0, "gga" = 0, "ggg" = 0, "ggt" = 0, "att" = 0, "atc" = 0, 
                     "ata" = 0, "ttc" = 0, "ttt" = 0, "aac" = 0, "aat" = 0, "cac" = 0, "cat" = 0, "gac" = 0, "gat" = 0, "tac" = 0, 
                     "tat" = 0, "tgc" = 0, "tgt" = 0, "gaa" = 0, "gag" = 0, "caa" = 0, "cag" = 0, "aaa" = 0, "aag" = 0, "tgg" = 0, 
                     "atg" = 0, "tag" = 0, "taa" = 0, "tga" = 0)
  
  E_lyst = list("tct" = 0, "tca" = 0, "tcg" = 0, "agc" = 0, "agt" = 0, "tcc" = 0, "cgt" = 0, "cga" = 0, "cgg" = 0, "aga" = 0, 
                "agg" = 0, "cgc" = 0, "ctt" = 0, "cta" = 0, "ctg" = 0, "tta" = 0, "ttg" = 0, "ctc" = 0, "gct" = 0, "gca" = 0, 
                "gcg" = 0, "gcc" = 0, "cct" = 0, "cca" = 0, "ccg" = 0,"ccc" = 0, "gtt" = 0, "gta" = 0, "gtg" = 0,"gtc" = 0, "act" = 0, 
                "aca" = 0, "acg" = 0, "acc" = 0, "ggc" = 0, "gga" = 0, "ggg" = 0, "ggt" = 0, "att" = 0, "atc" = 0, "ata" = 0, 
                "ttc" = 0, "ttt" = 0, "aac" = 0, "aat" = 0, "cac" = 0, "cat" = 0, "gac" = 0, "gat" = 0, "tac" = 0, "tat" = 0, 
                "tgc" = 0, "tgt" = 0, "gaa" = 0, "gag" = 0, "caa" = 0, "cag" = 0, "aaa" = 0, "aag" = 0, "tgg" = 0, "atg" = 0, "tag" = 0, "taa" = 0, "tga" = 0)
  
  P_lyst = list("tct" = 0, "tca" = 0, "tcg" = 0, "agc" = 0, "agt" = 0, "tcc" = 0, "cgt" = 0, "cga" = 0, "cgg" = 0, "aga" = 0, 
                "agg" = 0, "cgc" = 0, "ctt" = 0, "cta" = 0, "ctg" = 0, "tta" = 0, "ttg" = 0, "ctc" = 0, "gct" = 0, "gca" = 0, 
                "gcg" = 0, "gcc" = 0, "cct" = 0, "cca" = 0, "ccg" = 0,"ccc" = 0, "gtt" = 0, "gta" = 0, "gtg" = 0,"gtc" = 0, 
                "act" = 0, "aca" = 0, "acg" = 0, "acc" = 0, "ggc" = 0, "gga" = 0, "ggg" = 0, "ggt" = 0, "att" = 0, "atc" = 0, 
                "ata" = 0, "ttc" = 0, "ttt" = 0, "aac" = 0, "aat" = 0, "cac" = 0, "cat" = 0, "gac" = 0, "gat" = 0, "tac" = 0, 
                "tat" = 0, "tgc" = 0, "tgt" = 0, "gaa" = 0, "gag" = 0, "caa" = 0, "cag" = 0, "aaa" = 0, "aag" = 0, "tgg" = 0, 
                "atg" = 0, "tag" = 0, "taa" = 0, "tga" = 0)
  
  A_lyst = list("tct" = 0, "tca" = 0, "tcg" = 0, "agc" = 0, "agt" = 0, "tcc" = 0, "cgt" = 0, "cga" = 0, "cgg" = 0, "aga" = 0, 
                "agg" = 0, "cgc" = 0, "ctt" = 0, "cta" = 0, "ctg" = 0, "tta" = 0, "ttg" = 0, "ctc" = 0, "gct" = 0, "gca" = 0, 
                "gcg" = 0, "gcc" = 0, "cct" = 0, "cca" = 0, "ccg" = 0,"ccc" = 0, "gtt" = 0, "gta" = 0, "gtg" = 0,"gtc" = 0, "act" = 0, 
                "aca" = 0, "acg" = 0, "acc" = 0, "ggc" = 0, "gga" = 0, "ggg" = 0, "ggt" = 0, "att" = 0, "atc" = 0, "ata" = 0, 
                "ttc" = 0, "ttt" = 0, "aac" = 0, "aat" = 0, "cac" = 0, "cat" = 0, "gac" = 0, "gat" = 0, "tac" = 0, "tat" = 0, 
                "tgc" = 0, "tgt" = 0, "gaa" = 0, "gag" = 0, "caa" = 0, "cag" = 0, "aaa" = 0, "aag" = 0, "tgg" = 0, "atg" = 0, "tag" = 0, "taa" = 0, "tga" = 0)
  
  plus1_lyst =  list("tct" = 0, "tca" = 0, "tcg" = 0, "agc" = 0, "agt" = 0, "tcc" = 0, "cgt" = 0, "cga" = 0, "cgg" = 0, "aga" = 0, 
                     "agg" = 0, "cgc" = 0, "ctt" = 0, "cta" = 0, "ctg" = 0, "tta" = 0, "ttg" = 0, "ctc" = 0, "gct" = 0, "gca" = 0, 
                     "gcg" = 0, "gcc" = 0, "cct" = 0, "cca" = 0, "ccg" = 0,"ccc" = 0, "gtt" = 0, "gta" = 0, "gtg" = 0,"gtc" = 0, 
                     "act" = 0, "aca" = 0, "acg" = 0, "acc" = 0, "ggc" = 0, "gga" = 0, "ggg" = 0, "ggt" = 0, "att" = 0, "atc" = 0, 
                     "ata" = 0, "ttc" = 0, "ttt" = 0, "aac" = 0, "aat" = 0, "cac" = 0, "cat" = 0, "gac" = 0, "gat" = 0, "tac" = 0, 
                     "tat" = 0, "tgc" = 0, "tgt" = 0, "gaa" = 0, "gag" = 0, "caa" = 0, "cag" = 0, "aaa" = 0, "aag" = 0, "tgg" = 0, 
                     "atg" = 0, "tag" = 0, "taa" = 0, "tga" = 0)
  
  plus2_lyst = list("tct" = 0, "tca" = 0, "tcg" = 0, "agc" = 0, "agt" = 0, "tcc" = 0, "cgt" = 0, "cga" = 0, "cgg" = 0, "aga" = 0, 
                    "agg" = 0, "cgc" = 0, "ctt" = 0, "cta" = 0, "ctg" = 0, "tta" = 0, "ttg" = 0, "ctc" = 0, "gct" = 0, "gca" = 0, 
                    "gcg" = 0, "gcc" = 0, "cct" = 0, "cca" = 0, "ccg" = 0,"ccc" = 0, "gtt" = 0, "gta" = 0, "gtg" = 0,"gtc" = 0, 
                    "act" = 0, "aca" = 0, "acg" = 0, "acc" = 0, "ggc" = 0, "gga" = 0, "ggg" = 0, "ggt" = 0, "att" = 0, "atc" = 0, 
                    "ata" = 0, "ttc" = 0, "ttt" = 0, "aac" = 0, "aat" = 0, "cac" = 0, "cat" = 0, "gac" = 0, "gat" = 0, "tac" = 0, 
                    "tat" = 0, "tgc" = 0, "tgt" = 0, "gaa" = 0, "gag" = 0, "caa" = 0, "cag" = 0, "aaa" = 0, "aag" = 0, "tgg" = 0, 
                    "atg" = 0, "tag" = 0, "taa" = 0, "tga" = 0)
  totals_lyst = list("tct" = 0, "tca" = 0, "tcg" = 0, "agc" = 0, "agt" = 0, "tcc" = 0, "cgt" = 0, "cga" = 0, "cgg" = 0, "aga" = 0, 
                     "agg" = 0, "cgc" = 0, "ctt" = 0, "cta" = 0, "ctg" = 0, "tta" = 0, "ttg" = 0, "ctc" = 0, "gct" = 0, "gca" = 0, 
                     "gcg" = 0, "gcc" = 0, "cct" = 0, "cca" = 0, "ccg" = 0,"ccc" = 0, "gtt" = 0, "gta" = 0, "gtg" = 0,"gtc" = 0, 
                     "act" = 0, "aca" = 0, "acg" = 0, "acc" = 0, "ggc" = 0, "gga" = 0, "ggg" = 0, "ggt" = 0, "att" = 0, "atc" = 0, 
                     "ata" = 0, "ttc" = 0, "ttt" = 0, "aac" = 0, "aat" = 0, "cac" = 0, "cat" = 0, "gac" = 0, "gat" = 0, "tac" = 0, 
                     "tat" = 0, "tgc" = 0, "tgt" = 0, "gaa" = 0, "gag" = 0, "caa" = 0, "cag" = 0, "aaa" = 0, "aag" = 0, "tgg" = 0, 
                     "atg" = 0, "tag" = 0, "taa" = 0, "tga" = 0)
  
  
  
  for(x in series){
    
    
    count <- spliced_counts[[{{ID}}]][x]
    sequence <- spliced_sequences[[{{ID}}]][(x - 9): (x + 11)]
    
    codons <- 
      unlist(lapply(list(sequence[1:3], sequence[4:6], sequence[7:9], sequence[10:12], sequence[13:15], 
                         sequence[16:18], sequence[19:21]), paste0, collapse = ""))
    
    
    
    minus4_lyst[[codons[1]]] <- minus4_lyst[[codons[1]]] + count
    minus3_lyst[[codons[2]]] <- minus3_lyst[[codons[2]]] + count
    E_lyst[[codons[3]]] <- E_lyst[[codons[3]]] + count
    P_lyst[[codons[4]]] <- P_lyst[[codons[4]]] + count
    A_lyst[[codons[5]]] <- A_lyst[[codons[5]]] + count
    plus1_lyst[[codons[6]]] <- plus1_lyst[[codons[6]]] + count
    plus2_lyst[[codons[7]]] <- plus2_lyst[[codons[7]]] + count
    
    totals_lyst[[codons[1]]] <- totals_lyst[[codons[1]]] + count
    totals_lyst[[codons[2]]] <- totals_lyst[[codons[2]]] + count
    totals_lyst[[codons[3]]] <- totals_lyst[[codons[3]]] + count
    totals_lyst[[codons[4]]] <- totals_lyst[[codons[4]]] + count
    totals_lyst[[codons[5]]] <- totals_lyst[[codons[5]]] + count
    totals_lyst[[codons[6]]] <- totals_lyst[[codons[6]]] + count
    totals_lyst[[codons[7]]] <- totals_lyst[[codons[7]]] + count
    
    
    
  }
  
  df = data.frame(transcript = {{ID}},
                  codon = names(minus4_lyst),
                  minus4 = unlist(minus4_lyst, use.names = F),
                  minus3 = unlist(minus3_lyst, use.names = F),
                  E_site = unlist(E_lyst, use.names = F),
                  P_site = unlist(P_lyst, use.names = F),
                  A_site = unlist(A_lyst, use.names = F),
                  plus1 = unlist(plus1_lyst, use.names = F),
                  plus2 = unlist(plus2_lyst, use.names = F),
                  total = unlist(totals_lyst, use.names = F))
  
  return(df)
  
}

# define home
home = "\\\\data.beatson.gla.ac.uk/data"
#home = "/home/local/BICR/jettles/data"

# define parent directory
parent_dir = file.path(home, "R11/James/CNOT3_Sel_RiboSeq/FINAL")

# read in codons vector

codons_vector <- 
  read_csv(file.path(home, "R11/bioinformatics_resources/useful_tables/codon_box_types.csv")) %>% 
  pull(codon) %>% 
  tolower() %>% 
  str_replace_all("u", "t")

# add stop codons
codons_vector <- c(codons_vector, "taa", "tag", "tga")

# check to see if all 64 codons present
length(codons_vector)

# read in master and define region lengths
master <- as_tibble(fread(file.path(home, "R11/James/sequences/master.csv"), drop = "V1", header = T))

master$fputr_length = nchar(master$nucleotide_sequence_fpUTR)
master$cds_length = nchar(master$nucleotide_sequence_CDS)
master$tputr_length = nchar(master$nucleotide_sequence_tpUTR)
master$sequence = paste0(master$nucleotide_sequence_fpUTR,
                         master$nucleotide_sequence_CDS,
                         master$nucleotide_sequence_tpUTR)
master$full_length = master$fputr_length + master$cds_length + master$tputr_length

master <- 
  master %>% select(ENST, fputr_length, cds_length, tputr_length, full_length, sequence)

# read in counts file

#counts_dir = file.path(parent_dir, "Counts_files")
counts_dir = file.path(parent_dir, "ribosome_occupancy")

# here want to create a filter so that transcripts are at least 21 nts long
# 21 nts = 7 codons (-4, -3, E, P, A, +1, +2)
counts <- 
  read_tsv(file.path(counts_dir, "CNOT3_SelRibo_counts.tsv"))%>% 
  inner_join(master, by = "ENST") %>% 
  filter(cds_length >= 300)

# define regions table
CNOT3_regulated_transcripts <- 
  read_csv(file.path(counts_dir, "CNOT3_monosomes_transcripts.csv")) %>% 
  pull(x)

signal_peptides <- 
  fread(file.path(home, "R11/James/sequences/features.csv")) %>% 
  select(ENST, cds_length, signal_sequence) %>% 
  filter(cds_length > 300) %>% 
  filter(ENST %in% CNOT3_regulated_transcripts) %>% 
  select(-cds_length)

counts <- 
  counts %>% 
  filter(ENST %in% signal_peptides$ENST) %>% 
  inner_join(signal_peptides, by = "ENST")

counts$SP_start <- 10
counts$SP_end <- 90

# split counts
counts_split <- split(counts, f = counts$sample)


data_lyst <- list()



for(le_sample in counts_split){
  
  
  print(unique(le_sample$sample))
  
  # extract a sample for development
  counts_table <- le_sample 
  
  print(counts_table)
  
  counts_table$ID <- counts_table$ENST
  
  # create vector of transcripts to analyse
  transcripts <- counts_table$ID
  
  #print(transcripts)
  
  # create a list of counts
  counts_lyst <- extract_counts(counts_table)
  
  #print(counts_lyst[1:100])
  
  # create a list of sequences
  sequences_lyst <- extract_sequences(counts_table)
  
  # splice CDS section for both counts and sequences list
  spliced_counts <- lapply(1:length(transcripts), splice_cds, counts_lyst)
  spliced_counts <- lapply(spliced_counts, function(x){x[1:90]})
  
  spliced_sequences <- lapply(1:length(transcripts), splice_cds, sequences_lyst)
  spliced_sequences <- lapply(spliced_sequences, function(x){x[1:90]})
  
  names(spliced_counts) <- transcripts
  names(spliced_sequences) <- transcripts
  
  
  # Assign counts to each site in the ribosome for each sliding window and sum across all SLs
  one_sample <- lapply(1:length(transcripts), count_codon_occupancy)
  
  final_output <- do.call("rbind", one_sample)
  
  final_output <- 
    final_output %>% 
    as_tibble() %>% 
    mutate(sample_id = unique(counts_table$sample))
  
  
  data_lyst[[unique(le_sample$sample)]] <- final_output
  
}


final_data <- do.call("rbind", data_lyst)

write.csv(final_data, "C:/Users/Jettles/OneDrive - University of Glasgow/Documents/get_shit_done/corrections/CNOT3SelRP_RDTs_1_30.csv")

final_data <- read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/get_shit_done/corrections/CNOT3SelRP_RDTs_1_30.csv") %>% select(-`...1`)

SPs <- read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/features.csv") %>% select(ENST, signal_sequence)


dts_summarised <- 
  final_data  %>% 
  group_by(sample_id, signal_sequence, codon) %>% 
  summarise(`-4` = sum(minus4),
            `-3` = sum(minus3),
            `E` = sum(E_site),
            `P` = sum(P_site),
            `A` = sum(A_site),
            `+1` = sum(plus1),
            `+2` = sum(plus2),
            total = sum(total))


dts_summarised_split <- split(dts_summarised, dts_summarised$sample_id)

downstream_processing <- function(df){
  
  df <-
    df %>%
    filter(codon != "tga" & codon != "tag" & codon != "taa") %>%
    mutate(minus_4 = (log2(`-4`)) - log2(total / 7),
           minus_3 = (log2(`-3`)) - log2(total / 7),
           E_site = (log2(E)) - log2(total / 7),
           P_site = (log2(P)) - log2(total / 7),
           A_site = (log2(A)) - log2(total / 7),
           plus_1 = (log2(`+1`)) - log2(total / 7),
           plus_2 = (log2(`+2`)) - log2(total / 7)) %>%
    select(codon, minus_4, minus_3, E_site, P_site, A_site, plus_1, plus_2, sample_id) %>%
    dplyr::rename(sample = sample_id)
  
  return(df)
}

data <- 
  do.call("rbind", lapply(dts_summarised_split, downstream_processing)) %>% 
  ungroup()

data$replicate <- unlist(lapply(str_split(data$sample, "_"), function(x){x[[1]]}))
data$condition <- unlist(lapply(str_split(data$sample, "_"), function(x){x[[3]]}))

ctrl_data <- 
  data %>% 
  select(-sample) %>% 
  gather(key = "ribosome_position", value = "log2(count_actual/count_expected)", minus_4:plus_2) 

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
            axis.text.x = element_text(size = 16, color = "black"),
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
            plot.margin=unit(c(2,2,2,10),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

ctrl_RDTs <- 
  ctrl_data %>% 
  filter(condition == "Tot" & signal_sequence == T) %>% 
  group_by(codon, ribosome_position) %>% 
  summarise(mean_RDT = mean(`log2(count_actual/count_expected)`))%>% 
  filter(ribosome_position == "A_site") %>% 
  arrange(-mean_RDT)

"#377EB8"
"#984EA3"
"#4DAF4A"
"#BF5B17"
"#FF7F00"
"#F0027F"
"#999999"

my_colour_vector <- 
  case_when(ctrl_RDTs$codon == "ctg" ~ "#377EB8",
            ctrl_RDTs$codon == "ctc" ~ "#984EA3", 
            ctrl_RDTs$codon == "ctt" ~ "#4DAF4A", 
            ctrl_RDTs$codon == "cta" ~ "#FF7F00", 
            ctrl_RDTs$codon == "ttg" ~ "#BF5B17", 
            ctrl_RDTs$codon == "tta" ~ "#F0027F", 
            .default = "#999999")

my_colour_vector2 <- c("#377EB8",
                       "#984EA3",
                       "#4DAF4A",
                       "#BF5B17",
                       "#FF7F00",
                       "#F0027F",
                       "#999999")

#### here ####
ctrl_RDTs_gplot <- 
  ctrl_RDTs %>%   
  mutate(leu_codons = factor(
    case_when(codon == "ctg" ~ "ctg",
              codon == "ctc" ~ "ctc", 
              codon == "ctt" ~ "ctt", 
              codon == "cta" ~ "cta", 
              codon == "ttg" ~ "ttg", 
              codon == "tta" ~ "tta",
              .default = "other"),
    levels = c("ctg", "ctc", "ctt", "cta", "ttg", "tta", "other"))) %>% 
  ggplot(aes(x = reorder(codon, -mean_RDT), y = mean_RDT, fill = leu_codons)) +
  geom_col() +
  scale_fill_manual(values = my_colour_vector2) +
  ggtitle("A Site Ribosome Dwell Time") +
  publication_theme() +
  ylab("Log2FC (Observed/Expected Counts)") +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.25,
                                   color = my_colour_vector, 
                                   face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "none")




deltas_final <- 
  data %>% 
  select(-sample) %>% 
  gather(key = "ribosome_position", value = "log2(count_actual/count_expected)", minus_4:plus_2) %>% 
  spread(key = condition, value = "log2(count_actual/count_expected)") %>% 
  mutate(delta = CNOT3 - Tot) %>% 
  select(signal_sequence, codon, replicate, ribosome_position, delta) 


deltas_SPs_ctrl <- 
  deltas_final %>% 
  filter(signal_sequence == T & ribosome_position == "A_site") %>% 
  select(-ribosome_position, -signal_sequence) %>% 
  group_by(codon) %>% 
  summarise(mean_delta = mean(delta)) %>% 
  arrange(-mean_delta)

my_colour_vector3 <- 
  case_when(deltas_SPs_ctrl$codon == "ctg" ~ "#377EB8",
            deltas_SPs_ctrl$codon == "ctc" ~ "#984EA3", 
            deltas_SPs_ctrl$codon == "ctt" ~ "#4DAF4A", 
            deltas_SPs_ctrl$codon == "cta" ~ "#FF7F00", 
            deltas_SPs_ctrl$codon == "ttg" ~ "#BF5B17", 
            deltas_SPs_ctrl$codon == "tta" ~ "#F0027F", 
            .default = "#999999")

### and here ###
deltas_plot <- 
  deltas_SPs_ctrl %>%   
  mutate(leu_codons = factor(
    case_when(codon == "ctg" ~ "ctg",
              codon == "ctc" ~ "ctc", 
              codon == "ctt" ~ "ctt", 
              codon == "cta" ~ "cta", 
              codon == "ttg" ~ "ttg", 
              codon == "tta" ~ "tta",
              .default = "other"),
    levels = c("ctg", "ctc", "ctt", "cta", "ttg", "tta", "other"))) %>% 
  ggplot(aes(x = reorder(codon, -mean_delta), y = mean_delta, fill = leu_codons)) +
  geom_col() +
  scale_fill_manual(values = my_colour_vector2) +
  ggtitle("CNOT3 Enrichment in A site RDT") +
  publication_theme() +
  ylab("(CNOT3 IP'd - Total RDTs)")+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.25,
                                   color = my_colour_vector3, 
                                   face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "none")

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/new_RDT_figure.png", res = 300, height = 3000, width = 3000)
cowplot::plot_grid(ctrl_RDTs_gplot,deltas_plot, nrow = 2, ncol = 1, align = "v")
dev.off()