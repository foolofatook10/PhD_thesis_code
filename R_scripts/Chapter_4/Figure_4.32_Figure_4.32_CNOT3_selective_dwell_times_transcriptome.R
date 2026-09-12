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
  
  print(paste("% complete:", round(int/length(transcripts) * 100, digits = 2)))
  
  ID = transcripts[[int]]
  #print(ID)
  
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
parent_dir = file.path(home, "R11/James/CNOT3_RiboSeq")

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



# read in merged data
data <- read_csv(file.path(parent_dir, "Analysis/DESeq2_output/merged_DESeq2.csv"))

master <- master %>% filter(ENST %in% unique(data$transcript))


# read in counts file

#counts_dir = file.path(parent_dir, "Counts_files")
#counts_dir = "C:/Users/Jettles/OneDrive - University of Glasgow/Documents"

# here want to create a filter so that transcripts are at least 21 nts long
# 21 nts = 7 codons (-4, -3, E, P, A, +1, +2)
counts <- 
  read_tsv(file.path(home, "R11/James/CNOT3_Riboseq/Counts_files/counts_combined.tsv"))%>% 
  inner_join(master, by = "ENST") %>% 
  filter(cds_length >= 21)


# split counts
counts_split <- split(counts, f = counts$sample)

data_lyst <- list()


for(le_sample in counts_split){
  
  print(unique(le_sample$sample))
  
  # extract a sample for development
  counts_table <- le_sample 
  
  
  
  counts_table$ID <- counts_table$ENST
  
  
  # create vector of transcripts to analyse
  transcripts <- counts_table$ID
  
  # create a list of counts
  counts_lyst <- extract_counts(counts_table)
  
  #print(counts_lyst[1:100])
  
  # create a list of sequences
  sequences_lyst <- extract_sequences(counts_table)
  
  #length(sequences_lyst[[2]]) == length(counts_lyst[[2]])
  
  # splice CDS section for both counts and sequences list
  spliced_counts <- lapply(1:length(transcripts), splice_cds, counts_lyst)
  spliced_sequences <- lapply(1:length(transcripts), splice_cds, sequences_lyst)
  
  #length(spliced_sequences[[2]]) == length(spliced_counts[[2]])
  
  names(spliced_counts) <- transcripts
  names(spliced_sequences) <- transcripts
  
  
  # Assign counts to each site in the ribosome for each sliding window and sum across all SLs
  one_sample <- lapply(1:length(transcripts), count_codon_occupancy)
  
  final_output <- 
    do.call("rbind", one_sample) %>% 
    mutate(sample_id = unique(counts_table$sample))
  
  data_lyst[[unique(counts_table$sample)]] <- final_output
  
}

final <- do.call("rbind", data_lyst)
rownames(final) <- NULL

write.csv(final, file = file.path(home, "JETTLES/thesis_figures/dwell_times_CNOT3SelRP_entire_transcriptome.csv"))
# 
#   final_output <- 
#     final_output %>% 
#     as_tibble() %>% 
#     group_by(codon) %>% 
#     summarise(`-4` = sum(minus4),
#               `-3` = sum(minus3),
#               `E` = sum(E_site),
#               `P` = sum(P_site),
#               `A` = sum(A_site),
#               `+1` = sum(plus1),
#               `+2` = sum(plus2),
#               total = sum(total)) %>% 
#     mutate(sample_id = unique(counts_table$sample))
#   
#   
#   data_lyst[[unique(le_sample$sample)]] <- final_output
#   print(final_output)
#   
# }
# 
# 
# 
# downstream_processing <- function(df){
#   
#   df <- 
#   df %>%
#     filter(codon != "tga" & codon != "tag" & codon != "taa") %>%
#     mutate(minus_4 = (log2(`-4`)) - log2(total / 7),
#            minus_3 = (log2(`-3`)) - log2(total / 7),
#            E_site = (log2(E)) - log2(total / 7),
#            P_site = (log2(P)) - log2(total / 7),
#            A_site = (log2(A)) - log2(total / 7),
#            plus_1 = (log2(`+1`)) - log2(total / 7),
#            plus_2 = (log2(`+2`)) - log2(total / 7)) %>% 
#     select(codon, minus_4, minus_3, E_site, P_site, A_site, plus_1, plus_2, sample_id) %>% 
#     dplyr::rename(sample = sample_id)
#   
#   return(df)
# }
# 
# 
# 
# data <- 
# do.call("rbind", lapply(data_lyst, downstream_processing))
# 
# #gather the data into tidy format----
# data %>%
#   pull(codon) %>%
#   unique() -> codons
# 
# RPF_sample_names <- unique(data$sample)
# 
# gathered_list <- list()
# for (sample in RPF_sample_names) {
#   for (codon in codons) {
#     data[data$codon == codon & data$sample == sample,] %>%
#       select(-sample) %>%
#       gather(key = codon, value = freq) %>%
#       rename(position = codon) %>%
#       mutate(codon = rep(codon),
#              sample = factor(rep(sample))) -> gathered_list[[paste(codon, sample, sep = "_")]]
#   }
# }
# gathered_data <- do.call("rbind", gathered_list)
# 
# gathered_data %>%
#   mutate(position = as.numeric(case_when(position == "minus_4" ~ -4,
#                                          position == "minus_3" ~ -3,
#                                          position == "E_site" ~ -2,
#                                          position == "P_site" ~ -1,
#                                          position == "A_site" ~ 0,
#                                          position == "plus_1" ~ 1,
#                                          position == "plus_2" ~ 2)),
#          wobble = factor(str_sub(codon, 3,3))) -> plot_data
# 
# for (sample in RPF_sample_names) {
#   plot_title <- str_replace_all(sample, "_", " ")
#   plot_title <- str_remove(plot_title, " RPFs")
#   
#   plot_data[plot_data$sample == sample,] %>%
#     ggplot(aes(x = position, y = freq, colour = wobble))+
#     geom_point(size = 2)+
#     stat_summary(fun=mean, geom="line", size = 1)+
#     ylab("normalised codon frequency")+
#     scale_x_continuous(limits = c(-4,2), breaks = -4:2)+
#     xlab("Codon position (0 = A-site)")+
#     ggtitle(plot_title) -> codon_plot
#   
#   print(codon_plot)
#   
# }
# 
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "-3") %>% arrange(-freq)})
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "-2") %>% arrange(-freq)})
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "-1") %>% arrange(-freq)})
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "0") %>% arrange(-freq)})
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "1") %>% arrange(-freq)})
# 
# 
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "-1") %>% arrange(-freq)})
# 
#        
# mean_data <- 
# plot_data %>% 
#   mutate(condition = unlist(lapply(str_split(sample, "_"),function(x){x[[3]]}))) %>% 
#   group_by(condition, position, codon, wobble) %>% 
#   summarise(mean_freq = mean(freq)) 
# 
# cca <- mean_data %>% filter(codon == "cca" & position == -2)
# cct <- mean_data %>% filter(codon == "cct"& position == -2)
# gaa <- mean_data %>% filter(codon == "gaa"& position == 0)
# tac <- mean_data %>% filter(codon == "tac"& position == 0)
# cgg <- mean_data %>% filter(codon == "cgg"& position == -1)
# cga <- mean_data %>% filter(codon == "cga"& position == -1)
# agg <- mean_data %>% filter(codon == "agg"& position == -1)
# gat <- mean_data %>% filter(codon == "gat"& position == -3)
# cgt <- mean_data %>% filter(codon == "cgt"& position == 1)
# 
# # read in publication theme
# publication_theme <- function(base_size=14, base_family="helvetica") {
#   library(grid)
#   library(ggthemes)
#   (theme_foundation(base_size=base_size, base_family=base_family)
#     + theme(plot.title = element_text(face = "bold",
#                                       size = rel(1.2), hjust = 0.5),
#             text = element_text(),
#             panel.background = element_rect(colour = NA),
#             plot.background = element_rect(colour = NA),
#             panel.border = element_rect(colour = NA),
#             axis.title = element_text(face = "bold",size = rel(1)),
#             axis.title.y = element_text(angle=90,vjust =2),
#             axis.title.x = element_text(vjust = -0.2),
#             axis.text.x = element_text(size = 16, color = "black"),
#             axis.text.y = element_text(size = 14, color = "black"),
#             axis.text = element_text(), 
#             axis.line = element_line(colour="black"),
#             axis.ticks = element_line(),
#             panel.grid.major = element_line(colour="#f0f0f0"),
#             panel.grid.minor = element_blank(),
#             legend.text = element_text(size = 16),
#             legend.key = element_rect(colour = NA),
#             legend.position = "bottom",
#             legend.direction = "horizontal",
#             legend.key.size= unit(0.2, "cm"),
#             legend.margin = unit(0, "cm"),
#             legend.title = element_text(face="italic", size = 16),
#             plot.margin=unit(c(2,2,2,2),"mm"),
#             strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
#             strip.text = element_text(face="bold")
#     ))
#   
# }
# 
# mean_data %>%   
# ggplot(aes(x = position, y = mean_freq, colour = wobble)) +
#   geom_point(size = 2) +
#   geom_text(data = cca, aes(x = position, y = mean_freq, label = codon)) +
#   geom_text(data = cct, aes(x = position, y = mean_freq, label = codon)) +
#   geom_text(data = gaa, aes(x = position, y = mean_freq, label = codon)) +
#   geom_text(data = tac, aes(x = position, y = mean_freq, label = codon)) +
#   geom_text(data = cgg, aes(x = position, y = mean_freq, label = codon)) +
#   geom_text(data = cga, aes(x = position, y = mean_freq, label = codon)) +
#   geom_text(data = agg, aes(x = position, y = mean_freq, label = codon)) +
#   geom_text(data = gat, aes(x = position, y = mean_freq, label = codon)) +
#   geom_text(data = cgt, aes(x = position, y = mean_freq, label = codon)) +
#   facet_wrap(~condition, nrow = 2) +
#   stat_summary(fun=mean, geom="line", size = 1)+
#   scale_x_continuous(limits = c(-4,2), breaks = -4:2)+
#   xlab("Codon position (0 = A-site)") +
#   
#   publication_theme()
