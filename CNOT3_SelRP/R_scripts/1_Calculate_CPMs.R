#load packages----
library(tidyverse)
library(data.table)
library(purrr)

#home <- "/home/local/BICR/jettles/data"
home <- "\\\\data.beatson.gla.ac.uk/data"

setwd(file.path(home, "R11/James/CNOT3_Sel_RiboSeq/FINAL/R_scripts"))

parent_dir2 <- file.path(home, "R11/James/CNOT3_Sel_RiboSeq/FINAL")

#read in common variables----
source("common_variables.R")

#set the threshold for the average CDS counts a transcript has to have across all samples for it to be included
min_counts <- 50

#set the thresholds for region lengths
UTR5_min_len <- 25
CDS_min_len <- 300
UTR3_min_len <- 0

region_cutoffs <- c(UTR5_min_len, CDS_min_len, UTR3_min_len)

#read in data----
region_lengths <- read_csv(file = file.path(home, "R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_region_lengths.csv"), col_names = c("transcript", "UTR5_len", "CDS_len", "UTR3_len"))

#read in the CDS counts to use to filter the data
#the following for loop reads in each final CDS counts file and renames the counts column by the sample name and saves each data frame to a list
data_list <- list()

RPF_sample_names

  
  for(sample_id in RPF_sample_names){
    
    # generate path top table
    path_to_counts <- file.path(parent_dir2,
                                "Analysis/CDS_counts",
                                paste0(sample_id, "_pc_final_counts_all_frames.csv"))
    
    # read in data
    CDS_counts <- read_csv(path_to_counts) 
    
    #rename column names
    colnames(CDS_counts) <- c("transcript", paste0(sample_id))
    
    # place table in list
    data_list[[sample_id]] <- CDS_counts
    
  }
  


#extract transcript IDs to keep----
#merge all data within the above list using reduce
#remove any transcripts with an NA in any sample (these will have had 0 counts)
#remove counts without at least 10 counts in each sample
#remove any transcripts with a mean count (across all samples) less than the min_counts threshold
#remove any transcripts with UTRs/CDSs shorter than thresholds defined above
#extract the transcript IDs
data_list %>%
  purrr::reduce(full_join, by = "transcript") %>%
  column_to_rownames("transcript") %>%
  drop_na() %>%
  filter_all(all_vars(. > 10)) %>% 
  filter(rowMeans(.) >= min_counts) %>%
  rownames_to_column("transcript") %>%
  inner_join(region_lengths, by = "transcript") %>%
  filter(UTR5_len >= UTR5_min_len & CDS_len >= CDS_min_len & UTR3_len >= UTR3_min_len) %>% # Removing smaller transcripts for binning functions later on. If a transcript has a 5'UTR of only 10 nucleotides, it doesn't really make sense to bin this!
  pull(transcript) -> filtered_transcripts


print("reading in positional counts data")

#read in counts data----
counts_list <- list()
UTR5_list <- list()

extract_CDS <- function(data){
  
  fputr_len <- unique(data$UTR5_len)
  cds_len <- unique(data$CDS_len)
  
  a <- data %>% 
    filter(Position > fputr_len) %>% 
    slice_head(n = cds_len) %>%
    select(transcript, Position, Nucleotide, Counts, CDS_len, CPM) 
  
  return(a)
}

add_codon_vector <- function(data) {
  
  data <- 
    data %>%
    mutate(codon = rep(1:(CDS_len[1] / 3), each = 3))
  
  return(data)
  
}

extract_UTR5 <- function(data){
  
  fputr_len <- unique(data$UTR5_len)
  
  a <- data %>% 
    filter(Position <= fputr_len) %>% 
    select(transcript, Position, Nucleotide, Counts, CPM) 
  
  return(a)
}



sample_id = "REP2_M_Tot"

  for(sample_id in RPF_sample_names){
    
    #replicate
    REP <- str_extract(sample_id, "REP\\d")
    
    #condition1
    condition1 <- str_remove_all(str_extract(sample_id, "_.*_"),"_")
    
    #condition2
    condition2 <- str_remove(str_extract(sample_id,"_.*$"),"_.*_")
    
    # generate path to csv
    path2csv <-  file.path(parent_dir2,
                           "Counts_files/csv_files",
                           paste0(sample_id, "_pc_final_counts.csv"))
    
    # print message to the command line
    print(paste("reading in file", paste0(sample_id)))
    
    # read in data using fread
    dat_temp <- fread(path2csv, header = T)
    
    # Make this table into a tibble
    dat_temp <- as_tibble(dat_temp)
    
    # filter transcripts 
    dat_temp <-
      dat_temp %>%
      filter(transcript %in% filtered_transcripts)
    
    # join with region lengths
    dat_temp <- inner_join(dat_temp,
                           region_lengths[c("transcript", "UTR5_len", "CDS_len")],
                           by = "transcript")
    
    # Calculate total counts per sample
    summed_counts <- sum(dat_temp$Counts)
    
    # generate counts per million by dividing each individual count per nucleotide by total counts of that sample * 1,000,000
    dat_temp_cpm <- dat_temp %>% mutate(CPM = Counts/summed_counts * 1000000)
    
    # Split this big data table into a list where each element in the list is a transcript
    dat_temp_split <- split(dat_temp_cpm, f = dat_temp_cpm$transcript)
    
    print("Extracting CDS")
    
    extracted_CDSs <- lapply(dat_temp_split, extract_CDS)
    
    print("Creating codon vector")
    
    added_codon_vectors <- lapply(extracted_CDSs, add_codon_vector)
    
    print("Summarise nucleotide counts to codon counts")
    
    big_table <- do.call("rbind", added_codon_vectors)
    
    # sum counts per codon and CPMs per codon
    codon_counts <-
      big_table %>% 
      group_by(transcript, codon) %>%
      summarise(summed_counts_codon = sum(Counts),
                summed_CPM_codon = sum(CPM)) %>%
      mutate(IP = rep(condition2),
             MD30 = rep(condition1),
             replicate = REP) %>% 
      ungroup()
    
    # get total counts per transcript
    codon_counts %>% 
      group_by(replicate, MD30, IP, transcript) %>% 
      summarise(total_counts_per_transcript = sum(summed_counts_codon)) %>% 
      ungroup() -> counts_per_transcript
    
    # inner join
    int <- 
      codon_counts %>% 
      inner_join(counts_per_transcript,
                 by = c("replicate", "MD30", "IP", "transcript"))
    
    #internally normalise the counts
    int_normalised<- 
      int %>% mutate(perc_counts_across_CDS = 
                       int$summed_counts_codon/int$total_counts_per_transcript * 100)
    
    # tidy up data and place in list
    int_normalised %>% 
      select(replicate, 
             MD30,
             IP,
             transcript, 
             codon, 
             summed_counts_codon, 
             summed_CPM_codon,
             perc_counts_across_CDS) -> counts_list[[paste0(sample_id)]]
    
    
    
    # extract 5'UTR
    extracted_UTR5s <- lapply(dat_temp_split, extract_UTR5)
    
    UTR5_table <- 
      do.call("rbind", extracted_UTR5s) %>% 
      mutate(MD30 = condition1,
             IP = condition2,
             replicate = REP)
    
    UTR5_list[[sample_id]] <- UTR5_table
    
  }


all_data <- do.call("rbind", counts_list)

print("writing CDS file")

fwrite(all_data, file = file.path(home, "R11/James/CNOT3_Sel_RiboSeq/FINAL/CPMs/CDS_CPMs.csv"))

UTR5_big_table <- do.call("rbind", UTR5_list)

print("writing 5'UTR file")

fwrite(UTR5_big_table, file = file.path(home, "R11/James/CNOT3_Sel_RiboSeq/FINAL/CPMs/UTR5_CPMs.csv"))