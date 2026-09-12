library(tidyverse)
library(data.table)

TM_data <- 
read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/get_shit_done/TM_intermediate_table.csv") %>% 
  select(ENST:codon_end)

TM_data <- 
  TM_data %>% 
  mutate(tag = paste0(ENST, "_", class, "_", region, "_", location, "_", ID, "_", codon_start, "_", codon_end))

TM_data_split <- split(TM_data, TM_data$tag)

TM_data_split[1:10]

append_raw_codon_positions <- 
  function(x){
    
    start_codon = unique(x$codon_start)
    end_codon = (start_codon + nrow(x) -1)
    
    y <- 
      x %>% 
      mutate(raw_codon = start_codon:end_codon) %>% 
      select(-codon_start, -codon_end, -tag)
    
    return(y)
    
  }


raw_codons <- 
do.call("rbind", lapply(TM_data_split, append_raw_codon_positions))

raw_codons_split_arranged <- lapply(split(raw_codons, raw_codons$ENST), arrange, raw_codon)

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

master <- read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/master.csv")

# sanity checks
master %>% 
  select(ENST, nucleotide_sequence_CDS) %>% 
  filter(ENST == "ENST00000005284.4") %>% 
  pull(nucleotide_sequence_CDS) %>% 
  split_into_codons()

raw_codons$class <- factor(raw_codons$class) 
summary(raw_codons$class)  

final <- 
raw_codons %>% 
  rename(rel_position = position,
         codon_position = raw_codon)

write.csv(final, "C:/Users/Jettles/OneDrive - University of Glasgow/Documents/TM_domains_FINAL.csv")
