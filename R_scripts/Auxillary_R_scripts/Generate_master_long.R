library(tidyverse)

#paths
home <- "\\\\data.beatson.gla.ac.uk/data"

master <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/master.csv") %>% 
  select(ENST, nucleotide_sequence_CDS)

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

tidy_data <- function(x){
  
  transcript <- master[x,1]
  nucleotide_sequence <- split_into_codons(master[x,2])
  
  return(data.frame(transcript = transcript, 
                    codon_position = 1:length(nucleotide_sequence),
                    sequence = nucleotide_sequence))
}

tidied_data <- do.call("rbind", lapply(1:nrow(master), tidy_data))

properties <- 
read_csv(file.path(home,"R11/bioinformatics_resources/useful_tables/codon_amino_acid_pairs_and_properties.csv")) %>% 
  #select(-Symbol) %>% 
  mutate(codon =  tolower(str_replace_all(codon, "U", "T")))


properties <- 
rbind(properties,
tibble(codon = c("taa", "tag", "tga"),
       AA = c("Stp", "Stp", "Stp"),
       Symbol = c("*", "*", "*"),
       Properties = rep("Stp", 3)))
  

all_data <- inner_join(tidied_data, properties, by = c("sequence" = "codon"))

write.csv(all_data, file = file.path(home, "R11/James/sequences/master_long.csv"))