library(data.table)
library(tidyverse)

# path to signal sequence information
transcripts_annotated_by_SigP <- "\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/canis_lupis.csv"
#transcripts_nucleotide_sequences <- "\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/mouse_master.csv"

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



# read in master table
master <- as_tibble(fread(file = transcripts_annotated_by_SigP, header = T, drop = "V1"))

# select the longest transcript per gene
longest <- master %>% 
  mutate(polypeptide_length = nchar(master$sequence)/3) %>%
  group_by(gene) %>%
  mutate(the_rank  = rank(-polypeptide_length, ties.method = "first")) %>%
  filter(the_rank == 1) %>% 
  select(-the_rank)



#longest <- longest %>% 
  #dplyr::rename(signal_peptide = signal_sequence) %>% 
  #dplyr::rename(sequence = nucleotide_sequence_CDS)

#Import codons table
path = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/useful_tables"
codons <- read_csv(paste0(path, "/codon_box_types.csv"))

#Make the capitalised codons lower case
codons$codon <- factor(tolower(codons$codon))

#replace letters "u" with "t"
codons$codon <- gsub("u","t", codons$codon)

# replace Ile 4-box with 3-box
codons <- codons %>% mutate(box = as.numeric(str_remove(case_when(AA == "Ile" ~ "3-box",
                                                                  .default = box), 
                                                        pattern = "-box")))

extract_pos <- function(x,y){
  x[[y]]
}

extract_codons_and_enumerate <- function(z, sequences) {
  
  df <- as.data.frame(table(unlist(lapply(sequences, extract_pos, y = z))))
  return(df %>% mutate(position = z))
}


process_data <- function(vector_of_sequences){
  
  sequences <- lapply(vector_of_sequences, split_into_codons)
  
  original_in <- length(sequences)
  
  print(paste("Total number of sequences in:", original_in))
  
  over_100_codons <- unlist(lapply(sequences, length)) > 100
  sequences <- sequences[over_100_codons]
  
  print(paste("Lost sequences (not at least 100 codons):", original_in - length(sequences)))
  
  
  near_final <- do.call("rbind", lapply(1:30, extract_codons_and_enumerate, sequences = sequences))
  colnames(near_final) <- c("codon", "codon_count", "codon_position")
  
  joined <- inner_join(near_final, codons, by = c("codon"))
  
  codon_frequency_table <- 
    inner_join(joined,
               joined %>% 
                 group_by(codon_position) %>% 
                 summarise(sum_count = sum(codon_count)),
               by = "codon_position") %>% 
    mutate(codon_percentage = codon_count/sum_count * 100)
  
  total_AAs_per_position <- 
    joined %>% 
    group_by(codon_position, AA) %>% 
    summarise(AA_count = sum(codon_count))
  
  final <- 
    inner_join(codon_frequency_table, 
               total_AAs_per_position, 
               by = c("codon_position", "AA")) %>% 
    mutate(RSCU = codon_count/AA_count * box)
  
  return(final)
  
}


number_of_signal_peptides <- 
  longest %>% 
  filter(signal_peptide == T) %>% 
  filter(nchar(sequence) > 300) %>% 
  nrow()

number_of_cytosolics <- 
  longest %>% 
  filter(signal_peptide == F) %>% 
  filter(nchar(sequence) > 300) %>% 
  nrow()

signal_sequence_information <- 
  longest %>% 
  filter(signal_peptide == T) %>% 
  pull(sequence) %>% 
  process_data() %>% 
  mutate(mRNAtype = paste0("Signal sequence\n = ", number_of_signal_peptides))

cytosolic_information <- 
  longest %>% 
  filter(signal_peptide == F) %>% 
  pull(sequence) %>% 
  process_data() %>% 
  mutate(mRNAtype = paste0("No signal sequence\n = ", number_of_cytosolics))

rbind(signal_sequence_information, cytosolic_information) %>% 
  filter(codon_position != 1) %>% 
  ggplot(aes(x = codon_position, y = RSCU, group = mRNAtype, colour = mRNAtype)) +
  geom_line() +
  facet_wrap(~codon)

final_information <- 
  rbind(signal_sequence_information, cytosolic_information) %>% 
  filter(codon_position != 1) 

final_information <- 
  final_information %>% 
  mutate(codons_to_colour = case_when(codon =="ctg" ~ "ctg",
                                      codon =="ctc" ~ "ctc",
                                      codon =="ctt" ~ "ctt",
                                      codon =="cta" ~ "cta",
                                      codon =="ttg" ~ "ttg",
                                      codon =="tta" ~ "tta",
                                      .default = "other"))

background <- final_information %>% filter(codons_to_colour == "other")
ctg <- final_information %>% filter(codons_to_colour == "ctg")
ctc <- final_information %>% filter(codons_to_colour == "ctc")
ctt <- final_information %>% filter(codons_to_colour == "ctt")
cta <- final_information %>% filter(codons_to_colour == "cta")
ttg <- final_information %>% filter(codons_to_colour == "ttg")
tta <- final_information %>% filter(codons_to_colour == "tta")

gplot <- 
  ggplot() +
  geom_line(data = background, 
            aes(x = codon_position, y = codon_percentage, group = codon),
            size = 1,
            colour = "#C8C8C8") +
  geom_line(data = ctg, 
            aes(x = codon_position, y = codon_percentage, group = codon),
            size = 1,
            colour = "#377EB8") +
  geom_line(data = ctc, 
            aes(x = codon_position, y = codon_percentage, group = codon),
            size = 1,
            colour = "#984EA3") +
  geom_line(data = ctt, 
            aes(x = codon_position, y = codon_percentage, group = codon),
            size = 1,
            colour = "#4DAF4A") +
  geom_line(data = cta, 
            aes(x = codon_position, y = codon_percentage, group = codon),
            size = 1,
            colour = "#FF7F00") +
  geom_line(data = ttg, 
            aes(x = codon_position, y = codon_percentage, group = codon),
            size = 1,
            colour = "#BF5B17") +
  geom_line(data = tta, 
            aes(x = codon_position, y = codon_percentage, group = codon),
            size = 1,
            colour = "#F0027F") +
  facet_wrap(~mRNAtype) +
  ylab("Codon frequency (%)") +
  xlab("Codon Position") +
  publication_theme()

gplot_RSCU <- 
  ggplot() +
  geom_line(data = background, 
            aes(x = codon_position, y = RSCU, group = codon),
            size = 1,
            colour = "#C8C8C8") +
  geom_line(data = ctg, 
            aes(x = codon_position, y = RSCU, group = codon),
            size = 1,
            colour = "#377EB8") +
  geom_line(data = ctc, 
            aes(x = codon_position, y = RSCU, group = codon),
            size = 1,
            colour = "#984EA3") +
  geom_line(data = ctt, 
            aes(x = codon_position, y = RSCU, group = codon),
            size = 1,
            colour = "#4DAF4A") +
  geom_line(data = cta, 
            aes(x = codon_position, y = RSCU, group = codon),
            size = 1,
            colour = "#FF7F00") +
  geom_line(data = ttg, 
            aes(x = codon_position, y = RSCU, group = codon),
            size = 1,
            colour = "#BF5B17") +
  geom_line(data = tta, 
            aes(x = codon_position, y = RSCU, group = codon),
            size = 1,
            colour = "#F0027F") +
  facet_wrap(~mRNAtype) +
  ylab("RSCU") +
  xlab("Codon Position") +
  publication_theme()


SP_containing_transcripts <- longest %>% 
  filter(signal_peptide == T) %>% 
  pull(sequence)

noSP_containing_transcripts <- longest %>% 
  filter(signal_peptide == F) %>% 
  pull(sequence)

RSCUs <- function(vec_of_seqs, tag){
  
  all_RSCUs <- lapply(str_split(vec_of_seqs, ""),
                      function(x){
                        
                        vec <- x[1:90]
                        return(seqinr::uco(vec, as.data.frame = T, NA.rscu = 0))
                        
                      })
  
  all_RSCUs2 <-  
    do.call("rbind", lapply(1:length(all_RSCUs), function(a){all_RSCUs[[a]] %>% mutate(transcript = a)})) %>% 
    mutate(mRNA_type = tag)
  
  return(all_RSCUs2)
  
}

RSCU_data_final <- 
  rbind(RSCUs(SP_containing_transcripts, tag = "Transcripts\nwith SP"),
        RSCUs(noSP_containing_transcripts, tag = "Transcripts\nwithout SP")) 


RSCU_data_final$mRNA_type <- factor(RSCU_data_final$mRNA_type,
                                    levels = c("Transcripts\nwithout SP" ,"Transcripts\nwith SP"))
ctg_RSCU <- 
  RSCU_data_final %>% 
  filter(codon == "ctg") %>% 
  ggplot(aes(x = mRNA_type, y = RSCU, fill = mRNA_type)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  #ggtitle("Gallus gallus") +
  scale_fill_manual(values = c("#FF7F0EFF", "#1F77B4FF")) +
  stat_summary(fun.y=mean, geom="point", shape=15, size=1.5, color="black", fill="black") +
  ylab("CTG RSCU (Codons 0-30)") +
  publication_theme() +theme(axis.title.x = element_blank(), legend.position = "none")

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/RSCU_dog.png",
    res = 300, height = 1000, width = 1000)
print(ctg_RSCU)
dev.off()

wilcox.test(RSCU~mRNA_type,
            RSCU_data_final %>% 
              filter(codon == "ctg"))

RSCU_data_final %>% 
  filter(!codon %in% c("atg", "tga", "tag", "taa")) %>% 
  ggplot(aes(x = noSP, y = SP, label = codon)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel() +
  geom_abline(lty = "dashed", size = 1) +
  publication_theme()
