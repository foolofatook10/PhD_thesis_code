library(data.table)
library(tidyverse)

# path to signal sequence information
transcripts_annotated_by_SigP <- "\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/features.csv"
transcripts_nucleotide_sequences <- "\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/master.csv"

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
            axis.text.x = element_text(size = 12, color = "black"),
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
master <- as_tibble(fread(file = transcripts_nucleotide_sequences, header = T, drop = "V1"))

# select the longest transcript per gene
longest <- 
  master %>% 
  mutate(polypeptide_length = nchar(master$nucleotide_sequence_CDS)/3) %>%
  group_by(ENSG) %>%
  mutate(the_rank  = rank(-polypeptide_length, ties.method = "first")) %>%
  filter(the_rank == 1) %>% 
  select(-the_rank)

features <- read_csv(transcripts_annotated_by_SigP) %>% 
  select(ENSG, ENST, signal_sequence) 

longest <- inner_join(longest, features, by = c("ENSG", "ENST"))


longest <- longest %>% 
  dplyr::rename(signal_peptide = signal_sequence) %>% 
  dplyr::rename(sequence = nucleotide_sequence_CDS)

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

RSCU_comparisons <- 
  rbind(signal_sequence_information, cytosolic_information) %>% 
  filter(codon_position != 1) %>%
  filter(codon != "atg")  


codons_order <- c("ctg", "ctc", "cta", "ctt", "ttg", "tta",
                  "gtg", "gtc", "gta", "gtt",
                  "gcg", "gcc", "gct", "gca",
                  "atc", "ata", "att",
                  "ttc", "ttt",
                  "tgg",
                  "tac", "tat",
                  "tcg", "tcc", "tca", "tct", "agc", "agt",
                  "acg", "acc", "aca", "act",
                  "aac", "aat",
                  "cag", "caa",
                  "tgc", "tgt",
                  "ggg", "ggc", "gga", "ggt",
                  "ccg", "ccc", "cca", "cct",
                  "cgg", "cgc", "cga", "cgt", "agg", "aga",
                  "cac", "cat",
                  "aag", "aaa",
                  "gac", "gat",
                  "gaa", "gag")

#codons %>% filter(codon == "cac")


labels_for_codons <- c("Leu-ctg", "Leu-ctc", "Leu-cta", "Leu-ctt", "Leu-ttg", "Leu-tta",
                       "Val-gtg", "Val-gtc", "Val-gta", "Val-gtt",
                       "Ala-gcg", "Ala-gcc", "Ala-gct", "Ala-gca",
                       "Ile-atc", "Ile-ata", "Ile-att",
                       "Phe-ttc", "Phe-ttt",
                       "Trp-tgg",
                       "Tyr-tac", "Tyr-tat",
                       "Ser-tcg", "Ser-tcc", "Ser-tca", "Ser-tct", "Ser-agc", "Ser-agt",
                       "Thr-acg", "Thr-acc", "Thr-aca", "Thr-act",
                       "Asn-aac", "Asn-aat",
                       "Gln-cag", "Gln-caa",
                       "Cys-tgc", "Cys-tgt",
                       "Gly-ggg", "Gly-ggc", "Gly-gga", "Gly-ggt",
                       "Pro-ccg", "Pro-ccc", "Pro-cca", "Pro-cct",
                       "Arg-cgg", "Arg-cgc", "Arg-cga", "Arg-cgt", "Arg-agg", "Arg-aga",
                       "His-cac", "His-cat",
                       "Lys-aag", "Lys-aaa",
                       "Asp-gac", "Asp-gat",
                       "Glu-gaa", "Glu-gag")


RSCU_comparisons$codon <- factor(RSCU_comparisons$codon, 
                                 levels = codons_order,
                                 labels = labels_for_codons)

RSCU_comparisons_gplot <- 
  RSCU_comparisons %>% 
  ggplot(aes(x = codon_position, y = RSCU, group = mRNAtype, colour = mRNAtype)) +
  geom_line(size = 1) +
  facet_wrap(~codon, nrow = 6, ncol = 10) +
  scale_color_manual(values = c("#FF7F0EFF", "#1F77B4FF")) +
  xlab("Codon position") +
  publication_theme() +
  theme(legend.title = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/all_codons_in SPs_RSCU.png",
    res= 300, height = 2500, width = 2500)
print(RSCU_comparisons_gplot)
dev.off()


## Figure B
SP_containing_transcripts <- longest %>% 
  filter(polypeptide_length > 100) %>% 
  filter(signal_peptide == T) %>% 
  pull(sequence)

noSP_containing_transcripts <- longest %>% 
  filter(polypeptide_length > 100) %>% 
  filter(signal_peptide == F) %>% 
  pull(sequence)

RSCUs <- function(vec_of_seqs, tag, nt_start, nt_end){
  
  all_RSCUs <- lapply(str_split(vec_of_seqs, ""),
                      function(x){
                        
                        vec <- x[nt_start:nt_end]
                        return(seqinr::uco(vec, as.data.frame = T, NA.rscu = 0))
                        
                      })
  
  all_RSCUs2 <-  
    do.call("rbind", lapply(1:length(all_RSCUs), function(a){all_RSCUs[[a]] %>% mutate(transcript = a)})) %>% 
    mutate(mRNA_type = tag)
  
  return(all_RSCUs2)
  
}

RSCU_data_final <- 
  rbind(RSCUs(SP_containing_transcripts, tag = "Transcripts\nwith SP", nt_start = 1, nt_end = 90),
        RSCUs(noSP_containing_transcripts, tag = "Transcripts\nwithout SP", nt_start = 1, nt_end = 90)) 


RSCU_data_final$mRNA_type <- factor(RSCU_data_final$mRNA_type,
                                    levels = c("Transcripts\nwithout SP" ,"Transcripts\nwith SP"))

val_gtg <-
  RSCU_data_final %>%
  filter(codon == "gtg") %>%
  ggplot(aes(x = mRNA_type, y = RSCU, fill = mRNA_type)) +
  geom_violin() +
  geom_boxplot(width = 0.05) +
  #facet_wrap(~codon) +
  #ggtitle("Gallus gallus") +
  scale_fill_manual(values = c("#FF7F0EFF", "#1F77B4FF")) +
  #ggtitle("CGG") +
  stat_summary(fun.y=mean, geom="point", shape=15, size=1.5, color="black", fill="black", position = position_dodge(0.5)) +
  ylab("RSCU (Codons 0-30)") +
  publication_theme() +theme(axis.title.x = element_blank(), legend.position = "none",
                             axis.text.x = element_text(size = 16))
# 
val_gtc <-
  RSCU_data_final %>%
  filter(codon == "gtc") %>%
  ggplot(aes(x = mRNA_type, y = RSCU, fill = mRNA_type)) +
  geom_violin() +
  geom_boxplot(width = 0.05) +
  #facet_wrap(~codon) +
  #ggtitle("Gallus gallus") +
  scale_fill_manual(values = c("#FF7F0EFF", "#1F77B4FF")) +
  #ggtitle("CGG") +
  stat_summary(fun.y=mean, geom="point", shape=15, size=1.5, color="black", fill="black", position = position_dodge(0.5)) +
  ylab("RSCU (Codons 0-30)") +
  publication_theme() +theme(axis.title.x = element_blank(), legend.position = "none",
                             axis.text.x = element_text(size = 16))