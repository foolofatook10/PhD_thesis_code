library(tidyverse)
library(data.table)
library(seqinr)
library(Biostrings)

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
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

aspergillus <- 
  read.fasta("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/aspergillus_nidulans/Aspergillus_nidulans.ASM1142v1.cds.all.fa", whole.header = T)

sequences <- 
  unlist(getSequence(aspergillus, as.string = T))
names <- unlist(getName.list(aspergillus))

get_transcripts <- function(x){
  
  #str_extract(x, pattern = "^CAJ\\d+\\.\\d+")
  lapply(str_split(x, " "), function(x){x[[1]]})
  
}

get_genes <- function(x){
  
  str_extract(x, pattern = "ANIA_\\d+")
  #lapply(str_split(x, " "), function(x){x[[1]]})
}


transcripts <- unlist(lapply(names, get_transcripts))
genes <- unlist(lapply(names, get_genes))

length(names)
length(unique(names))

aspergillus_sequences <- 
  data.frame(gene = genes,
             transcript = transcripts,
             sequence = sequences)


aspergillus_sequences_atg <- 
  aspergillus_sequences %>% 
  filter(str_detect(sequence, pattern = "^atg"))

aspergillus_sequences_stops <- 
  aspergillus_sequences_atg %>% 
  filter(str_detect(sequence, pattern = "(taa|tag|tga)$"))

aspergillus_sequences_correct_length <- 
  aspergillus_sequences_stops %>% 
  filter(nchar(sequence)%%3 == 0)

final_polished_sequences <- 
  aspergillus_sequences_correct_length %>% 
  filter(!str_detect(sequence, pattern = "[^atgc]"))

length(final_polished_sequences$gene %>% unique())

longest <-
  final_polished_sequences %>% 
  mutate(cds_length = nchar(sequence)) %>%
  group_by(gene) %>%
  mutate(the_rank  = rank(-cds_length, ties.method = "first")) %>%
  filter(the_rank == 1) %>% 
  select(-the_rank)%>% 
  filter(cds_length > 300)

AA_sequences <- as.vector(translate(Biostrings::DNAStringSet(longest$sequence)))

longest$AA_sequence <- tolower(AA_sequences)

sequences_split <- str_split(longest$AA_sequence, "")
names_fasta <- longest$transcript
length(names_fasta)

write.fasta(sequences_split[1:5000], names = names_fasta[1:5000], file.out = "C:/Users/Jettles/OneDrive - University of Glasgow/Documents/aspergillus_1_5000.fa")
write.fasta(sequences_split[5001:9923], names = names_fasta[5001:9923], file.out = "C:/Users/Jettles/OneDrive - University of Glasgow/Documents/aspergillus_5000_10000.fa")



fyles <- paste0("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/aspergillus_nidulans/output_matureAN", 
                c("1", "2"),
                ".fasta")

transcripts_with_SP <- unlist(lapply(lapply(fyles, read.fasta),getName))

transcripts_annotated <- 
  longest %>% 
  select(gene, transcript, sequence) %>% 
  mutate(signal_peptide = case_when(transcript %in% transcripts_with_SP ~ T,
                                    .default = F))

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

# make a list with the nucleotide sequences chopped into codons
my_lyst <- sapply(transcripts_annotated$sequence, split_into_codons)

# Name each element of the list by its ENST tag
names(my_lyst) <- pull(transcripts_annotated, var = transcript)

# obtain logical vector of CDSs that are greater than 100 codons
over_100_codons <- lapply(my_lyst, length) > 100

# filter my_lyst for CDSs that are greater than 100 codons
my_lyst <-  my_lyst[c(over_100_codons)]

# only consider the first 50 codons per sequence
codonsN <-lapply(my_lyst, head, n=100)

# condense list into a tbale format to make it easier to work with
codonsN_table <- data.frame(do.call("rbind", codonsN))

max_n <- max(as.numeric(str_remove(colnames(codonsN_table), "X")))

# label columns with numbers from 1 to 50
colnames(codonsN_table) <- c(paste0("codon",1:max_n))

# convert rownames to column with name transcript
codonsN_table <- 
  codonsN_table %>% rownames_to_column("transcript")

# convert table to long format
codonsN_table_long <- 
  codonsN_table %>% 
  gather(key = "position", value =  "codon", codon1:paste0("codon", max_n))


codons <- 
  rbind(codons,
        data.frame(codon = c("taa", "tga", "tag"), AA = rep("Stp",3), box = "3-box"))

# add amino acid information and signal sequence information to the table
codonsN_table_long <- 
  codonsN_table_long %>%
  inner_join(codons %>% select(codon, AA), by = "codon") %>%
  inner_join(transcripts_annotated %>% select(transcript, signal_peptide), by = "transcript")

presence_SP <- 
  codonsN_table_long %>% 
  select(transcript, signal_peptide) %>% 
  unique() %>% 
  group_by(signal_peptide) %>% 
  tally()

signal_peptides_number = presence_SP %>% 
  filter(signal_peptide == TRUE) %>% 
  pull(n) 
no_signal_peptides_number = presence_SP %>% 
  filter(signal_peptide == FALSE) %>% 
  pull(n) 

# convert positions to numerics
codonsN_table_long$position <- as.numeric(str_remove(codonsN_table_long$position, pattern = "codon"))

raw_counts <- 
  codonsN_table_long %>%
  group_by(signal_peptide, position) %>%
  dplyr::count(vars = codon) %>% 
  ungroup()

# calculate frequency of codon usage per position
percentages <- cbind(raw_counts,
                     raw_counts %>%
                       group_by(signal_peptide, position) %>%
                       summarise(percent = n/sum(n) * 100) %>%
                       ungroup() %>%
                       select(percent))

# add column for wobble base position
percentages <- percentages %>%
  dplyr::rename(codon = vars) %>%
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A/U", "A/U", "G/C","G/C")))

# add tags
temp <-
  percentages %>% 
  mutate(tag = case_when(signal_peptide == F ~ paste0("No signal sequence", 
                                                      "\nn = ", 
                                                      no_signal_peptides_number),
                         .default = paste0("Signal sequence", "\nn = ", signal_peptides_number)))

# vector of leucyl codons
leucyl_codons <- c(codons %>% filter(AA == "Leu") %>% pull(var = codon))

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
                                   .default = "other")) %>% 
  ungroup()

# make leucyl codons column a factor
lines_processed$leucyl <- factor(lines_processed$leucyl, levels = c("ctg", "ctc", "ctt", "ttg", "cta", "tta","other"))

background <- lines_processed %>% filter(leucyl == "other" & position != 1 &codon != "gcc")
background_gcc <- lines_processed %>% filter(codon == "gcc" & position != 1)
leucyl2plot <- lines_processed %>% filter(leucyl %in% c("ctg", "ctc", "ctt", "ttg", "cta", "tta") & position != 1)

leucyl_colours = c("#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F")

final_gplot <- 
  ggplot() +
  geom_line(data = background,
            aes(x = position, y = percent, group = codon),
            size = 1,
            colour = "#C8C8C8") +
  geom_line(data = background_gcc,
            aes(x = position, y = percent, group = codon),
            size = 1,
            colour = "#4DBBD5FF") +
  geom_line(data = leucyl2plot,
            aes(x = position, y = percent, group = codon, colour = leucyl),
            size = 1) +
  ylim(0,20) +
  xlab("Codon Position") +
  ylab("Percent") +
  scale_color_manual(values = leucyl_colours) +
  ggtitle("Aspergillus nidulans") +
  facet_wrap(~tag) +
  coord_cartesian(xlim = c(0,50)) +
  scale_x_continuous(limits = c(2,50), breaks = c(2,10,20,30,40,50)) +
  publication_theme() +
  theme(legend.position = "none", 
        legend.direction = "vertical", 
        legend.title = element_blank())

# ggsci::pal_npg("nrc")(4)
# "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF"

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/A_nidulans.png",
    res = 300,
    height = 1000,
    width = 1500)
print(final_gplot)
dev.off()


