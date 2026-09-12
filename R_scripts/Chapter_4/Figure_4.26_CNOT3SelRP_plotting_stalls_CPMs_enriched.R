library(tidyverse)
library(data.table)

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
            legend.title = element_text(face="italic", size = 14),
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

codons_deseq <- 
read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/Analysis/DESeq2_output/M_RPFs_28nts.csv") 

master <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/master.csv") %>% 
  select(Gene, ENSG, ENST, nucleotide_sequence_CDS)
  
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

master <- master %>% filter(ENST %in% unique(codons_deseq$transcript))

# split the sequences from master
split_sequences <- lapply(master$nucleotide_sequence_CDS, split_into_codons)

# apply transcript names to list
transcript_names <- master$ENST

names(split_sequences) <- transcript_names

tidy_data <- function(x){
  
  
  transcript <- transcript_names[[x]]
  codon_sequence <- split_sequences[[x]]
  
  df <- 
  data.frame(transcript = transcript,
             codon = codon_sequence,
             position = 1:length(codon_sequence))
  
  return(df)
}

tidy_data <- lapply(1:length(split_sequences), tidy_data)
names(tidy_data) <- transcript_names

codons_deseq <- 
  codons_deseq %>% filter(transcript %in% c(master$ENST))

enriched <- codons_deseq %>% filter(log2FoldChange > 0 & padj < 0.05)


extract_regions <- function(a){
  
  ENST <- enriched[a,]$transcript 
  
  print(a)
  
  position <- enriched[a,]$codon
  
  raw <- tidy_data[[{{ENST}}]]
  
  total_length <- nrow(raw)
  
  
  
  start <- (-1 * position) + 1 
  end <- total_length - position
  
  df <- 
  raw %>% 
    mutate(relative_position = start:end)
  
  df$initial_position <- position
  
  
  #print(df)
  
  return(df)
}




# enriched_split <- split(enriched, enriched$transcript)
# 
# # filter enriched split for the greatest log2FoldChange enrichment
# filter_greatest_log2 <- function(y){
#   
#   to_return <- y %>% filter(log2FoldChange == max(log2FoldChange))
#   return(to_return)
# }
# 
# enriched_split <- lapply(enriched_split, filter_greatest_log2)
#enriched <- do.call("rbind", enriched_split)


test <- lapply(1:nrow(enriched), extract_regions)

filter_slice <- function(b){
  
  b %>% 
    filter(relative_position >= -60 & relative_position <= 60)
  
}

final_data <- 
do.call("rbind", lapply(test, filter_slice))

codons_raw <- 
final_data %>% 
  group_by(relative_position, codon) %>% 
  tally() %>% 
  ungroup()

codons_raw_sum <- 
  final_data %>% 
  group_by(relative_position) %>% 
  tally() %>% 
  dplyr::rename(sum_n = n) %>% 
  ungroup()


codons_table <- read_csv("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/useful_tables/codon_amino_acid_pairs_and_properties.csv")

#Make the capitalised codons lower case
codons_table$codon <- factor(tolower(codons_table$codon))

#replace letters "u" with "t"
codons_table$codon <- gsub("u","t", codons_table$codon)


annotated_data <- 
inner_join(codons_raw, codons_raw_sum, by = "relative_position") %>% 
  mutate(percentage = n/sum_n * 100)%>% 
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A", "U", "G","C"))) %>% 
  inner_join(codons_table %>% select(codon, AA, Properties),
             by = "codon")

nucleotides_gplot <- 
annotated_data %>% 
  group_by(relative_position, wobble) %>% 
  summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = wobble)) +
  geom_line() +
  publication_theme() +
  geom_point() +
  ylab("Codon percentage") +
  xlab("Relative position (P site = 0)") +
  #coord_cartesian(xlim = c(-20,20)) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0,NA))

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/wobble.png",
    res = 300, height = 1000, width = 1000)
print(nucleotides_gplot)
dev.off()


AA_properties <- 
annotated_data %>% 
  group_by(relative_position, Properties) %>% 
  summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = Properties)) +
  geom_point(alpha = 0.5) +
  geom_line(size = 1, alpha = 0.5) +
  scale_y_continuous(limits = c(0,NA))  +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_color_brewer(palette = "Set1") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  publication_theme() +
  theme(legend.position = "right", legend.direction = "vertical")

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/codons_deseq_AAproperties.png",
    res = 300, height = 1000, width = 1200)
print(AA_properties)
dev.off()

AA_frequencies <- 
annotated_data %>% 
  group_by(relative_position, AA) %>% 
  summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage)) +
  geom_line(alpha = 0.5) +
  #geom_point() +
  scale_y_continuous(limits = c(0,NA))  +
  #coord_cartesian(xlim = c(-10,10)) +
  facet_wrap(~AA) +
  scale_color_brewer(palette = "Set1") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)")

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/codons_deseq_AA_individual.png",
    res = 300, height = 1500, width = 2000)
print(AA_frequencies)
dev.off()

# all codons
annotated_data %>% 
  ggplot(aes(x = relative_position, y = percentage)) +
  geom_line(alpha = 0.5) +
  #geom_point() +
  scale_y_continuous(limits = c(0,NA))  +
  #coord_cartesian(xlim = c(-10,10)) +
  facet_wrap(~codon) +
  scale_color_brewer(palette = "Set1") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)")





# Proline
Proline <- 
annotated_data %>%
  filter(codon %in% c("ccg", "ccc", "cca", "cct")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  ggtitle("Proline") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  publication_theme() +
  theme(legend.title = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/codons_deseq_Proline.png", res = 300, height = 1000, width = 1000)
print(Proline)
dev.off()

# Glutamate
Glu <- 
annotated_data %>%
  filter(codon %in% c("gaa", "gag")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Glutamate") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/codons_deseq_GLU.png", res = 300, height = 1000, width = 1000)
print(Glu)
dev.off()


# Arginine
Arg <- 
  annotated_data %>%
  filter(codon %in% c("cgg", "cga", "cgt", "cgc", "agg", "aga")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Arginine") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/codons_deseq_ARG.png", res = 300, height = 1000, width = 1000)
print(Arg)
dev.off()

# Leucine
Leu <- 
  annotated_data %>%
  filter(codon %in% c("ctg", "ctc", "ctt", "cta", "ttg", "tta")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Leucine") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/codons_deseq_LEU.png", res = 300, height = 1000, width = 1000)
print(Leu)
dev.off()

# Serine
Ser <- 
  annotated_data %>%
  filter(codon %in% c("tct", "tcc", "tca", "tcg", "agt", "agc")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Serine") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/codons_deseq_LEU.png", res = 300, height = 1000, width = 1000)
print(Ser)
dev.off()


# Valine
Val <- 
  annotated_data %>%
  filter(codon %in% c("gtt", "gtc", "gta", "gtg")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Valine") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/codons_deseq_LEU.png", res = 300, height = 1000, width = 1000)
print(Val)
dev.off()


# Aspartate
Asp <- 
  annotated_data %>%
  filter(codon %in% c("gat", "gac")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Aspartate") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/codons_deseq_Asp.png", res = 300, height = 1000, width = 1000)
print(Asp)
dev.off()



# Lysine
Lys <- 
  annotated_data %>%
  filter(codon %in% c("aag", "aaa")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Lysine") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/codons_deseq_Asp.png", res = 300, height = 1000, width = 1000)
print(Lys)
dev.off()

# Isoleucine
Ile <- 
  annotated_data %>%
  filter(codon %in% c("att", "atc", "ata")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Isoleucine") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())

Phe <- 
  annotated_data %>%
  filter(codon %in% c("ttt", "ttc")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Phenylalanine") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())

Thr <- 
  annotated_data %>%
  filter(codon %in% c("act", "acc", "acg", "aca")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Threonine") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())


His <- 
  annotated_data %>%
  filter(codon %in% c("cat", "cac")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Histidine") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())



Ala <- 
  annotated_data %>%
  filter(codon %in% c("gct", "gcg", "gcc", "gca")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Alanine") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())

Gly <- 
  annotated_data %>%
  filter(codon %in% c("ggt", "ggc", "ggg", "gga")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Glycine") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())

Tyr <- 
  annotated_data %>%
  filter(codon %in% c("tat", "tac")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Tyrosine") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())

Trp <- 
  annotated_data %>%
  filter(codon %in% c("tgg", "atg")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Trp/Met") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())


Gln <- 
  annotated_data %>%
  filter(codon %in% c("caa", "cag")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Glutamine") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())

Asn <- 
  annotated_data %>%
  filter(codon %in% c("aat", "aac")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Asparagine") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())

Cys <- 
  annotated_data %>%
  filter(codon %in% c("tgt", "tgc")) %>% 
  #group_by(relative_position, wobble) %>% 
  #summarise(percentage = sum(percentage)) %>% 
  ggplot(aes(x = relative_position, y = percentage, colour = codon)) +
  geom_line() +
  geom_point() +
  #coord_cartesian(xlim = c(-10,10)) +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme()  +
  ggtitle("Cysteine") +
  ylab("Percentage") +
  xlab("Relative position (P site = 0)") +
  theme(legend.title = element_blank())

tidied_up_data <- do.call("rbind", tidy_data)
rownames(tidied_up_data) <- NULL  

A_site_codons <- codons_deseq %>% mutate(codon = codon + 1)

codons_deseq_tidy <- 
inner_join(tidied_up_data, 
           A_site_codons, 
           by = c("transcript", "position" = "codon")) %>% 
  select(gene_sym, transcript, position, codon, log2FoldChange, padj) %>% 
  na.omit() %>% 
  inner_join(codons_table, by = "codon")
  
split_deseq <- split(codons_deseq_tidy, codons_deseq_tidy$AA)

plot_log2s <- 
function(x){
  
  x %>%
    filter(log2FoldChange > 0 & padj < 0.05) %>% 
    ggplot(aes(x = reorder(codon, -log2FoldChange), y = log2FoldChange)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.2)) +
    #geom_boxplot() +
    ggtitle(paste0(unique(x$AA)))
}

lapply(split_deseq, plot_log2s)


codons_deseq_tidy %>% 
  filter(log2FoldChange > 0 & padj < 0.05) %>% 
  #filter(AA == "Pro") %>% 
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A/U", "A/U", "G/C","G/C"))) %>% 
  select(-Symbol) %>% 
  ggplot(aes(x = wobble, y = log2FoldChange)) +
  #geom_violin() +
  geom_boxplot(width = 0.8) 
  #geom_point(position = position_jitter(width = 0.))

  


#codons_deseq %>% filter(transcript == "ENST00000370873.9")

CPMs <- 
fread("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/CDS_CPMs.csv") %>% 
  filter(MD30 == "M") %>% 
  mutate(CPM = (summed_CPM_codon/3)) %>% 
  group_by(IP, transcript, codon) %>% 
  summarise(mean_CPM = mean(CPM))

CPMs_split <- split(CPMs, CPMs$IP)

join_data <- function(x){
  
  x %>% 
    dplyr::rename(position = codon) %>% 
    inner_join(final_data, 
               by = c("transcript", "position"))
  
}

CPMs_joined <- do.call("rbind", lapply(CPMs_split, join_data))

enriched

CPMs <- 
CPMs_joined %>% 
  group_by(IP, relative_position) %>% 
  summarise(CPM = mean(mean_CPM), sd_CPM = sd(mean_CPM)) %>% 
  ggplot(aes(x = relative_position, y = CPM, colour = IP)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#0C7BDC", "#FFC20A")) +
  ylab("CPM") +
  xlab("Relative position (P site = 0)") +
  scale_y_continuous(limits = c(0,NA)) +
  ggtitle("Codon enriched (n = 2954)") +
  publication_theme() 

#### codons per position ####

codons_per_position <- 
final_data %>% 
  group_by(relative_position) %>% 
  tally() %>% 
  ggplot(aes(x = relative_position, y = n)) +
  geom_line(size = 1) +
  ylab("Number of codons") +
  scale_y_continuous(limits = c(0,NA)) +
  xlab("Relative position\n(P site = 0)") +
  publication_theme()

CPMs_gplot <- 
gridExtra::grid.arrange(CPMs, codons_per_position, nrow =1)

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/codons_deseq_aligned_CPMS.png", res = 300, height = 1000, width = 1000)
print(CPMs)
dev.off()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/codons_deseq_codons_per_position.png", res = 300, height = 1000, width = 1000)
print(codons_per_position)
dev.off()


CPMs <- 
  fread("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/CDS_CPMs.csv") %>% 
  filter(MD30 == "M") %>% 
  mutate(CPM = (summed_CPM_codon/3)) %>% 
  group_by(IP, transcript, codon) %>% 
  summarise(mean_CPM = mean(CPM))

CPMs_transcript <- 
CPMs %>% 
  inner_join(master %>% select(Gene, ENST),
             by = c("transcript" = "ENST"))

CPMs_transcript <- 
CPMs_transcript %>% filter(transcript %in% unique(enriched$transcript))

splitted <- 
split(CPMs_transcript, CPMs_transcript$transcript)

names(splitted) <- unique(CPMs_transcript$Gene)

plot_CPM_gplots <- function(x){
  
  df <- splitted[[x]]
  ENST <- unique(df$transcript)
  
  enriched_position <- 
  enriched %>% 
    filter(transcript == {{ENST}}) %>% 
    pull(codon)
  
  final <- 
  df %>% ggplot(aes(x = codon, y = mean_CPM, colour = IP))+
    geom_line() +
    geom_vline(xintercept = enriched_position, alpha = 0.5) +
    ggtitle(paste0(unique(df$Gene))) +
    publication_theme()
  
  return(final)
  
}

gplots <- lapply(1:length(splitted), plot_CPM_gplots)

names(gplots) <- unique(CPMs_transcript$Gene)

unique(CPMs_transcript$Gene)

gplots$MYC

gplots$CAPRIN1
gplots$EIF2AK1
gplots$TRIM28
gplots$FLRT2
gplots$MLF2
gplots$TRAM2
gplots$AAAS
gplots$NOTCH2
gplots$EEF1B2
gplots$EEF2
gplots$ATF4

gplots$RCN1
gplots$WNK1

gplots$ENO1

gplots$COPE
gplots$IFNGR1
gplots$ITGB1
gplots$G3BP1
gplots$AARS1
gplots$GPI
gplots$SELENOS
gplots$SELENOT

gplots$TMEM30A
gplots$SNAI2
gplots$CHPF2
gplots$GRN
gplots$RCN1
gplots$ADCK2
gplots$CDIPT
gplots$NOMO1
gplots$NDUFB7
gplots$HSP90AA1
gplots$SLC7A5
gplots$EXOSC5

gplots$NCLN
gplots$PSENEN
gplots$CCDC47
gplots$TMCO1
gplots$TRAM1
gplots$TRAM2

gplots$TMEM109
gplots$GORASP2
gplots$SMAD3


