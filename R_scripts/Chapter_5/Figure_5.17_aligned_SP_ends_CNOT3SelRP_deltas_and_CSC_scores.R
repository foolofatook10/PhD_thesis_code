library(tidyverse)

# load in script defining the border for the ~800 SP containing transcripts in CNOT3 SelRP dataset
borders <- 
  read_tsv("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_SPs_regions.tsv", col_names = F) %>% 
  select(X1 : X5)

# tidy up the columns
colnames(borders) <- c("ENST", "SigP", "region", "start", "end")

# extract the very end codons of SPs
borders_end <- 
  borders %>% filter(region == "c-region") %>% 
  select(ENST, end)

# Load in CPMs
CPMS <- as_tibble(fread("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/CDS_CPMs.csv"))

# Calculate mean CPM per codon position
CPMS$summed_CPM_codon <- (CPMS$summed_CPM_codon/3) 

# Calculate mean CPM per codon position across replicates
mean_CPMs <- 
CPMS %>% 
  filter(MD30 == "M") %>% 
  group_by(IP, transcript, codon) %>% 
  summarise(CPM = mean(summed_CPM_codon))

# C3_enriched <- read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/Analysis/DESeq2_output/RPFs_CNOT3_DEseq2_apeglm_LFC_shrinkage_M_28nts.csv")
# 
# C3_enriched_transcripts <- 
#   C3_enriched %>% mutate(enrichment = case_when(padj < 0.05 & log2FoldChange > 0 ~ "enriched",
#                                                 .default = "not_enriched")) %>% 
#   filter(enrichment == "enriched") %>% 
#   pull(transcript)

# Calculate the deltas
deltas <- 
  CPMS %>% filter(MD30 == "M") %>% 
  select(replicate, IP, transcript, codon, summed_CPM_codon) %>% 
  spread(key = IP, value = summed_CPM_codon) %>% 
  mutate(delta = CNOT3 - Tot) 

# calculate the mean deltas
mean_deltas <- 
deltas %>% 
  group_by(transcript, codon) %>% 
  summarise(mean_delta = mean(delta))

# filter deltas for jsut the 800 transcripts we're interested in
SP_borders_deltas <- 
mean_deltas %>% 
  filter(transcript %in% borders_end$ENST) 

# join with end codons of SPs
SP_borders_deltas <-   
SP_borders_deltas %>% 
  inner_join(borders_end, by = c("transcript" = "ENST"))

SP_borders_deltas_split <- split(SP_borders_deltas, SP_borders_deltas$transcript)

add_positions <- function(x){
  
  data <- SP_borders_deltas_split[[x]]
  
  end <- data$end[1]
  relative_position_start <- (-1*end)
  relative_position_end <- (max(data$codon) - end - 1)
  relative_position_vector <- relative_position_start:relative_position_end
  
  data$relative_position <- relative_position_vector
  return(data)
}

# centre the data at the end codon of the SP
data_centred <- 
do.call("rbind",
        lapply(1:length(SP_borders_deltas_split), add_positions))

#Load theme
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
            axis.title = element_text(size = 16),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = 14), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

aligned_ends <- 
data_centred %>% 
  filter(relative_position >= -10 & relative_position <= 300) %>% 
  group_by(relative_position) %>% 
  summarise(delta = mean(mean_delta)) %>% 
  ggplot(aes( x = relative_position, y = delta)) +
  geom_line(size = 1, colour = "#0C7BDC") +
  xlab("Relative position to SP end (codons)") +
  ylab("Delta CPM\n(CNOT3 IP - Total)") +
  geom_vline(xintercept = 41, lty = "dashed", size = 0.75, colour = "#606060") +
  scale_x_continuous(limits = c(-10, 300),
                     breaks = c(1, 41, 100, 150, 200, 250, 300),
                     labels = c("+1", "+41", "+100", "+150", "+200", "+250", "+300")) +
  publication_theme()

data_centred$transcript %>% unique()

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/aligned_ends_SP.png",
    res = 300, height =  1000, width = 1500)
print(aligned_ends)
dev.off()

CSCs <- 
read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/get_shit_done/chapter5/SGs_CSC_scores.csv") %>% 
  select(codon, CSC = ctrl_cor_estimate) 

master_long <- as_tibble(fread("\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/master_long.csv"))


proline_tryout <- 
data_centred %>%
  filter(relative_position >= 1 & relative_position <= 41) %>%
  inner_join(master_long %>% select(ENST, codon_position, sequence, AA, Properties),
             by = c("transcript" = "ENST", "codon" = "codon_position")) %>%
  dplyr::rename(codon_position = codon) 

proline_transcripts <- 
proline_tryout %>% 
  filter(relative_position < 5 & AA == "Pro") %>% pull(transcript) %>% unique()

high_pro <- 
proline_tryout %>% 
  filter(transcript %in% proline_transcripts) %>% 
  filter(relative_position <= 10) %>% 
  group_by(transcript, AA) %>% 
  summarise(n_pro = n()) %>% filter(AA == "Pro") %>% 
  arrange(-n_pro) %>% 
  filter(n_pro >= 3) %>% pull(transcript)

low_pro <- 
  proline_tryout %>% 
  filter(transcript %in% proline_transcripts) %>% 
  filter(relative_position <= 10) %>% 
  group_by(transcript, AA) %>% 
  summarise(n_pro = n()) %>% filter(AA == "Pro") %>% 
  arrange(-n_pro) %>% 
  filter(n_pro < 3) %>% pull(transcript)

no_pro <- 
proline_tryout %>% 
  filter(!transcript %in% high_pro) %>% 
  filter(!transcript %in% low_pro) %>% 
  pull(transcript) %>% unique()

some_pro <- c(low_pro, high_pro)

length(no_pro) + length(low_pro) +length(high_pro)

data_centred %>% 
  mutate(proline_content = case_when(transcript %in% high_pro ~ ">= 3 residues",
                                     transcript %in% low_pro ~ "1-2 residues",
                                     .default = "0 residues")) %>% 
  filter(relative_position >= -10 & relative_position <= 300) %>% 
  group_by(relative_position, proline_content) %>% 
  summarise(delta = mean(mean_delta)) %>% 
  ggplot(aes( x = relative_position, y = delta, colour = proline_content)) +
  geom_line(size = 1) +
  xlab("Relative position to SP end (codons)") +
  ylab("Delta CPM\n(CNOT3 IP - Total)") +
  geom_vline(xintercept = 41, lty = "dashed", size = 0.75, colour = "#606060") +
  scale_x_continuous(limits = c(-10, 300),
                     breaks = c(1, 41, 100, 150, 200, 250, 300),
                     labels = c("+1", "+41", "+100", "+150", "+200", "+250", "+300")) +
  publication_theme() 


length(CSCs$codon)

CSC_final <-
data_centred %>% 
  filter(transcript %in% some_pro) %>% 
  filter(relative_position >= -10 & relative_position <= 300) %>% 
  inner_join(master_long %>% select(ENST, codon_position, sequence, AA, Properties),
             by = c("transcript" = "ENST", "codon" = "codon_position")) %>% 
  dplyr::rename(codon_position = codon) %>% 
  inner_join(CSCs %>% dplyr::rename(sequence = "codon"), by = "sequence") %>% 
  group_by(relative_position) %>% 
  summarise(mean_CSC = mean(CSC))  

CSC_scores_gplot <- 
CSC_final %>% 
  filter(relative_position >= 1 & relative_position <= 81) %>% 
ggplot(aes(x = relative_position, y = mean_CSC)) +
  geom_line() +
  geom_point() +
  geom_smooth() +
  publication_theme() +
  xlab("Codon position\nrelative to SP end") +
  ylab("Mean CSC score") +
  scale_x_continuous(breaks = c(1,20,41,60,80),
                     labels = c("+1", "+20", "+41", "+60", "+80")) 

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/aligned_SPs_ends_CSC_scores.png",
    res = 300, height = 1000, width = 1000)
print(CSC_scores_gplot)
dev.off()

bxplot_data <- 
CSC_final %>% 
  filter(relative_position >= 1 & relative_position <= 81) %>% 
  mutate(region = case_when(relative_position <= 41 ~ "before",
                            #relative_position >20 & relative_position <= 40 ~ "middle",
                            relative_position > 40 ~ "after"))  

bxplot_data$region <- factor(bxplot_data$region, 
                             levels = c("before", "after"),
                             labels = c("Codons\n0-40", "Codons\n41-81"))

boxplot_CSC_gplot <- 
bxplot_data %>% 
  ggplot(aes(x = region, y = mean_CSC)) +
  geom_violin() +
  geom_point(position = position_jitter(width = 0.1)) +
  geom_boxplot(width = 0.1) +
  ylab("Codon position CSC score") +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 18))

t.test(mean_CSC ~ region, data = bxplot_data)
wilcox.test(mean_CSC ~ region, data = bxplot_data)

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/boxplot_CSC_aligned_SPs.png",
    res = 300, height = 1000, width = 1000)
print(boxplot_CSC_gplot)
dev.off()

sequence_features <- 
data_centred %>%
  #filter(transcript %in% some_pro) %>% 
  filter(relative_position >= 1 & relative_position <= 100) %>% 
  inner_join(master_long %>% select(ENST, codon_position, sequence, AA, Properties),
             by = c("transcript" = "ENST", "codon" = "codon_position")) %>% 
  dplyr::rename(codon_position = codon)

data_centred%>%
  filter(transcript %in% some_pro) %>% pull(transcript) %>% unique() %>% length()

sequence_features %>% 
  group_by(relative_position, Properties) %>% 
  summarise(n_property = n()) %>% 
  mutate(property_percentage = n_property/778 * 100) %>% 
  ggplot(aes(x = relative_position, y = property_percentage, colour = Properties)) +
  geom_line()  +
  geom_vline(xintercept = 41)

sequence_features %>% 
  group_by(relative_position, AA) %>% 
  summarise(n_AA = n()) %>% 
  mutate(AA_percentage = n_AA/778 * 100) %>% 
  ggplot(aes(x = relative_position, y = AA_percentage)) +
  geom_line() +
  facet_wrap(~AA)+
  geom_vline(xintercept = 41)

sequence_features %>% 
  mutate(wobble = case_when(str_detect(sequence, pattern = "a$") ~ "a",
                            str_detect(sequence, pattern = "t$") ~ "t",
                            str_detect(sequence, pattern = "g$") ~ "g",
                            str_detect(sequence, pattern = "c$") ~ "c")) %>% 
  group_by(relative_position, wobble) %>% 
  summarise(n_codons = n()) %>% 
  mutate(wobble_percentage = n_codons/334 * 100) %>% 
  ggplot(aes(x = relative_position, y = wobble_percentage, colour = wobble)) +
  geom_line() +
  geom_smooth()+
  geom_vline(xintercept = 41) +
  scale_y_continuous(limits = c(0,NA))

sequence_features %>% 
  filter(relative_position >=20 & relative_position <= 60) %>% 
  group_by(relative_position, sequence) %>% 
  summarise(n_codons = n()) %>% 
  mutate(codon_percentage = n_codons/334 * 100) %>% 
  ggplot(aes(x = relative_position, y = codon_percentage)) +
  facet_wrap(~sequence) +
  geom_line() 

box_types <- 
read_csv("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/useful_tables/codon_box_types.csv") %>% 
  mutate(codon = str_replace_all(tolower(codon), "u", "t")) %>% 
  mutate(box = as.numeric(str_remove(case_when(AA == "Ile" ~ "3-box", .default = box), "-box")))

RSCU_data <- 
inner_join(
sequence_features %>% 
  group_by(relative_position, AA, sequence) %>% 
  summarise(n_codons = n()),
sequence_features %>% 
  group_by(relative_position, AA) %>% 
  summarise(total_residues = n()),
by = c("relative_position", "AA")) %>% 
  mutate(codon_freq = n_codons/total_residues) %>% 
  inner_join(box_types, by = c("sequence" = "codon", "AA")) %>% 
  mutate(RSCU = codon_freq * box) %>% 
  mutate(wobble = case_when(str_detect(sequence, pattern = "a$") ~ "a",
                            str_detect(sequence, pattern = "t$") ~ "t",
                            str_detect(sequence, pattern = "g$") ~ "g",
                            str_detect(sequence, pattern = "c$") ~ "c"))

 
RSCU_data %>% 
  group_by(relative_position, wobble) %>% 
  summarise(mean_RSCU = mean(RSCU)) %>% 
  filter(relative_position >= 20 & relative_position <= 60) %>% 
  ggplot(aes(x = relative_position, y = mean_RSCU, colour = wobble)) +
  geom_line()

RSCU_data %>% 
  filter(AA == "Leu") %>% 
  ggplot(aes(x = relative_position, y = RSCU, colour = sequence)) +
  geom_line()

RSCU_data %>% 
  filter(AA == "Arg") %>% 
  ggplot(aes(x = relative_position, y = RSCU, colour = sequence)) +
  geom_line()

RSCU_data %>% 
  filter(AA == "Pro") %>% 
  ggplot(aes(x = relative_position, y = RSCU, colour = sequence)) +
  geom_line()

RSCU_data %>% 
  filter(AA == "Glu") %>% 
  ggplot(aes(x = relative_position, y = RSCU, colour = sequence)) +
  geom_line()

RSCU_data %>% 
  filter(AA == "Ala") %>% 
  ggplot(aes(x = relative_position, y = RSCU, colour = sequence)) +
  geom_line()




