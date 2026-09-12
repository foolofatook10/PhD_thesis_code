library(tidyverse)
library(data.table)

# CNOT3 regulated transcripts
CNOT3_transcripts <- read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/CNOT3_monosomes_transcripts.csv") %>% 
  pull(x)

C3_enriched <- read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/Analysis/DESeq2_output/RPFs_CNOT3_DEseq2_apeglm_LFC_shrinkage_M_28nts.csv")
C3_enriched_transcripts <- 
  C3_enriched %>% mutate(enrichment = case_when(padj < 0.05 & log2FoldChange > 0 ~ "enriched",
                                                .default = "not_enriched")) %>% 
  filter(enrichment == "enriched") %>% 
  pull(transcript)

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
            legend.text = element_text(size = 16),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

CPMs <- 
  fread("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/CDS_CPMs.csv") %>% 
  as_tibble() %>% 
  select(replicate:codon,summed_CPM_codon) %>% 
  dplyr::rename(CPM = summed_CPM_codon)
features <- as_tibble(fread("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/features.csv"))

# features 
# 
# mean_delta_per_transcript <- 
#   CPMs %>% 
#   filter(MD30 == "M") %>% 
#   spread(key = IP, value = CPM) %>% 
#   mutate(delta = CNOT3 - Tot) %>% 
#   select(replicate, transcript, codon, delta) %>% 
#   group_by(transcript, codon) %>% 
#   summarise(mean_delta = mean(delta))
# 
# mean_delta_per_transcript_tagged <- 
#   mean_delta_per_transcript %>% 
#   inner_join(features%>% select(Gene, ENST, signal_sequence),
#              by = c("transcript" = "ENST")) %>% 
#   filter(signal_sequence == T)
# 
# deltas_split <- 
#   split(mean_delta_per_transcript_tagged,
#         mean_delta_per_transcript_tagged$transcript)
# 
# plot_deltas <- function(x){
#   
#   gplot <- 
#     x %>% 
#     filter(codon < 400) %>% 
#     ggplot(aes(x = codon, y = mean_delta)) +
#     ggtitle(unique(x$Gene))+
#     geom_point() +
#     geom_line()
#   
#   return(gplot)
# }
# 
# gplot_deltas <- lapply(deltas_split, plot_deltas)
# 
# gplot_deltas[[50]]

mean_CPMs <- 
  CPMs %>% 
  filter(MD30 == "M") %>% 
  group_by(transcript, codon, IP) %>% 
  summarise(CPM = mean(CPM)) 

mean_CPMs_SPs <- 
  mean_CPMs %>% 
  inner_join(features %>% select(Gene, ENST, signal_sequence),
             by = c("transcript" = "ENST")) %>% 
  filter(signal_sequence == T)

mean_CPMs_SPs <- mean_CPMs_SPs %>% filter(transcript %in% C3_enriched_transcripts)

CPMs_split <- split(mean_CPMs_SPs, mean_CPMs_SPs$transcript)

plot_CPMs <- 
  function(x){
    
    gplot <- 
      x %>% 
      filter(codon >= 1 & codon <= 50) %>% 
      ggplot(aes(x = codon, y = CPM, colour = IP)) +
      ggtitle(unique(x$Gene))+
      scale_color_brewer(palette = "Set1") +
      geom_point() +
      geom_line() +
      publication_theme()
    
    return(gplot)
    
  }

names <- 
  as.vector(unlist(lapply(CPMs_split, function(x){x$Gene %>% unique()})))

#gplots_CPMs <- lapply(CPMs_split, plot_CPMs)

#gplots_CPMs[[13]]

# mean_CPMs%>% 
#   inner_join(features %>% select(Gene, ENST, signal_sequence),
#              by = c("transcript" = "ENST")) %>% 
#   ungroup() %>% 
#   filter(codon <= 50) %>% 
#   spread(key = IP, value = CPM) %>% 
#   mutate(delta = CNOT3 - Tot) %>% 
#   select(codon, signal_sequence, delta) %>% 
#   group_by(signal_sequence, codon) %>% 
#   summarise(mean_delta = mean(delta)) %>% 
#   ggplot(aes(x = codon, y = mean_delta, colour = signal_sequence)) +
#   geom_line() +
#   geom_point()

# setwd("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/CPMs/CNOT3_enriched")
# 
# 
# 
# for(i in 1:length(gplots_CPMs)){
#   png(paste0(names[[i]], "_CPMs_SPs_50_150_CNOT3.png"),
#       res = 300, height = 1000, width = 2000)
#   print(gplots_CPMs[[i]])
#   dev.off()
  
}



stall_sites <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/stall_sites.csv") %>% 
  select(-`...1`)

master <- fread("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/master.csv", drop = "V1", header = T)

stall_sites$classification %>% unique()

stall_sites <- 
  stall_sites %>% filter(classification == "enriched")

stall_sites_1_50 <- 
  stall_sites %>% 
  inner_join(features %>% select(ENST, signal_sequence), by = c("transcript" = "ENST")) %>% 
  filter(signal_sequence == T) %>% 
  filter(codon >= 6) %>% 
  filter(codon <= 50)

stall_sites_1_50 <- stall_sites_1_50 %>% mutate(tag = paste0(transcript, "_", codon))

stalls_split <- 
  split(stall_sites_1_50, stall_sites_1_50$tag)

extract_regions <- 
  function(x){
    
    tran <- x$transcript
    original_codon <- x$codon
    start <- original_codon - 5
    end <- original_codon + 5
    
    df <- 
      data.frame(transcript = tran,
                 stall = original_codon,
                 raw_codons = start:end,
                 relative_positions = -5:5)
    
    #df %>% inner_join(mean_delta_per_transcript, by = c("transcript", "raw_codons" = "codon"))
    
    return(df)
    
    
  }

regions_extracted <- do.call("rbind", lapply(stalls_split, extract_regions))
rownames(regions_extracted) <- NULL

regions_extracted %>% 
  inner_join(mean_delta_per_transcript, by = c("transcript", "raw_codons" = "codon")) %>% 
  group_by(relative_positions) %>% 
  summarise(delta = mean(mean_delta)) %>% 
  ggplot(aes(x = relative_positions, y = delta)) +
  geom_point()

aligned_CPMs <- 
  regions_extracted %>% 
  inner_join(mean_CPMs_SPs %>% 
               spread(key = IP, value = CPM), 
             by = c("transcript", "raw_codons" = "codon")) %>% 
  select(transcript, stall, raw_codons, relative_positions, CNOT3, Tot) %>% 
  gather(key = "IP", value = "CPM", CNOT3:Tot) %>% 
  group_by(relative_positions, IP) %>% 
  summarise(mean_CPM = mean(CPM)) 

CPMs_aligned_stalls <- 
  ggplot(aligned_CPMs, aes(x = relative_positions, y = mean_CPM, colour = IP)) +
  geom_line(size = 2) +
  #geom_point(size = 1) +
  ylab("CPM") +
  ggtitle("Aligned stalls") +
  scale_x_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5), 
                     labels = c("-5", "-4", "-3", "-2", "E", "P", "A", "+2", "+3", "+4", "+5")) +
  scale_color_manual(values = c("#0C7BDC", "#FFC20A")) +
  publication_theme() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank())

regions_extracted %>% 
  inner_join(mean_delta_per_transcript, by = c("transcript", "raw_codons" = "codon")) %>% 
  group_by(relative_positions) %>% 
  summarise(delta = mean(mean_delta)) %>% 
  ggplot(aes(x = relative_positions, y = delta)) +
  geom_point()

# what codon at position 0?

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
my_lyst <- sapply(master$nucleotide_sequence_CDS, split_into_codons)

# Name each element of the list by its ENST tag
transcript_names <- pull(master, var = ENST)

# make a function to attach 
add_positionality <- 
  function(x){
    
    codons <- my_lyst[[x]]
    ENST <- transcript_names[[x]]
    
    df_final <-
      data.frame(transcript = ENST,
                 codon = codons,
                 position = 1:length(codons)) %>% 
      as_tibble()
    
    return(df_final)
    
  }

# raw table of codons per transcript
codons_positions <- 
  do.call("rbind", lapply(1:length(transcript_names), add_positionality))

codons_relative_positions <- 
  regions_extracted %>% 
  inner_join(codons_positions, by = c("transcript", "raw_codons" = "position")) %>% 
  group_by(relative_positions, codon) %>% 
  tally() %>% 
  ungroup()

total_n <- 
  codons_relative_positions %>% 
  group_by(relative_positions) %>% 
  summarise(total = sum(n))

stall_sites_freq %>% 
  filter(percentage > 6) 

stall_sites_freq <- 
  codons_relative_positions %>% 
  inner_join(total_n, by = "relative_positions") %>% 
  mutate(percentage = n/total * 100) %>%
  mutate(codons2colour = case_when(codon == "ctc" & relative_positions %in% c(-4,-1,0, 3) ~ "Leu-CTC",
                                   codon == "ctg" & relative_positions %in% c(-5,-2,-1,1,2,3,4,5) ~ "Leu-CTG",
                                   codon == "gcc" & relative_positions %in% c(-5,-3,-1,0,2,4) ~ "Ala-GCC",
                                   codon == "ccg" & relative_positions %in% c(0) ~ "Pro-CCG",
                                   .default = "other")) 

stall_sites_freq$codons2colour <- factor(stall_sites_freq$codons2colour,
                                         levels = c("Leu-CTG", "Leu-CTC", "Pro-CCG", "Ala-GCC", "other"))

stall_sites_codons <- 
  stall_sites_freq %>% 
  ggplot(aes(x = relative_positions, y = percentage, colour = codons2colour)) +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_color_manual(values = c("#377EB8", "#984EA3", "#A65628", "#E6AB02", "#999999")) +
  scale_x_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5), 
                     labels = c("-5", "-4", "-3", "-2", "E", "P", "A", "+2", "+3", "+4", "+5")) +
  ylab("Codon percentage") +
  publication_theme() +
  #ggtitle("CNOT3 enrcihed sites:\ncodons 50-150") +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank(),
        axis.title.x = element_blank())

#Import codons table
path = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/useful_tables"
codons <- read_csv(paste0(path, "/codon_box_types.csv"))

#Make the capitalised codons lower case
codons$codon <- factor(tolower(codons$codon))

#replace letters "u" with "t"
codons$codon <- gsub("u","t", codons$codon)

#stall_sites_AAs %>% filter(sum_perc > 11) %>% filter(AA == "Ala")

# stall_sites_AAs <- 
#   stall_sites_freq %>%
#   ungroup() %>% 
#   inner_join(codons, by = "codon") %>% 
#   group_by(relative_positions,AA) %>% 
#   summarise(sum_perc = sum(percentage)) %>% 
#   mutate(AA_colour = factor(case_when(AA == "Leu" & relative_positions %in% c(-5,-4,-2,-1,0,1,2,3,4,5) ~"Leu",
#                                      # AA == "Gly" & relative_positions %in% c(-1,1) ~ "Gly",
#                                       AA == "Pro" & relative_positions %in% c(-2, -1,0) ~ "Pro",
#                                       #AA == "Arg" & relative_positions == 0 ~ "Arg",
#                                       AA == "Ala" & relative_positions %in% c(-5,-4,-3, -1,4) ~"Ala",
#                                       .default = "other"),
#                             levels = c("Leu", "Pro", "Ala" "other")))

# AA_gplot <- 
#   stall_sites_AAs %>% 
#   ggplot(aes(x = relative_positions, y = sum_perc, colour = AA_colour)) +
#   geom_point(position = position_jitter(width = 0.2)) +
#   scale_color_manual(values = c("#FF1F5B", "#00CD6C", "#009ADE", "#999999")) +
#   scale_x_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5), 
#                      labels = c("-5", "-4", "-3", "-2", "E", "P", "A", "+2", "+3", "+4", "+5")) +
#   ylab("AA Percentage") +
#   #ggtitle("CNOT3 enrcihed sites:\ncodons 50-150") +
#   publication_theme() +
#   theme(legend.position = "right",
#         legend.direction = "vertical",
#         legend.title = element_blank(),
#         axis.title.x = element_blank())

raw_numbers <- 
  stall_sites_freq %>%
  ungroup() %>% 
  inner_join(codons, by = "codon") %>% 
  select(relative_positions, codon, n, AA, box) %>% 
  mutate(box = as.numeric(str_remove(box, "-box")))

RSCU_data <- 
  raw_numbers %>% 
  group_by(relative_positions, AA) %>% 
  summarise(AA_tot = sum(n)) %>% 
  inner_join(raw_numbers, by = c("relative_positions", "AA")) %>% 
  mutate(RSCU = n/AA_tot * box)

interesting_RSCUs <- 
  RSCU_data %>% 
  filter(RSCU > 3) 

RSCU_gplot <- 
  RSCU_data %>% 
  mutate(RSCU_colour = factor(case_when(codon == "ctg" & relative_positions %in% c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5) ~ "Leu-CTG",
                                        codon == "agc" & relative_positions %in% c(-1, 1, 5) ~ "Ser-AGC",
                                        codon == "atc" & relative_positions %in% c(-5, -4, -3, 0, 1) ~ "Ile-ATC",
                                        codon == "gtg" & relative_positions %in% c(-4, -1, 0) ~ "Val-GTG",
                                        codon == "acc" & relative_positions %in% c(4) ~ "Thr-ACC",
                                        .default = "other"),
                              levels = c("Leu-CTG", "Ile-ATC", "Val-GTG", "Ser-AGC", "Thr-ACC", "other"))) %>% 
  ggplot(aes(x = relative_positions, y = RSCU, colour = RSCU_colour)) +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_color_manual(values = c("#377EB8", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#999999")) +
  scale_x_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5), 
                     labels = c("-5", "-4", "-3", "-2", "E", "P", "A", "+2", "+3", "+4", "+5")) +
  publication_theme() +
  theme(legend.position = "right", legend.direction = "vertical", legend.title = element_blank(), axis.title.x = element_blank())


png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/SP_1_50_stallsitesaligned_codonfreq.png",
    res = 300, height = 1000, width = 1500)
print(stall_sites_codons)
dev.off()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/SP_1_50_stallsitesaligned_AAfreq.png",
    res = 300, height = 1000, width = 1500)
print(AA_gplot)
dev.off()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/SP50150_stallsitesaligned_RSCU.png",
    res = 300, height = 1000, width = 1500)
print(RSCU_gplot)
dev.off()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/SP_1_50_stallsitesaligned_CPMs.png",
    res = 300, height = 1000, width = 1500)
print(CPMs_aligned_stalls)
dev.off()

CNOT3_SP_transcripts <- mean_CPMs_SPs %>% pull(transcript) %>% unique()

AA_sequences <- 
  master %>% 
  filter(ENST %in% CNOT3_SP_transcripts) %>% 
  as_tibble() %>% 
  select(ENST, amino_acid_sequence) %>% 
  filter(nchar(amino_acid_sequence) < 10000)

seqs <- str_split(AA_sequences$amino_acid_sequence, "")
names_for_fasta <- AA_sequences$ENST

library(seqinr)

write.fasta(sequences = seqs, names = names_for_fasta, file.out = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_SPs_AAs.fasta")

borders <- 
  read_tsv("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_SPs_regions.tsv", col_names = F) %>% 
  select(X1 : X5)

colnames(borders) <- c("ENST", "SigP", "region", "start", "end")

borders_end <- 
  borders %>% filter(region == "c-region") %>% 
  select(ENST, end)

table(borders_end$end <= 14)

centre_data <- 
  function(x){
    
    transcript <- x$ENST
    SP_end <- x$end
    start <- SP_end - 10
    end <- SP_end + 500
    
    df <- 
      data.frame(transcript = transcript,
                 SP_end = SP_end,
                 raw_codon = start:end,
                 rel_pos = -10:500)
    
    return(df)
    
  }

data_centred <- do.call("rbind", lapply(split(borders_end, borders_end$ENST), centre_data))
rownames(data_centred) <- NULL

aligned_ends <- 
  inner_join(data_centred, codons_positions, by = c("transcript", "raw_codon" = "position")) %>% 
  inner_join(mean_CPMs_SPs %>% 
               spread(key = IP, value = CPM), 
             by = c("transcript", "raw_codon" = "codon")) %>% 
  select(-Gene, -signal_sequence)

SP_aligned_ends_CPMs <- 
  aligned_ends %>% 
  gather(key = "IP", value = "CPM", CNOT3:Tot) %>% 
  group_by(IP, rel_pos) %>% 
  summarise(mean_CPM = mean(CPM)) %>% 
  ggplot(aes(x = rel_pos, y = mean_CPM, colour = IP)) +
  geom_line(size = 1) +
  ylab("CPM") +
  xlab("Codon") +
  scale_color_manual(values = c("#0C7BDC", "#FFC20A")) +
  publication_theme() 

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/SPendsaligned_CPMs.png",
    res = 300, height = 1000, width = 2500)
print(SP_aligned_ends_CPMs)
dev.off()


deltas_CPMS <- 
  aligned_ends %>% 
  mutate(delta = CNOT3 - Tot) %>% 
  group_by(rel_pos) %>% 
  summarise(mean_del = mean(delta)) %>% 
  ggplot(aes(x = rel_pos, y = mean_del)) +
  geom_line(colour = ) +
  geom_vline(xintercept  = 41, colour = "#666666", lty = "dashed") +
  scale_x_continuous(breaks = c(0, 41,100,200,300,400,500)) +
  xlab("Codon") +
  ylab("Delta CPM (CNOT3-Tot)") +
  publication_theme()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/SPendsaligned_deltaCPMs.png",
    res = 300, height = 1000, width = 2500)
print(deltas_CPMS)
dev.off()

new_split <- split(aligned_ends, aligned_ends$transcript)

new_split <- aligned_ends %>% 
  filter(rel_pos ==41)

append_positions <- 
  function(x){
    
    min_row <- 41 - 40
    max_row <- 41 + 40
    
    y <- 
      x %>% 
      dplyr::filter(rel_pos >= min_row & rel_pos <= max_row) %>% 
      mutate(new_rel_pos = -40:40)
    
    
    #y %>%  mutate(new_rel_pos = 1:max_row)
    
    return(y)
    
    
  }

#lapply(new_split, append_positions)

exit_site_centred <- do.call("rbind", lapply(new_split, append_positions))
rownames(exit_site_centred) <- NULL

exit_site_centred %>% 
  group_by(rel_pos, codon) %>% tally() %>% 
  ungroup() %>% 
  group_by(rel_pos) %>% 
  summarise(sum_n = sum(n)) %>% 
  ggplot(aes(x = rel_pos, y = sum_n)) +
  geom_line()

CSCs <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/CSCs.csv") %>% 
  mutate(codon = gsub("u", "t", tolower(codon)))

CSCs_secretory <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/CSC_scores_cytosolic_secretory.csv") %>% 
  mutate(codon = gsub("u", "t", tolower(codon))) %>% 
  filter(transcript_type == "secretory") %>% 
  select(-`...1`)

exit_site_centred_CSCs <- 
  exit_site_centred %>% 
  inner_join(CSCs %>% select(codon, ctrl_cor_estimate, delta_cor_estimate), by = "codon")

exit_site_centred_CSCs %>% 
  group_by(new_rel_pos) %>% 
  summarise(mean_CSC = mean(ctrl_cor_estimate)) %>% 
  ggplot(aes(x = new_rel_pos, y = mean_CSC)) +
  geom_line() +
  geom_smooth() +
  geom_point()

exit_site_centred_CSCs <- 
  exit_site_centred %>% 
  inner_join(CSCs_secretory %>% select(codon, ctrl_cor_estimate, delta_cor_estimate), by = "codon")

codon41centred <- 
  exit_site_centred_CSCs %>% 
  group_by(new_rel_pos) %>% 
  summarise(mean_CSC = mean(ctrl_cor_estimate)) %>% 
  ggplot(aes(x = new_rel_pos, y = mean_CSC)) +
  geom_line() +
  geom_point() +
  geom_smooth() +
  xlab("Codon") +
  ylab("Mean CSC") +
  publication_theme()

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/secretoryCSCs_codon41.png",
    res = 300, height = 750, width = 1250)
print(codon41centred)
dev.off()

enumerated_codons <- 
  exit_site_centred_CSCs %>% 
  group_by(new_rel_pos, codon) %>% 
  tally() %>% 
  ungroup()

codon_freqs <- 
  enumerated_codons %>% 
  group_by(new_rel_pos) %>% 
  summarise(tot_n = sum(n)) %>% 
  inner_join(enumerated_codons, by = "new_rel_pos") %>% 
  mutate(codon_perc = n/tot_n * 100) %>% 
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A", "U", "G","C"))) 

AUGC <- 
  codon_freqs %>% 
  group_by(new_rel_pos, wobble) %>% 
  summarise(sum_perc = sum(codon_perc)) %>% 
  ggplot(aes(x = new_rel_pos, y = sum_perc, colour = wobble)) +
  geom_line(size = 1) +
  geom_point() +
  scale_color_manual(values = c("#8EA8E6", "#A3E68E", "#E68EDE", "#E79F36")) +
  scale_y_continuous(limits = c(0,NA)) +
  xlab("Codon") +
  ylab("Percentage") +
  publication_theme()

codon_freqs %>% 
  group_by(new_rel_pos, wobble) %>% 
  summarise(sum_perc = sum(codon_perc)) %>% 
  ggplot(aes(x = new_rel_pos, y = sum_perc, colour = wobble)) +
  geom_line(size = 1) +
  geom_point() +
  scale_color_manual(values = c("#8EA8E6", "#A3E68E", "#E68EDE", "#E79F36")) +
  scale_y_continuous(limits = c(0,NA)) +
  xlab("Codon") +
  ylab("Percentage") +
  coord_cartesian(xlim = c(-20, 20)) +
  publication_theme()

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/AUGC_codon41.png",
    res = 300, height = 1000, width = 1500)
print(AUGC)
dev.off()

codon_freqs %>% 
  filter(new_rel_pos >= -5 & new_rel_pos <= 5) %>% 
  ggplot(aes(x = new_rel_pos, y = codon_perc, group = codon)) +
  geom_line() +
  geom_point() +
  facet_wrap(~codon)

AA_percentages_codon41 <- 
  codon_freqs %>% 
  filter(new_rel_pos >= -40 & new_rel_pos <= 40) %>% 
  inner_join(codons, by = "codon") %>% 
  group_by(new_rel_pos, AA) %>% 
  summarise(sum_perc = sum(codon_perc)) %>% 
  ggplot(aes(x = new_rel_pos, y = sum_perc, group = AA)) +
  geom_line() +
  facet_wrap(~AA, nrow = 9, ncol = 2) +
  ylab("AA Percentage") +
  xlab("Residues position") +
  publication_theme()

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/AA_percentages_codon41.png",
    res = 300, height = 3000, width = 1250)
print(AA_percentages_codon41)
dev.off()


codons <- 
  codons %>% mutate(property = case_when(AA %in% c("Ala", "Val", "Ile", "Leu", "Met", "Phe", "Tyr", "Trp") ~ "hydrophobic",
                                         AA %in% c("Ser", "Thr", "Asn", "Gln") ~ "polar",
                                         AA %in% c("Cys", "Gly", "Pro") ~"special",
                                         AA %in% c("Arg", "His", "Lys") ~ "positive",
                                         AA %in% c("Asp", "Glu") ~ "negative"))

AA_properites <- 
  codon_freqs %>% 
  filter(new_rel_pos >= -40 & new_rel_pos <= 40) %>% 
  inner_join(codons, by = "codon") %>% 
  group_by(new_rel_pos, property) %>% 
  summarise(sum_perc = sum(codon_perc)) %>% 
  ggplot(aes(x = new_rel_pos, y = sum_perc, group = property, colour = property)) +
  geom_line(size = 1) +
  scale_y_continuous(limits = c(0,NA)) +
  ylab("AA Percentage") +
  xlab("Residue position") +
  scale_color_manual(values = c("#7FC97F", "#BEAED4", "#FDC086", "#386CB0", "#BF5B17"))+
  publication_theme() +
  theme(legend.title = element_blank())


png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/AA_properties_codon41.png",
    res = 300, height = 1000, width = 1750)
print(AA_properites)
dev.off()

codon_freqs %>% 
  filter(new_rel_pos >= -40 & new_rel_pos <= 40) %>% 
  inner_join(codons, by = "codon") %>% 
  filter(wobble == "U") %>% 
  ggplot(aes(x = new_rel_pos, y = codon_perc, group = codon)) +
  geom_line() +
  facet_wrap(~codon) +
  scale_y_continuous(limits = c(0,NA))

prolines <- 
  codon_freqs %>% 
  filter(new_rel_pos >= -40 & new_rel_pos <= 40) %>% 
  inner_join(codons, by = "codon") %>% 
  filter(AA == "Pro") %>% 
  ggplot(aes(x = new_rel_pos, y = codon_perc, colour = codon)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#1B9E77", "#BF5B17", "#984EA3")) +
  xlab("Codon position") +
  ylab("Percentage") +
  publication_theme()

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/Prolines_codon41.png",
    res = 300, height = 1000, width = 1250)
print(prolines)
dev.off()

leucines <- 
  codon_freqs %>% 
  filter(new_rel_pos >= -40 & new_rel_pos <= 40) %>% 
  inner_join(codons, by = "codon") %>% 
  filter(AA == "Leu") %>% 
  ggplot(aes(x = new_rel_pos, y = codon_perc, colour = codon)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#1B9E77", "#BF5B17", "#984EA3", "#FF7F00", "#88CCEE")) +
  xlab("Codon position") +
  ylab("Percentage") +
  publication_theme()

Gln <- 
  codon_freqs %>% 
  filter(new_rel_pos >= -40 & new_rel_pos <= 40) %>% 
  inner_join(codons, by = "codon") %>% 
  filter(AA == "Gln") %>% 
  ggplot(aes(x = new_rel_pos, y = codon_perc, colour = codon)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#1B9E77", "#BF5B17", "#984EA3", "#FF7F00", "#88CCEE")) +
  xlab("Codon position") +
  ylab("Percentage") +
  publication_theme()

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/Gln_codon41.png",
    res = 300, height = 1000, width = 1250)
print(Gln)
dev.off()

Ala <- 
  codon_freqs %>% 
  filter(new_rel_pos >= -40 & new_rel_pos <= 40) %>% 
  inner_join(codons, by = "codon") %>% 
  filter(AA == "Ala") %>% 
  ggplot(aes(x = new_rel_pos, y = codon_perc, colour = codon)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#1B9E77", "#BF5B17", "#984EA3", "#FF7F00", "#88CCEE")) +
  xlab("Codon position") +
  ylab("Percentage") +
  publication_theme()

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/Ala_codon41.png",
    res = 300, height = 1000, width = 1250)
print(Ala)
dev.off()

Leu <- 
  codon_freqs %>% 
  filter(new_rel_pos >= -40 & new_rel_pos <= 40) %>% 
  inner_join(codons, by = "codon") %>% 
  filter(AA == "Leu") %>% 
  ggplot(aes(x = new_rel_pos, y = codon_perc, colour = codon)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#1B9E77", "#BF5B17", "#984EA3", "#FF7F00", "#88CCEE")) +
  xlab("Codon position") +
  ylab("Percentage") +
  publication_theme()

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/Leu_codon41.png",
    res = 300, height = 1000, width = 1250)
print(Leu)
dev.off()

Glu <- 
  codon_freqs %>% 
  filter(new_rel_pos >= -40 & new_rel_pos <= 40) %>% 
  inner_join(codons, by = "codon") %>% 
  filter(AA == "Glu") %>% 
  ggplot(aes(x = new_rel_pos, y = codon_perc, colour = codon)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#1B9E77", "#BF5B17", "#984EA3", "#FF7F00", "#88CCEE")) +
  xlab("Codon position") +
  ylab("Percentage") +
  publication_theme()

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/Glu_codon41.png",
    res = 300, height = 1000, width = 1250)
print(Glu)
dev.off()

Thr <- 
  codon_freqs %>% 
  filter(new_rel_pos >= -40 & new_rel_pos <= 40) %>% 
  inner_join(codons, by = "codon") %>% 
  filter(AA == "Thr") %>% 
  ggplot(aes(x = new_rel_pos, y = codon_perc, colour = codon)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#1B9E77", "#BF5B17", "#984EA3", "#FF7F00", "#88CCEE")) +
  xlab("Codon position") +
  ylab("Percentage") +
  publication_theme()

Gly <- 
  codon_freqs %>% 
  filter(new_rel_pos >= -40 & new_rel_pos <= 40) %>% 
  inner_join(codons, by = "codon") %>% 
  filter(AA == "Gly") %>% 
  ggplot(aes(x = new_rel_pos, y = codon_perc, colour = codon)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#1B9E77", "#BF5B17", "#984EA3", "#FF7F00", "#88CCEE")) +
  xlab("Codon position") +
  ylab("Percentage") +
  publication_theme()

Phe <- 
  codon_freqs %>% 
  filter(new_rel_pos >= -40 & new_rel_pos <= 40) %>% 
  inner_join(codons, by = "codon") %>% 
  filter(AA == "Phe") %>% 
  ggplot(aes(x = new_rel_pos, y = codon_perc, colour = codon)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#377EB8", "#1B9E77", "#BF5B17", "#984EA3", "#FF7F00", "#88CCEE")) +
  xlab("Codon position") +
  ylab("Percentage") +
  publication_theme()

enumerated_codons %>% filter()

RColorBrewer::brewer.pal(10, "Accent")
"#377EB8", "#984EA3","#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F", "#999999" 
"#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E"  "#A6761D" "