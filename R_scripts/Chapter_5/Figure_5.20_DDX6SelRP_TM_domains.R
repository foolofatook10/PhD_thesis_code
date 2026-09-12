library(tidyverse)
library(data.table)
library(smplot2)

#Load theme
publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = 28, hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.title.x = element_text(size = 16, vjust = -0.2),
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
            plot.margin=unit(c(5,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

# Load in CPMs per dataset
CPMS <- as_tibble(fread("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/riboseq_datasets/DDX6_CDS_CPMs.csv"))

# Calculate the mean CPM per codon by dividing summed CPM count by 3
#CPMS$summed_CPM_codon <- (CPMS$summed_CPM_codon/3) 

# Filter out for D30s and tidy headings
CPMS <-
  CPMS %>%
  #filter(MD30!= "D30") %>%
  select(replicate, transcript, codon, CPM, IP = condition)

# calculate the mean CPM at each codon position within each transcript
mean_CPMs <-
  CPMS %>%
  group_by(transcript, codon, IP) %>%
  summarise(
    #n = n(),
    mean_CPM = mean(CPM))


# As eventually merging this data to TM data, CPM data is split so that inner_join doesn't freak out (each transcript is represented twice, one for CNOT3 IP'd RPFs and the other for translatome RPFs) 
mean_CPMs_split <- split(mean_CPMs, mean_CPMs$IP)

# Load in preprocessed TM data

##ENST = Ensembl ID
## codon_position = the actual codon position within the CDS of the transcript
## codon = the codon ID of the codon position
## class = the class of transcript defined by the number of signal peptide (SP) or transmembrane (TM) domains it has.I have called a targeting sequence (TS) to mean either an SP or TM domain
## region = each codon position is tagged with whether it belongs to the region before the TM domain (prefix), the TS itself (TS) or the region after the TS (suffix)
## ID = the TS number in order of appearance. ID = 0 is an SP whilst ID > 0 is by definition a TM domain. The prefixes and suffixes of a TS are also ID tagged
## rel_position = the relative codon position approaching the start of the TS domain (region = prefix), from the start of the TS (region = TS) or from the end of the TM domain (suffix)
TM_domains <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/TM_domains_FINAL.csv") %>% 
  #mutate(tag = paste0(ENST, "_", ID)) %>% 
  select(-`...1`) %>% 
  filter(class != "sponly") %>% # Not interested in mRNAs with an N terminal SP and no TM domains
  filter(ID != 0) %>% # ID = 0 is an SP (which can occur in single pass type 1 transcripts and multipass type 1 transcripts). This analysis purely focuses on TM domains within the body of the CDS
  mutate(tag = paste0(ENST, "_", ID)) # Create a unique tag for each TM domain along each transcript

TM_domains$region <- factor(TM_domains$region, levels = c("prefix", "TS", "suffix")) # TM domains have an associated prefix (region before TM domain) and suffix (region after TM domain). None of these regions should overlap

# Problem 1: The relative_position of the TM domains is from the start of the TM domain. Here I wish to align the TM domains by their ends and only plot the suffix region afterwards, all on one graph. Therefore I need to reverse the relative positions such that they count down:
new_data <- TM_domains %>% filter(region == "TS") # filter for just TM domains
new_data_split <-  split(new_data, new_data$tag) # Split this data into unique TM domains across the transcriptome


# reverse relative position
TMs_ammended <- 
  do.call("rbind",lapply(new_data_split,
                         function(x){
                           
                           relstart <- ((-1* max(x$rel_position))) # get the max length of the TS and * -1 to count up from
                           relend <- -1
                           
                           return(x %>% 
                                    mutate(rel_position = relstart:relend))
                           
                           
                         }))

# Purely for graphical purposes, I want the TS to end at position 0 rather than -1 (as the suffix relative positions begin at 1 so it looks nice if it is continuous). THerefore I plus one to all relative positions.
TMs_ammended$rel_position <- TMs_ammended$rel_position + 1 

# recombine TMs with correct rel_positions to the suffixes
corrected_TMs <- 
  rbind(TMs_ammended,
        TM_domains %>% filter(region == "suffix"))

# remember, after this point this is just TM domains within the CNOT3 data.
all_data <-
  do.call("rbind", lapply(mean_CPMs_split, function(x){
    
    return(x %>% 
             inner_join(corrected_TMs %>% 
                          select(ENST, region, rel_position, codon_sequence = codon, codon_position, ID, class, tag),
                        by = c("transcript" = "ENST",
                               "codon" = "codon_position")))
    
    
  })) %>%
  select(transcript, codon, IP, mean_CPM, region, rel_position, codon_sequence, ID, class, tag)

# just consider -21 codons from the TM ends (vast majority of TM domains are ~21 codons in length -- see below)
all_data_filt <- 
  all_data %>% 
  filter(rel_position >= -21) #%>%

# check to see if any transcripts contain codon positions that are outrageously high. If so, just remove the entire transcript
very_high_CPMs <- 
  all_data_filt %>% 
  filter(mean_CPM > 5 | mean_CPM < -5) %>% 
  pull(transcript) %>% 
  unique() #Identifies 70 transcripts

# Exclude these transcripts
all_data <- 
  all_data %>% 
  filter(!transcript %in% very_high_CPMs)


# Number of unique transcripts (one transcript per gene) for each of the classes represented in the DDX6 SelRP dataset
all_data %>% ungroup() %>%  select(transcript, class) %>% unique() %>% pull(class) %>% table()

# Number of TM domains considered 
all_data %>% 
  ungroup() %>%  
  select(class, tag) %>% 
  unique() %>% 
  group_by(class) %>% 
  summarise(n_TMs = n())

# Understand the length distribution of these TM domains
## Just filter for the TM domains and one version of them (i.e the CNOT3 condition)
DDX6_TM_data <- 
  all_data %>% 
  filter(region == "TS") %>% 
  filter(IP == "DDX6")

# Split the data
DDX6_TM_data_split <- 
  split(DDX6_TM_data, DDX6_TM_data$tag)

# understand the lengths of all TM domains
TM_lengths_gplot <- 
  unlist(lapply(DDX6_TM_data_split, function(x){(min(x$rel_position) -1) * -1})) %>% 
  table() %>% 
  as.data.frame() %>% 
  mutate(percentage = Freq/sum(Freq) * 100) %>% 
  dplyr::rename(length_of_TM = ".") %>% 
  ggplot(aes(x = length_of_TM, y = percentage, label = Freq)) +
  geom_col() +
  geom_text(colour = "red", vjust = 0) +
  xlab("Length of TM domains (codons)") +
  ylab("Proportion of TM domains (%)") +
  publication_theme()

thesis_figures_path = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures"

png(file.path(thesis_figures_path, "TMs_in_DDX6_SelRP_TM_lengths.png"),
    res = 300, height = 1250, width = 2000)
print(TM_lengths_gplot)
dev.off()

# here is where you filter for which transcripts to consider for the resulting gplots
multipass_transcript_names <- c("multipassT1", "multipassT2", "multipassT3")
singlepass_transcript_names <- c("type1", "type2", "type3")

transcripts_considered <- all_data %>% filter(class == "type3")

# Calculate mean CPM
test4 <- transcripts_considered %>% 
  dplyr::rename(CPM = mean_CPM) %>% 
  group_by(region, rel_position, IP) %>%
  summarise(mean_CPM = mean(CPM),
            n = n(),
            sd_CPM = sd(CPM, na.rm = T),
            SE_CPM = sd_CPM/sqrt(n))  


actual_plot <- 
  test4 %>% 
  filter(rel_position >=-21) %>% 
  ggplot(aes(x = rel_position, y = mean_CPM, colour = IP, group = IP))   +
  geom_rect(aes(xmin=-21, xmax=0, ymin=-Inf, ymax=Inf),
            fill = alpha("#DECBE4", 0.05),
            colour = alpha("#DECBE4", 0.05),
            inherit.aes = F)+
  geom_rect(aes(xmin=0, xmax=100, ymin=-Inf, ymax=Inf),
            fill = alpha("#B3CDE3", 0.05),
            colour = alpha("#B3CDE3", 0.05),
            inherit.aes = F)+
  scale_color_manual(values = c("#00A087B2", "#FFC20A")) +
  scale_fill_manual(
    values = c("#00A087B2", "#FFC20A")
  ) +
  geom_vline(xintercept = 0, lty = "longdash", size = 1, colour = "black") +
  #geom_ribbon(aes(ymin = mean_CPM- mean_SE, ymax = mean_CPM+ mean_SE, fill = IP), alpha = 0.3) +
  geom_line(size = 1)  +
  geom_ribbon(aes(ymin = mean_CPM- SE_CPM, 
                  ymax = mean_CPM + SE_CPM, 
                  fill = IP), alpha = 0.4, colour = NA) +
  
  scale_x_continuous(limits = c(-21, 100),
                     breaks = seq(-20, 100, 20),
                     labels = c("-20", "0", "+20", "+40", "+60", "+80", "+100"))  +
  scale_y_continuous(limits = c(0,NA)) +
  xlab("Relative codon position to TM end") +
  ylab("Mean CPM\n(DDX6 IP'd - Total RPFs)") +
  ggtitle("Type 3\nTM domains") +
  publication_theme() +
  theme(legend.position = "none")

confidence_plot <- 
  test4 %>% 
  filter(rel_position >=-21) %>% 
  filter(IP == "DDX6") %>% 
  ggplot(aes(x = rel_position, y = n)) +
  geom_rect(aes(xmin=-21, xmax=0, ymin=-Inf, ymax=Inf),
            fill = alpha("#DECBE4", 0.05),
            colour = alpha("#DECBE4", 0.05),
            inherit.aes = F)+
  geom_rect(aes(xmin=0, xmax=100, ymin=-Inf, ymax=Inf),
            fill = alpha("#B3CDE3", 0.05),
            colour = alpha("#B3CDE3", 0.05),
            inherit.aes = F) +
  geom_line(size = 1, colour = "black") +
  scale_y_continuous(limits = c(0,NA), breaks = c(seq(0,7000,25))) +
  scale_x_continuous(limits = c(-21, 100),
                     breaks = seq(-20, 100, 20),
                     labels = c("-20", "0", "+20", "+40", "+60", "+80", "+100")) +
  geom_vline(xintercept = 0, lty = "longdash", size = 1, colour = "black") +
  xlab("Relative codon position to TM end") +
  ylab("Number of codons") +
  publication_theme()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/DDX6_SelRP_singlepassT3_TM_domains.png", 
    res = 300, 
    height = 2250, 
    width = 1500)
print(cowplot::plot_grid(actual_plot, confidence_plot, ncol = 1, align = "v", rel_heights = c(1,0.5)))
dev.off()



generate_aucs <- function(filter_term){
  
  if(filter_term== "all_TMs"){
    
    stats_analysis <- 
      all_data %>% 
      dplyr::rename(CPM = mean_CPM) %>%
      ungroup() %>% 
      select(region, IP, CPM, rel_position, tag)%>% 
      group_by(region, tag, IP) %>% 
      summarise(auc = sm_auc(rel_position, CPM)) %>% 
      mutate(subgroup = filter_term)
    
  } else if(filter_term == "single_pass"){
    
    stats_analysis <- 
      all_data %>% 
      filter(class %in% c("type1", "type2", "type3")) %>% 
      dplyr::rename(CPM = mean_CPM) %>%
      ungroup() %>% 
      select(region, IP, CPM, rel_position, tag)%>% 
      group_by(region, tag, IP) %>% 
      summarise(auc = sm_auc(rel_position, CPM)) %>% 
      mutate(subgroup = filter_term)
    
  } else if(filter_term == "multi_pass"){
    
    stats_analysis <- 
      all_data %>% 
      filter(class %in% c("multipassT1", "multipassT2", "multipassT3")) %>% 
      dplyr::rename(CPM = mean_CPM) %>%
      ungroup() %>% 
      select(region, IP, CPM, rel_position, tag)%>% 
      group_by(region, tag, IP) %>% 
      summarise(auc = sm_auc(rel_position, CPM)) %>% 
      mutate(subgroup = filter_term)
    
  } else{
    
    stats_analysis <- 
      all_data %>% 
      filter(class == filter_term) %>% 
      dplyr::rename(CPM = mean_CPM) %>%
      ungroup() %>% 
      select(region, IP, CPM, rel_position, tag)%>% 
      group_by(region, tag, IP) %>% 
      summarise(auc = sm_auc(rel_position, CPM)) %>% 
      mutate(subgroup = filter_term)
    
    
  }
  
  return(stats_analysis)
  
}

type1_aucs <- generate_aucs("type1")
type2_aucs <- generate_aucs("type2")
type3_aucs <- generate_aucs("type3")

multipass_type1_aucs <- generate_aucs("multipassT1")
multipass_type2_aucs <- generate_aucs("multipassT2")
multipass_type3_aucs <- generate_aucs("multipassT3")

all_TMs_aucs <- generate_aucs("all_TMs")
single_pass_aucs <- generate_aucs("single_pass")
multi_pass_aucs <- generate_aucs("multi_pass")

# multipass_type3_aucs %>% 
#   ggplot(aes(x = IP, y = log10(auc))) +
#   geom_violin() +
#   geom_boxplot(width = 0.2) +
#   facet_wrap(~region, scales = "free_y")

all_stats_data <- 
  do.call("rbind", list(#type1_aucs, 
    #type2_aucs, 
    #type3_aucs, 
    #multipass_type1_aucs, 
    multipass_type2_aucs, 
    #multipass_type3_aucs,
    all_TMs_aucs,
    single_pass_aucs,
    multi_pass_aucs))

stats_data_final <- 
  all_stats_data %>%
  #mutate(auc = log10(auc)) %>% 
  select(region, subgroup, tag, IP, auc) %>%
  pivot_wider(
    names_from = IP,
    values_from = auc
  ) %>%
  group_by(subgroup, region) %>%
  group_modify(~ {
    test <- t.test(.x$DDX6, .x$Total, paired = TRUE)
    
    tibble(
      #n = sum(complete.cases(.x$CNOT3, .x$Tot)),
      #mean_CNOT3 = mean(.x$CNOT3, na.rm = TRUE),
      #mean_Tot = mean(.x$Tot, na.rm = TRUE),
      #mean_difference = mean(.x$CNOT3 - .x$Tot, na.rm = TRUE),
      #t = unname(test$statistic),
      #df = unname(test$parameter),
      p = test$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    significant = p_adj < 0.05
  ) %>%
  arrange(subgroup) %>% 
  mutate(stars = case_when(p_adj < 0.05 & p_adj >= 0.01 ~ "*",
                           p_adj < 0.01 & p_adj >= 0.001 ~ "**",
                           p_adj < 0.001 ~ "***",
                           .default = "N.S"))

thesis_figures_path = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures"

write.csv(stats_data_final, file.path(thesis_figures_path, "DDX6SelRP_TM_domain_stats_all.csv"))
