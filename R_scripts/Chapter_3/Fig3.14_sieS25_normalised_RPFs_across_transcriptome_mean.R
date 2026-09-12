#load packages----
library(tidyverse)
library(grid)
library(gridExtra)
#library(viridis)

#read in common variables----
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/R_scripts/common_variables.R")

#set what you have called your control and treated samples. This can be a vector of strings if more than one treatment has been used.
control <- "NTC"
treatment <- "siRPS25"

#read in functions----
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/R_scripts/meta_plots/binning_RiboSeq_functions.R")

#create themes----
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
            axis.text.x = element_text(size = 18, color = "black"),
            axis.text.y = element_text(size = 18, color = "black"),
            #axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = 16),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(1.2, "cm"),
            legend.margin = unit(0.5, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

my_theme <- publication_theme()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_blank())

UTR5_theme <- my_theme+
  theme(legend.position="none",
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank())

CDS_theme <- my_theme+
  theme(legend.position="none",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line.y = element_blank())

UTR3_theme <- my_theme+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.line.y = element_blank())

plot_binned_lines <- 
  function(df, SD = F, conditions = NULL, mylabels = NULL, colours = NULL, CDS_only = F) {
    
    
    
    #CDS
    df %>%
      ggplot(aes(x = grouping, y = average_counts, colour = condition))+
      geom_line(size = 1)+
      geom_point() +
      facet_wrap(~class, scales = "free") +
      {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
      {if(!is.null(colours))scale_colour_manual(values=colours)}+
      {if(!is.null(colours))scale_fill_manual(values=colours)}+
      ggtitle("CDS: Enriched") +
      CDS_theme()-> CDS_plot
    
    
    return(CDS_plot)
    
  }




#read in data----
load(file = file.path(parent_dir, "Counts_files/R_objects/binned_list_siRPS25.Rdata"))
summary(binned_list[[1]])
head(binned_list[[1]])



region_lengths <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_region_lengths.csv", col_names = c("transcript", "UTR5_len", "CDS_len", "UTR3_len"))
most_abundant_transcripts <- read_csv(file = file.path("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))


# read in merged data
data <- read_csv(file.path("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/Analysis/DESeq2_output/merged_DESeq2.csv")) %>% 
  select(gene_sym, gene, transcript, RPFs_group)

data$RPFs_group <- factor(data$RPFs_group)
summary(data$RPFs_group)

RPFsup_transcripts <-  data %>% filter(RPFs_group == "RPFs up") %>% pull(transcript)
RPFsdown_transcripts <-  data %>% filter(RPFs_group == "RPFs down") %>% pull(transcript)
Totalsup_transcripts <-  data %>% filter(RPFs_group == "Totals up") %>% pull(transcript)
Totalsdown_transcripts <-  data %>% filter(RPFs_group == "Totals down") %>% pull(transcript)
bothup_transcripts <-  data %>% filter(RPFs_group == "both up") %>% pull(transcript)
bothdown_transcripts <-  data %>% filter(RPFs_group == "both down") %>% pull(transcript)


binned_list_RPFs_up <-  
  lapply(binned_list, function(x) {x %>% filter(transcript %in% RPFsup_transcripts)})

binned_list_RPFs_down <-  
  lapply(binned_list, function(x) {x %>% filter(transcript %in% RPFsdown_transcripts)})

binned_list_totals_up <-  
  lapply(binned_list, function(x) {x %>% filter(transcript %in% Totalsup_transcripts)})

binned_list_totals_down <-  
  lapply(binned_list, function(x) {x %>% filter(transcript %in% Totalsdown_transcripts)})

binned_list_both_up <-  
  lapply(binned_list, function(x) {x %>% filter(transcript %in% bothup_transcripts)})

binned_list_both_down <-  
  lapply(binned_list, function(x) {x %>% filter(transcript %in% bothdown_transcripts)})

unlist(lapply(binned_list_RPFs_up, function(x){x %>% pull(transcript) %>% unique() %>% length()})) %>% unique()
unlist(lapply(binned_list_RPFs_down, function(x){x %>% pull(transcript) %>% unique() %>% length()})) %>% unique()
unlist(lapply(binned_list_totals_up, function(x){x %>% pull(transcript) %>% unique() %>% length()})) %>% unique()
unlist(lapply(binned_list_totals_down, function(x){x %>% pull(transcript) %>% unique() %>% length()})) %>% unique()
unlist(lapply(binned_list_both_up, function(x){x %>% pull(transcript) %>% unique() %>% length()})) %>% unique()
unlist(lapply(binned_list_both_down, function(x){x %>% pull(transcript) %>% unique() %>% length()})) %>% unique()

summary(data$RPFs_group)




#all transcripts----
#summarise
#summarise within each sample
summarised_binned_list_RPFs_up <- lapply(binned_list_RPFs_up, summarise_data, value = "binned_cpm", grouping = "bin")
summarised_binned_list_RPFs_down <- lapply(binned_list_RPFs_down, summarise_data, value = "binned_cpm", grouping = "bin")

summarised_binned_list_totals_up <- lapply(binned_list_totals_up, summarise_data, value = "binned_cpm", grouping = "bin")
summarised_binned_list_totals_down <- lapply(binned_list_totals_down, summarise_data, value = "binned_cpm", grouping = "bin")

summarised_binned_list_both_up <- lapply(binned_list_both_up, summarise_data, value = "binned_cpm", grouping = "bin")
summarised_binned_list_both_down <- lapply(binned_list_both_down, summarise_data, value = "binned_cpm", grouping = "bin")

summarised_binned_list_all <- lapply(binned_list, summarise_data, value = "binned_cpm", grouping = "bin")


#summarise within each condition (across replicates)
do.call("rbind", summarised_binned_list_RPFs_up) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned_RPFs_up

do.call("rbind", summarised_binned_list_RPFs_down) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned_RPFs_down

do.call("rbind", summarised_binned_list_totals_up) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned_totals_up

do.call("rbind", summarised_binned_list_totals_down) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned_totals_down

do.call("rbind", summarised_binned_list_both_up) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned_both_up

do.call("rbind", summarised_binned_list_both_down) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned_both_down

do.call("rbind", summarised_binned_list_all) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned_all

summarised_binned_all$condition <- factor(summarised_binned_all$condition,
                                          levels = c("NTC", "siRPS25"),
                                          labels = c("NTC", "sieS25"))

summarised_binned_all %>%
  filter(region == "CDS") %>% 
  mutate(grouping = grouping * 2) %>% 
  ggplot(aes(x = grouping, y = average_counts, colour = condition))+
  geom_line(size = 1)+
  #geom_smooth() +
  geom_point() +
  scale_color_manual(values = c("#FFC20A", "#984EA3")) +
  #facet_wrap(~class) +
  scale_y_continuous(limits = c(0,NA)) +
  #ggtitle("All transcripts") +
  xlab("CDS%") +
  ylab("Mean nCPMs") +
  publication_theme() +
  theme(legend.position = "none")-> gplot_all

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/sieS25_CDS_binned_plot_all_transcripts_mean.png",
    height = 1000, width = 1500, res = 300)
print(gplot_all)
dev.off()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/sieS25_legend.png",
    height = 100, width = 1000, res = 300)
cowplot::ggdraw(cowplot::get_legend(gplot_all))
dev.off()

summarised_binned_RPFs_up %>%
  filter(region == "CDS") %>% 
  mutate(grouping = grouping * 2) %>% 
  ggplot(aes(x = grouping, y = average_counts, colour = condition))+
  geom_line(size = 1)+
  #geom_smooth() +
  geom_point() +
  scale_color_manual(values = c("#FFC20A", "#88CCEE")) +
  #facet_wrap(~class) +
  scale_y_continuous(limits = c(0,NA)) +
  ggtitle("RPFs up (n = 222)") +
  xlab("CDS%") +
  ylab("mean CPMs") +
  publication_theme() -> gplot_RPFsUP

summarised_binned_RPFs_down %>%
  filter(region == "CDS") %>% 
  mutate(grouping = grouping * 2) %>% 
  ggplot(aes(x = grouping, y = average_counts, colour = condition))+
  geom_line(size = 1)+
  geom_point() +
  scale_color_manual(values = c("#FFC20A", "#88CCEE")) +
  #facet_wrap(~class) +
  scale_y_continuous(limits = c(0,NA)) +
  ggtitle("RPFs down (n = 141)") +
  xlab("CDS%") +
  ylab("mean CPMs") +
  publication_theme() -> gplot_RPFsdown

summarised_binned_totals_up %>%
  filter(region == "CDS") %>% 
  mutate(grouping = grouping * 2) %>% 
  ggplot(aes(x = grouping, y = average_counts, colour = condition))+
  geom_line(size = 1)+
  geom_point() +
  scale_color_manual(values = c("#FFC20A", "#88CCEE")) +
  #facet_wrap(~class) +
  scale_y_continuous(limits = c(0,NA)) +
  ggtitle("Totals up (n = 46)") +
  xlab("CDS%") +
  ylab("mean CPMs") +
  publication_theme() -> gplot_totalsup

summarised_binned_totals_down %>%
  filter(region == "CDS") %>% 
  mutate(grouping = grouping * 2) %>% 
  ggplot(aes(x = grouping, y = average_counts, colour = condition))+
  geom_line(size = 1)+
  geom_point() +
  scale_color_manual(values = c("#FFC20A", "#88CCEE")) +
  #facet_wrap(~class) +
  scale_y_continuous(limits = c(0,NA)) +
  ggtitle("Totals down (n = 72)") +
  xlab("CDS%") +
  ylab("mean CPMs") +
  publication_theme() -> gplot_totalsdown

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CDS_binned_plot_stability_mean.png",
    height = 1000, width = 2000, res = 300)
print(CDS_plot)
dev.off()

summarised_binned_both_up %>%
  filter(region == "CDS") %>% 
  mutate(grouping = grouping * 2) %>% 
  ggplot(aes(x = grouping, y = average_counts, colour = condition))+
  geom_line(size = 1)+
  geom_point() +
  scale_color_manual(values = c("#FFC20A", "#88CCEE")) +
  #facet_wrap(~class) +
  scale_y_continuous(limits = c(0,NA)) +
  ggtitle("Both up (n = 1027)") +
  xlab("CDS%") +
  ylab("mean CPMs") +
  publication_theme() -> gplot_bothup

summarised_binned_both_down %>%
  filter(region == "CDS") %>% 
  mutate(grouping = grouping * 2) %>% 
  ggplot(aes(x = grouping, y = average_counts, colour = condition))+
  geom_line(size = 1)+
  geom_point() +
  scale_color_manual(values = c("#FFC20A", "#88CCEE")) +
  #facet_wrap(~class) +
  scale_y_continuous(limits = c(0,NA)) +
  ggtitle("Both down (n = 998)") +
  xlab("CDS%") +
  ylab("mean CPMs") +
  publication_theme() -> gplot_bothdown


