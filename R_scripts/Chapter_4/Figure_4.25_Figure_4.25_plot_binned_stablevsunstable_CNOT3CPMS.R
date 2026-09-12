#load packages----
library(tidyverse)
library(grid)
library(gridExtra)
#library(viridis)

#read in common variables----
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/R_scripts/common_variables.R")

#set what you have called your control and treated samples. This can be a vector of strings if more than one treatment has been used.
control <- "M_Tot"
treatment <- "M_CNOT3"

#read in functions----
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/R_scripts/meta_plots/binning_Sel-RiboSeq_functions.R")

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
            axis.text.x = element_text(size = 16, color = "black"),
            axis.text.y = element_text(size = 16, color = "black"),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = 10),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(1, "cm"),
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
load(file = file.path(parent_dir, "Counts_files/R_objects/M_binned_list.Rdata"))
summary(binned_list[[1]])
head(binned_list[[1]])

load(file = file.path(parent_dir, "Counts_files/R_objects/M_single_nt_list.Rdata"))
summary(single_nt_list[[1]])
head(single_nt_list[[1]])

region_lengths <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_region_lengths.csv", col_names = c("transcript", "UTR5_len", "CDS_len", "UTR3_len"))
most_abundant_transcripts <- read_csv(file = file.path("\\\\data.beatson.gla.ac.uk/data/R11/external_sequencing_data/Gillen_2021_RiboSeq/Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))

# read in deseq2
stability <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/stability_extremes.csv") 

stable_transcripts <-  stability %>% filter(stability == "stable") %>% pull(ENST)
unstable_transcripts <-  stability %>% filter(stability == "unstable") %>% pull(ENST)

binned_list_stable <-  
  lapply(binned_list, function(x) {x %>% filter(transcript %in% stable_transcripts)})
binned_list_unstable <-  
  lapply(binned_list, function(x) {x %>% filter(transcript %in% unstable_transcripts)})

lapply(binned_list, function(x){x %>% pull(transcript) %>% unique() %>% length()})
lapply(binned_list_stable, function(x){x %>% pull(transcript) %>% unique() %>% length()})
lapply(binned_list_unstable, function(x){x %>% pull(transcript) %>% unique() %>% length()})

tpms <- 
  read_csv("\\\\data.beatson.gla.ac.uk/data/R11/external_sequencing_data/Gillen_2021_RiboSeq/Analysis/DESeq2_output/tpms.csv")

tpms_long <- 
  tpms %>% 
  gather(key = "sample", value = tpm, Ctrl_1_Totals:CNOT1_3_Totals) 

samples_split <- str_split(tpms_long$sample, "_")
condition <- unlist(lapply(samples_split, function(x){x[[1]]}))

tpms_long$condition <- condition

ctrl_tpms <- tpms_long %>% 
  filter(condition == "Ctrl") %>% 
  group_by(transcript) %>% 
  summarise(mean_tpm = mean(tpm))

summarise_data <- function(df, value, grouping) {
  df %>%
    inner_join(ctrl_tpms, by = "transcript") %>% 
    mutate(nCPM = binned_cpm/mean_tpm) %>% 
    dplyr::rename(value = "nCPM") %>% 
    dplyr::rename(grouping = "bin") %>%
    group_by(region, condition, replicate, grouping) %>%
    summarise(mean_counts = mean(value),
              median_counts = median(value)) %>%
    ungroup() -> summarised_df
  
  return(summarised_df)
}




#all transcripts----
#summarise
#summarise within each sample
summarised_binned_list_stable <- lapply(binned_list_stable, summarise_data, value = "ncpm", grouping = "bin")
summary(summarised_binned_list_stable[[1]])
print(summarised_binned_list_stable[[1]])

#summarise within each condition (across replicates)
do.call("rbind", summarised_binned_list_stable) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned_stable
summary(summarised_binned_stable)
print(summarised_binned_stable)

stable_CDS <- 
  summarised_binned_stable %>% filter(region == "CDS") %>% 
  mutate(class = "stable")


#all transcripts----
#summarise
#summarise within each sample
summarised_binned_list_unstable <- 
  lapply(binned_list_unstable, 
         summarise_data, 
         value = "ncpm", 
         grouping = "bin")
summary(summarised_binned_list_unstable[[1]])
print(summarised_binned_list_unstable[[1]])

#summarise within each condition (across replicates)
do.call("rbind", summarised_binned_list_unstable) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned_unstable

summarised_binned_unstable %>% filter(region == "CDS") %>% pull(average_counts)

summary(summarised_binned_unstable)
print(summarised_binned_unstable)


unstable_CDS <- 
  summarised_binned_unstable %>% filter(region == "CDS") %>% 
  mutate(class = "unstable")

df <- 
  rbind(stable_CDS, unstable_CDS)

df$class <- factor(df$class, levels = c("unstable", "stable"),
                   labels = c("Unstable transcripts\n(n = 1059)", "Stable transcripts\n(n = 1207)")) 

length(unstable_transcripts)
length(stable_transcripts)
df %>%
  mutate(grouping = grouping * 2) %>% 
  ggplot(aes(x = grouping, y = average_counts, colour = condition))+
  geom_line(size = 1)+
  geom_point() +
  scale_color_manual(values = c("#0C7BDC", "#FFC20A")) +
  facet_wrap(~class) +
  scale_y_continuous(limits = c(0,NA)) +
  ggtitle("CDS") +
  xlab("CDS%") +
  ylab("nCPMs") +
  publication_theme() +
  theme(legend.title = element_blank(), legend.position = "none", legend.direction = "vertical")-> CDS_plot

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CDS_binned_plot_stability.png",
    height = 1250, width = 2000, res = 300)
print(CDS_plot)
dev.off()

#plot lines
binned_line_plots <- 
  plot_binned_lines(df = df, 
                    SD = T, 
                    conditions = c(control, treatment), 
                    mylabels = c("Tot", "CNOT3"), 
                    colours = c(, ))



png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3SelRP_metagene_deseq_depleted.png", res = 300, height = 1000, width = 2000) 
grid.arrange(
  binned_line_plots_enriched[[2]], 
  binned_line_plots_dep[[2]], 
  #binned_line_plots_dep[[3]], 
  nrow = 2, 
  # ncol = 3,
  heights = c(32,1),
  widths = c(2,2))
dev.off()
