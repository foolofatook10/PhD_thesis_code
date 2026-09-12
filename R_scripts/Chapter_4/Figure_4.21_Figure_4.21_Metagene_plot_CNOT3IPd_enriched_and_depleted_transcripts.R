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
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            #axis.title = element_text(face = "bold",size = rel(1)),
            #axis.title.y = element_text(angle=90,vjust =2),
            #axis.title.x = element_text(vjust = -0.2),
            #axis.text.x = element_text(size = 16, color = "black"),
            #axis.text.y = element_text(size = 16, color = "black"),
            #axis.text = element_text(), 
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
deseq <- 
  read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/Analysis/DESeq2_output/RPFs_CNOT3_DEseq2_apeglm_LFC_shrinkage_M_28nts.csv") %>% 
  mutate(grouping = case_when(padj < 0.05 & log2FoldChange > 0  ~ "enriched",
                              padj < 0.05 & log2FoldChange < 0 ~ "depleted",
                              .default = "unchanged"))

enriched_transcripts <-  deseq %>% filter(grouping == "enriched") %>% pull(transcript)
depleted_transcripts <-  deseq %>% filter(grouping == "depleted") %>% pull(transcript)

binned_list_enriched <-  
  lapply(binned_list, function(x) {x %>% filter(transcript %in% enriched_transcripts)})
binned_list_depleted <-  
  lapply(binned_list, function(x) {x %>% filter(transcript %in% depleted_transcripts)})

lapply(binned_list, function(x){x %>% pull(transcript) %>% unique() %>% length()})
lapply(binned_list_enriched, function(x){x %>% pull(transcript) %>% unique() %>% length()})
lapply(binned_list_depleted, function(x){x %>% pull(transcript) %>% unique() %>% length()})




#all transcripts----
#summarise
#summarise within each sample
summarised_binned_list_enriched <- lapply(binned_list_enriched, summarise_data, value = "binned_cpm", grouping = "bin")
summary(summarised_binned_list_enriched[[1]])
print(summarised_binned_list_enriched[[1]])

#summarise within each condition (across replicates)
do.call("rbind", summarised_binned_list_enriched) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned_enriched
summary(summarised_binned_enriched)
print(summarised_binned_enriched)

enriched_CDS <- 
summarised_binned_enriched %>% filter(region == "CDS") %>% 
  mutate(class = "enriched")


#all transcripts----
#summarise
#summarise within each sample
summarised_binned_list_depleted <- 
  lapply(binned_list_depleted, 
         summarise_data, 
         value = "binned_cpm", 
         grouping = "bin")
summary(summarised_binned_list_depleted[[1]])
print(summarised_binned_list_depleted[[1]])

#summarise within each condition (across replicates)
do.call("rbind", summarised_binned_list_depleted) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned_depleted

summarised_binned_depleted %>% filter(region == "CDS") %>% pull(average_counts)

summary(summarised_binned_depleted)
print(summarised_binned_depleted)


depleted_CDS <- 
  summarised_binned_depleted %>% filter(region == "CDS") %>% 
  mutate(class = "depleted")

df <- 
rbind(enriched_CDS, depleted_CDS)

df$class <- factor(df$class, levels = c("depleted", "enriched"),
                   labels = c("Deseq depleted", "Deseq enriched")) 

df %>%
  mutate(grouping = grouping * 2) %>% 
  ggplot(aes(x = grouping, y = average_counts, colour = condition))+
  geom_line(size = 1)+
  geom_point() +
  scale_color_manual(values = c("#0C7BDC", "#FFC20A")) +
  facet_wrap(~class, scales = "free") +
  scale_y_continuous(limits = c(0,NA)) +
  ggtitle("CDS") +
  xlab("CDS%") +
  ylab("CPMs") +
  publication_theme() +
  theme(legend.title = element_blank(), legend.position = "none", legend.direction = "vertical")-> CDS_plot

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CDS_binned_plot_deseq.png",
    height = 1000, width = 2000, res = 300)
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
