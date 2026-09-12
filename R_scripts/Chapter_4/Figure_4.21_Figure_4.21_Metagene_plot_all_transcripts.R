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
                                      size = 26, hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            #axis.title = element_text(face = "bold",size = rel(1)),
            #axis.title.y = element_text(angle=90,vjust =2),
            #axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(size = 20, color = "black"),
            axis.text.y = element_text(size = 20, color = "black"),
            #axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = 10),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
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
  theme(legend.position="none")
        #axis.ticks.x = element_blank(),
        #axis.text.y = element_text(size = 18))
        #axis.text.x = element_blank())

CDS_theme <- my_theme+
  theme(legend.position="none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank())

UTR3_theme <- my_theme+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.line.y = element_blank())

#read in data----
load(file = file.path(parent_dir, "Counts_files/R_objects/M_binned_list.Rdata"))
summary(binned_list[[1]])
head(binned_list[[1]])

load(file = file.path(parent_dir, "Counts_files/R_objects/M_single_nt_list.Rdata"))
summary(single_nt_list[[1]])
head(single_nt_list[[1]])

region_lengths <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_region_lengths.csv", col_names = c("transcript", "UTR5_len", "CDS_len", "UTR3_len"))
most_abundant_transcripts <- read_csv(file = file.path("\\\\data.beatson.gla.ac.uk/data/R11/external_sequencing_data/Gillen_2021_RiboSeq/Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))

TM_containing_transcripts <- read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/extracted_TMs.csv")

TM_containing_transcripts <- 
  TM_containing_transcripts %>% 
  #filter(class != "sponly" & class != "type1" & class != "unknown") %>% 
  pull(var = "ENST") %>% unique()

binned_list2 <-  
  lapply(binned_list, function(x) {x %>% filter(!transcript %in% TM_containing_transcripts)})
binned_list3 <-  
  lapply(binned_list, function(x) {x %>% filter(transcript %in% TM_containing_transcripts)})

lapply(binned_list, function(x){x %>% pull(transcript) %>% unique() %>% length()})
lapply(binned_list2, function(x){x %>% pull(transcript) %>% unique() %>% length()})
lapply(binned_list3, function(x){x %>% pull(transcript) %>% unique() %>% length()})


# binned_list <-  
# lapply(binned_list, function(x) {x %>% filter(transcript %in% TM_containing_transcripts)})


#all transcripts----
#summarise
#summarise within each sample
summarised_binned_list <- lapply(binned_list, summarise_data, value = "binned_cpm", grouping = "bin")
summary(summarised_binned_list[[1]])
print(summarised_binned_list[[1]])

#summarise within each condition (across replicates)
do.call("rbind", summarised_binned_list) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned
summary(summarised_binned)
print(summarised_binned)

# "#E41A1C"  "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF"
#  "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494" "#B3B3B3"


plot_binned_lines <- 
  function(df, SD = F, conditions = NULL, mylabels = NULL, colours = NULL, CDS_only = F) {
  
  #order or filter conditions (if conditions argument used)
  if (!is.null(conditions)) {
    df %>%
      filter(condition %in% conditions) %>%
      mutate(condition = factor(condition, levels = conditions, labels = mylabels, ordered = T)) -> df
  }
  
  # #calculate axis limits
  # if (SD == F) {
  #   if (CDS_only == F) {
  #     ylims <- c(0,max(df$average_counts))
  #   } else {
  #     ylims <- c(0,max(df$average_counts[df$region == "CDS"]))
  #   }
  #   
  # } else {
  #   ylims <- c(min(df$average_counts - df$sd_counts),
  #              max(df$average_counts + df$sd_counts))
  # }
  # 
    
    ylims <- c(0,0.2)
    
  #5'UTR
  df[df$region == "UTR5",] %>%
    ggplot(aes(x = (grouping*4), y = average_counts, colour = condition))+
    geom_line(size = 1)+
    geom_point() +
    {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
    {if(!is.null(colours))scale_colour_manual(values=colours)}+
    {if(!is.null(colours))scale_fill_manual(values=colours)}+
    ggtitle("5'UTR") +
    UTR5_theme+
    #{if(!is.null(title))ggtitle("")}+
    ylim(ylims) -> UTR5_plot
  
  #CDS
  df[df$region == "CDS",] %>%
    ggplot(aes(x = (grouping*2), y = average_counts, colour = condition))+
    geom_line(size = 1)+
    geom_point() +
    {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
    {if(!is.null(colours))scale_colour_manual(values=colours)}+
    {if(!is.null(colours))scale_fill_manual(values=colours)}+
    ggtitle("CDS") +
    {if(CDS_only == F)CDS_theme}+
    {if(CDS_only)CDS_only_theme}+
    #{if(!is.null(title))ggtitle(title)}+
    ylim(ylims) -> CDS_plot
  
  #3'UTR
  df[df$region == "UTR3",] %>%
    ggplot(aes(x = grouping*4, y = average_counts, colour = condition))+
    geom_line(size = 1)+
    geom_point() +
    {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
    {if(!is.null(colours))scale_colour_manual(values=colours)}+
    {if(!is.null(colours))scale_fill_manual(values=colours)}+
    ggtitle("3'UTR") +
    #{if(!is.null(title))ggtitle("")}+
    ylim(ylims)+
    UTR3_theme -> UTR3_plot
  
  if (CDS_only) {
    return(CDS_plot) 
  } else {
    return(list(UTR5_plot, CDS_plot, UTR3_plot))
  }
}

plot_binned_lines2 <- 
  function(df, SD = F, conditions = NULL, mylabels = NULL, colours = NULL, CDS_only = F) {
    
    #order or filter conditions (if conditions argument used)
    if (!is.null(conditions)) {
      df %>%
        filter(condition %in% conditions) %>%
        mutate(condition = factor(condition, levels = conditions, labels = mylabels, ordered = T)) -> df
    }
    
    # #calculate axis limits
    # if (SD == F) {
    #   if (CDS_only == F) {
    #     ylims <- c(0,max(df$average_counts))
    #   } else {
    #     ylims <- c(0,max(df$average_counts[df$region == "CDS"]))
    #   }
    #   
    # } else {
    #   ylims <- c(min(df$average_counts - df$sd_counts),
    #              max(df$average_counts + df$sd_counts))
    # }
    # 
    
    ylims <- c(0,0.2)
    
    #5'UTR
    df[df$region == "UTR5",] %>%
      ggplot(aes(x = grouping, y = average_counts, colour = condition))+
      geom_line(size = 1)+
      geom_point() +
      {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
      {if(!is.null(colours))scale_colour_manual(values=colours)}+
      {if(!is.null(colours))scale_fill_manual(values=colours)}+
      #ggtitle("5'UTR") +
      UTR5_theme+
      #{if(!is.null(title))ggtitle("")}+
      ylim(ylims) -> UTR5_plot
    
    #CDS
    df[df$region == "CDS",] %>%
      ggplot(aes(x = grouping, y = average_counts, colour = condition))+
      geom_line(size = 1)+
      geom_point() +
      {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
      {if(!is.null(colours))scale_colour_manual(values=colours)}+
      {if(!is.null(colours))scale_fill_manual(values=colours)}+
      #ggtitle("CDS") +
      {if(CDS_only == F)CDS_theme}+
      {if(CDS_only)CDS_only_theme}+
      #{if(!is.null(title))ggtitle(title)}+
      ylim(ylims) -> CDS_plot
    
    #3'UTR
    df[df$region == "UTR3",] %>%
      ggplot(aes(x = grouping, y = average_counts, colour = condition))+
      geom_line(size = 1)+
      geom_point() +
      {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
      {if(!is.null(colours))scale_colour_manual(values=colours)}+
      {if(!is.null(colours))scale_fill_manual(values=colours)}+
      #ggtitle("3'UTR") +
      #{if(!is.null(title))ggtitle("")}+
      ylim(ylims)+
      UTR3_theme -> UTR3_plot
    
    if (CDS_only) {
      return(CDS_plot) 
    } else {
      return(list(UTR5_plot, CDS_plot, UTR3_plot))
    }
  }


#plot lines
binned_line_plots <- 
  plot_binned_lines(df = summarised_binned, 
                    SD = T, 
                    conditions = c(control, treatment), 
                    mylabels = c("Tot", "CNOT3"), 
                    colours = c("#FFC20A", "#0C7BDC"))

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts binned lines.png"),
    res = 300,
    width = 3000, 
    height = 1000)

CNOT3_RPFs_all <- 
grid.arrange(
  #top=textGrob("Monosomes",gp=gpar(fontsize=20,font=3)),
  binned_line_plots[[1]], 
  binned_line_plots[[2]], 
  binned_line_plots[[3]], 
  nrow = 2, 
  # ncol = 3,
  heights = c(32,1),
  widths = c(1.5,2,1.5))
             #)
dev.off()

### testing something out


#all transcripts----
#summarise
#summarise within each sample
summarised_binned_list2 <- lapply(binned_list2, summarise_data, value = "binned_cpm", grouping = "bin")
summary(summarised_binned_list2[[1]])
print(summarised_binned_list2[[1]])

#summarise within each condition (across replicates)
do.call("rbind", summarised_binned_list2) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned2
summary(summarised_binned2)
print(summarised_binned2)

# "#E41A1C"  "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF"
#  "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494" "#B3B3B3"

#plot lines
binned_line_plots2 <- 
  plot_binned_lines2(df = summarised_binned2, 
                    SD = T, 
                    conditions = c(control, treatment), 
                    mylabels = c("Tot", "CNOT3"), 
                    colours = c("#377EB8", "#66C2A5"))


CNOT3_RPFs_noTMs <- 
grid.arrange(
  binned_line_plots2[[1]], 
  binned_line_plots2[[2]], 
  binned_line_plots2[[3]], 
  nrow = 2, 
  # ncol = 3,
  heights = c(32,1),
  widths = c(1,2,1.5))
#)


### testing 3rd thing out

#all transcripts----
#summarise
#summarise within each sample
summarised_binned_list3 <- lapply(binned_list3, summarise_data, value = "binned_cpm", grouping = "bin")
summary(summarised_binned_list3[[1]])
print(summarised_binned_list3[[1]])

#summarise within each condition (across replicates)
do.call("rbind", summarised_binned_list3) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned3
summary(summarised_binned3)
print(summarised_binned3)

# "#E41A1C"  "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF"
#  "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494" "#B3B3B3"

#plot lines
binned_line_plots3 <- 
  plot_binned_lines2(df = summarised_binned3, 
                    SD = T, 
                    conditions = c(control, treatment), 
                    mylabels = c("Tot", "CNOT3"), 
                    colours = c("#377EB8", "#66C2A5"))


CNOT3_RPFs_TMs <- 
grid.arrange(
  binned_line_plots3[[1]], 
  binned_line_plots3[[2]], 
  binned_line_plots3[[3]], 
  nrow = 2, 
  # ncol = 3,
  heights = c(32,1),
  widths = c(1,2,1.5))
#)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(CNOT3_RPFs_all)

grid.arrange(
  CNOT3_RPFs_all, 
  CNOT3_RPFs_TMs,
  CNOT3_RPFs_noTMs,
  nrow = 3,
  ncol = 1)
#)





#######
#back to original script
#####


binned_line_plots_all_replicates <- plot_binned_all_replicates(summarised_binned_list, 
                                                               conditions = c(control, treatment), 
                                                               mylabels = c("Tot", "CNOT3"),
                                                               colours = c("#377EB8", "#66C2A5"))

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts binned lines all replicates.png"), width = 1000, height = 300)
grid.arrange(
  top=textGrob("Monosomes",gp=gpar(fontsize=20,font=3)),
  binned_line_plots_all_replicates[[1]], 
  binned_line_plots_all_replicates[[2]], 
  binned_line_plots_all_replicates[[3]], 
  nrow = 2, 
  # ncol = 3,
  heights = c(64,1),
  widths = c(1,2,1.5))
dev.off()

#calculate and plot delta
binned_delta_data <- calculate_binned_delta(binned_list, value = "binned_cpm", control = control, treatment = treatment, paired_data = F)


binned_delta_plots <- plot_binned_delta(binned_delta_data)

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts binned delta.png"), width = 1000, height = 200)
grid.arrange(
  top=textGrob("Monosomes",gp=gpar(fontsize=20,font=3)),
  binned_delta_plots[[1]], 
  binned_delta_plots[[2]], 
  binned_delta_plots[[3]],
  nrow = 2,
  heights = c(64,1),
  widths = c(1,2,1))
dev.off()

#positional----
#normalise within each transcript
positional_list <- lapply(binned_list, calculate_positional_counts)

#summarise within each sample
summarised_positional_list <- lapply(positional_list, summarise_data, value = "positional_counts", grouping = "bin")

#summarise within each condition (across replicates)
do.call("rbind", summarised_positional_list) %>%
  group_by(grouping, condition) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_positional

#plot
positional_line_plots <- plot_positional_lines(df = summarised_positional, 
                                               SD = T, 
                                               conditions = c(control, treatment),
                                               mylabels = c("Tot", "CNOT3"),
                                               title = "CDS",
                                               colours = c("#377EB8", "#66C2A5"))

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts binned positional lines_M.png"), width = 500, height = 200)
print(positional_line_plots)
dev.off()

#calculate and plot delta
binned_positional_delta <- calculate_positional_delta(positional_list, control = control, treatment = treatment, paired_data = F)
positional_binned_delta_plots <- plot_positional_delta(binned_positional_delta)

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts binned positional delta_M.png"), width = 500, height = 200)
print(positional_binned_delta_plots)
dev.off()

#single nt----
#summarise
summarised_single_nt_list <- lapply(single_nt_list, summarise_data, value = "single_nt_cpm", grouping = "window")
summary(summarised_single_nt_list[[1]])
print(summarised_single_nt_list[[1]])

do.call("rbind", summarised_single_nt_list) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_single_nt
summary(summarised_single_nt)
print(summarised_single_nt)

#plot
single_nt_line_plots <- plot_single_nt_lines(summarised_single_nt, 
                                             SD=T, 
                                             plot_ends=F, 
                                             conditions = c(control, treatment),
                                             mylabels = c("Tot", "CNOT3"),
                                             colours = c("#377EB8", "#66C2A5"))

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts single nt lines_M.png"), width = 1300, height = 300)
grid.arrange(
  top=textGrob("Monosomes",gp=gpar(fontsize=20,font=3)),
  single_nt_line_plots[[1]], 
  single_nt_line_plots[[2]], 
  single_nt_line_plots[[3]], 
  single_nt_line_plots[[4]], 
  nrow = 2, 
  heights = c(64,1),
  widths = c(1,2,2,1.5))
dev.off()

#calculate and plot delta
single_nt_delta_data <- calculate_single_nt_delta(single_nt_list, 
                                                  value = "single_nt_cpm", 
                                                  control = control, 
                                                  treatment = treatment, 
                                                  paired_data = F)
single_nt_delta_plots <- plot_single_nt_delta(single_nt_delta_data, SD = T)

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts single nt delta_M.png"), width = 1300, height = 200)
grid.arrange(
  top=textGrob("Monosomes",gp=gpar(fontsize=20,font=3)),
  single_nt_delta_plots[[1]], 
  single_nt_delta_plots[[2]], 
  single_nt_delta_plots[[3]], 
  single_nt_delta_plots[[4]], 
  nrow = 2, 
  heights = c(64,1),
  widths = c(1,2,2,1))
dev.off()

#Dep vs Antidep vs Indep----
#read in DESeq2 output
DESeq2_data <- 
  read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output/RPFs_CNOT3_DEseq2_apeglm_LFC_shrinkage_M.csv")) %>% 
  mutate(RPFs_group = factor(case_when(log2FoldChange <= 0 & padj < 0.05 ~ "downregulated",
                                       log2FoldChange >= 0 & padj < 0.05 ~ "upregulated",
                                       TRUE ~ "unchanged")))

summary(DESeq2_data)


#extract transcript IDs
RPFs_down_IDs <- DESeq2_data$transcript[DESeq2_data$RPFs_group == "downregulated" & !(is.na(DESeq2_data$RPFs_group))]
RPFs_up_IDs <- DESeq2_data$transcript[DESeq2_data$RPFs_group == "upregulated" & !(is.na(DESeq2_data$RPFs_group))]
no_change_IDs <- DESeq2_data$transcript[DESeq2_data$RPFs_group == "unchanged" & !(is.na(DESeq2_data$TE_group))]

#plot binned
plot_subset(IDs = RPFs_up_IDs, subset = "upregulated", sub_dir = "Upregulated",
            control = control, treatment = treatment,
            binned_value = "binned_cpm", single_nt_value = "single_nt_cpm",
            plot_binned = T, plot_single_nt = T, plot_positional = T,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = T, mylabels = c("Tot", "CNOT3"))




#plot heatmaps
TE_down_heatmap <- plot_binned_heatmaps(IDs = RPFs_up_IDs, col_lims = c(-0.02, 0.01), control = control, treatment = treatment, value = "binned_cpm")
TE_up_heatmap <- plot_binned_heatmaps(IDs = TE_up_IDs, col_lims = c(-0.01, 0.02), control = control, treatment = treatment, value = "binned_cpm")

png(filename = file.path(parent_dir, "plots/binned_plots/Dep", paste(treatment, "TE-down binned heatmap.png")), width = 1000, height = 1000)
print(TE_down_heatmap)
dev.off()

png(filename = file.path(parent_dir, "plots/binned_plots/Dep", paste(treatment, "TE-up binned heatmap.png")), width = 1000, height = 1000)
print(TE_up_heatmap)
dev.off()

#GSEA pathways----
library(fgsea)

#read in pathways
source("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/read_human_GSEA_pathways.R")

#hallmark
#read in fgsea output
load(file = file.path(parent_dir, "Analysis/fgsea/hallmark_results.Rdata"))

#extract transcript IDs
hallmark_pathways <- hallmark_results[[3]]$pathway[hallmark_results[[3]]$padj < 0.05]

lapply(hallmark_pathways, plot_GSEA_binned,
       GSEA_set = pathways.hallmark, sub_dir = "hallmark",
       control = control, treatment = treatment,
       binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
       plot_binned = T, plot_single_nt = F, plot_positional = F,
       plot_delta = T, SD = T, paired_data = F)



