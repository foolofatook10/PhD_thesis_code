#load packages----
library(tidyverse)
library(grid)
library(gridExtra)

#read in common variables----
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/R_scripts/common_variables.R")

#read in data----
region_lengths <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_region_lengths.csv", col_names = c("transcript", "UTR5_len", "CDS_len", "UTR3_len"))

#create themes----
my_theme <- theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 16),
        legend.position="none")

#create themes----
publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = 36, hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = 20),
            axis.title.y = element_text(angle=90,vjust =2, size = 24),
            axis.title.x = element_text(face = "bold", vjust = -0.2, size = 25),
            axis.text.x = element_text(size = 24, color = "black"),
            axis.text.y = element_text(size = 24, color = "black"),
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
            legend.title = element_text(face="italic", size = 24),
            plot.margin=unit(c(2,4,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

#Read in the counts csvs, calculate Counts per Million (CPM) and splice into first and last 25/50nt of CDS/UTRs----
spliced_list <- list()
for (sample in RPF_sample_names) {
  
  #read in counts csv
  df <- read_csv(file = file.path(parent_dir, "Counts_files/csv_files", paste0(sample, "_pc_final_counts.csv")))
  
  total_counts <- sum(df$Counts)
  
  df %>%
    mutate(CPM = (Counts / total_counts) * 1000000)  %>%
    inner_join(region_lengths, by = "transcript") -> merged_data
  
  #splice the transcript
  #UTR5 end
  merged_data %>%
    filter(Position <= UTR5_len) %>%
    group_by(transcript) %>%
    top_n(n = 25, wt = Position) %>% #extracts the 3' most nts
    ungroup() %>%
    mutate(nt = (Position - UTR5_len) - 1) %>%
    group_by(nt) %>%
    summarise(mean_cpm = mean(CPM)) %>%
    ungroup() %>%
    mutate(region = rep("UTR5")) -> UTR5_end
  
  #CDS start
  merged_data %>%
    filter(Position > UTR5_len & Position <= (UTR5_len + CDS_len)) %>%
    mutate(nt = Position - UTR5_len) %>%
    filter(nt <= 50) %>%
    group_by(nt) %>%
    summarise(mean_cpm = mean(CPM)) %>%
    ungroup() %>%
    mutate(region = rep("CDS")) -> CDS_start
  
  #CDS end
  merged_data %>%
    filter(Position > UTR5_len & Position <= (UTR5_len + CDS_len)) %>%
    group_by(transcript) %>%
    top_n(n = 50, wt = Position) %>% #extracts the 3' most nts
    ungroup() %>%
    mutate(nt = (Position - (UTR5_len + CDS_len)) - 1) %>%
    group_by(nt) %>%
    summarise(mean_cpm = mean(CPM)) %>%
    ungroup() %>%
    mutate(region = rep("CDS")) -> CDS_end
  
  #UTR3 start
  merged_data %>%
    filter(Position > (UTR5_len + CDS_len)) %>%
    mutate(nt = Position - (UTR5_len + CDS_len)) %>%
    filter(nt <= 25) %>%
    group_by(nt) %>%
    summarise(mean_cpm = mean(CPM)) %>%
    ungroup() %>% 
    mutate(region = rep("UTR3")) -> UTR3_start
  
  bind_rows(UTR5_end, CDS_start, CDS_end, UTR3_start) %>%
    mutate(sample = rep(sample)) -> spliced_list[[sample]]
}

all_data <- do.call("rbind", spliced_list)
  
names <- str_split(all_data$sample, "_")

all_data$replicate <- as.numeric(str_remove(unlist(lapply(names, function(x){x[[1]]})), "REP"))
all_data$condition <- unlist(lapply(names, function(x){x[[2]]}))

all_data <- all_data %>% 
  mutate(condition = case_when(condition == "siRPS25" ~"sieS25",
                               .default = condition))

samples_split <- split(all_data, all_data$condition)

samples_split_processing <- function(x){
  
  y <- 
  x %>% 
    select(-sample) %>% 
    group_by(condition, region, nt) %>% 
    summarise(cpm = mean(mean_cpm),
              sd_cpm = sd(mean_cpm),
              se = sd_cpm/sqrt(5))
  
  return(y)
    
}

df <- samples_split_processing(samples_split[[1]])
table(is.na(samples_split[[1]]$mean_cpm))

plot_gplots <- function(df){
  
  cond <- unique(df$condition)
  
  if(cond == "NTC"){
    line_col = "#FFC20A"
  } else if(cond == "siCNOT3"){
    line_col = "#88CCEE"
  } else{
    line_col = "#984EA3"
  }
  
  ylims <- c(0,max(df$mean_cpm))
  
  # tytle <- 
  #   ggplot() +
  #   ggtitle(unique(df$condition)) +
  #   publication_theme() +
  #   theme(plot.margin = unit(c(1,0,0,0), "mm"))
  
  UTR5_end_plot <- 
    df %>%
    filter(region == "UTR5" & nt < 0) %>% 
    ggplot(aes(x = nt, y = mean_cpm)) +
    ylim(ylims)+
    xlab("nt\n(rel. to start codon)")+
    #ylab("mean CPM") +
    stat_summary(fun = mean, geom = "line", colour = line_col, size = 1) +
    stat_summary(fun = mean, geom = "point", colour = line_col) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5, colour = "black") +
    ggtitle("5'UTR") +
    publication_theme() +
    theme(axis.title.y = element_blank())
  
  #CDS
  df %>% 
    filter(region == "CDS" & nt > 0) %>%
    ggplot(aes(x = nt, y = mean_cpm))+
    #geom_line(size = 1, colour = line_col)+
    ylim(ylims)+
    xlab("nt\n(rel. to start codon)")+
    ggtitle("5'CDS") +
    stat_summary(fun = mean, geom = "line", colour = line_col, size = 1) +
    stat_summary(fun = mean, geom = "point", colour = line_col) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5, colour = "black") +
    publication_theme()+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) -> CDS_start_plot
  
  df %>% 
    filter(region == "CDS" & nt < 0) %>%
    ggplot(aes(x = nt, y = mean_cpm))+
    #geom_line(size = 1, colour = line_col)+
    ylim(ylims)+
    xlab("nt\n(rel. to stop codon)")+
    ggtitle("3'CDS") +
    stat_summary(fun = mean, geom = "line", colour = line_col, size = 1) +
    stat_summary(fun = mean, geom = "point", colour = line_col) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5, colour = "black") +
    publication_theme()+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) -> CDS_end_plot
  
  #3'UTR
  df%>% 
    filter(region == "UTR3" & nt > 0) %>% 
    ggplot(aes(x = nt, y = mean_cpm))+
    #geom_line(size = 1, colour = line_col)+
    ylim(ylims)+
    xlab("nt\n(rel. to stop codon)")+
    ggtitle("3'UTR") +
    stat_summary(fun = mean, geom = "line", colour = line_col, size = 1) +
    stat_summary(fun = mean, geom = "point", colour = line_col) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5, colour = "black") +
    publication_theme()+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) -> UTR3_start_plot
  
  gplot1 <- cowplot::plot_grid(UTR5_end_plot, 
                               CDS_start_plot, 
                               CDS_end_plot, 
                               UTR3_start_plot,
                               nrow = 1,
                               rel_widths = c(0.6,1,1,0.4))
  
  
  #gplot2 <- cowplot::plot_grid(tytle, gplot1, nrow = 2, rel_heights = c(0.05, 1))
  return(gplot1)
  
}

gplots <- lapply(samples_split, plot_gplots)

vec_names <- c("NTC", "siCNOT3", "sieS25")

for(i in 1:length(vec_names)){
  
  png(filename = file.path(parent_dir, 
                           "plots/offset_aligned_single_nt_plots",
                           paste0(vec_names[[i]],".png")),
      res = 300, height = 1250, width = 7000)
  print(gplots[[i]])
  dev.off()
  
}

gplot_for_legend <- 
all_data %>% 
  filter(region == "UTR5" & nt < 0) %>% 
  group_by(condition, nt) %>% 
  summarise(CPM = mean(mean_cpm)) %>% 
  mutate(condition = factor(condition, levels = c("NTC",
                                                  "siCNOT3",
                                                  "sieS25"))) %>% 
  ggplot(aes(x = nt, y = CPM, colour = condition)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#FFC20A", "#88CCEE", "#984EA3")) +
  publication_theme() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size= unit(0.8, "cm"),
        legend.margin = unit(0.1, "cm"))

leg <- cowplot::get_legend(gplot_for_legend)

png(filename = file.path(parent_dir, "plots/offset_aligned_single_nt_plots", "legend.png"), width = 1000, height = 100, res = 300)
print(cowplot::ggdraw(leg))
dev.off()

gplots[[3]]

png(filename = file.path(parent_dir, "plots/offset_aligned_single_nt_plots", paste(sample, "offset aligned single nt plot.png")), width = 5500, height = 900, res = 300)
plot_gplots(samples_split[[1]])


df[df$region == "UTR5" & df$nt < 0,] %>%
  ggplot(aes(x = nt, y = cpm))+
  geom_line(size = 1, colour = "#FFC20A")+
  geom_errorbar(aes(ymin = cpm - se, ymax = cpm + se), width = 0.1)
  #ylim(ylims)+
  xlab("nt\n(relative to start codon)")+
  ylab("mean CPM")+
  publication_theme() -> UTR5_end_plot


sample = "REP4_NTC_M"
#plot samples individually----
for (sample in RPF_sample_names) {
  
  df <- spliced_list[[sample]]
  
  cond <- unlist(str_split(sample, "_"))[2]
  
  if(cond == "NTC"){
    line_col = "#FFC20A"
  } else if(cond == "siCNOT3"){
    line_col = "#88CCEE"
  } else{
    line_col = "#984EA3"
  }
  
  ylims <- c(0,max(df$mean_cpm))
  
  tytle <- 
    ggplot() +
    ggtitle(str_remove(str_replace_all(unique(paste0(df$sample)), "_", " "), "\\sM$")) +
    publication_theme()
  
  print(tytle)
  
  #5'UTR
  df[df$region == "UTR5" & df$nt < 0,] %>%
    ggplot(aes(x = nt, y = mean_cpm))+
    geom_line(size = 1, colour = line_col)+
    ylim(ylims)+
    xlab("nt\n(relative to start codon)")+
    ylab("mean CPM")+
    publication_theme() -> UTR5_end_plot
  
  #CDS
  df[df$region == "CDS" & df$nt > 0,] %>%
    ggplot(aes(x = nt, y = mean_cpm))+
    geom_line(size = 1, colour = line_col)+
    ylim(ylims)+
    xlab("nt\n(relative to start codon)")+
    publication_theme()+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) -> CDS_start_plot
  
  df[df$region == "CDS" & df$nt < 0,] %>%
    ggplot(aes(x = nt, y = mean_cpm))+
    geom_line(size = 1, colour = line_col)+
    ylim(ylims)+
    xlab("nt\n(relative to stop codon)")+
    publication_theme()+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) -> CDS_end_plot
  
  #3'UTR
  df[df$region == "UTR3" & df$nt > 0,] %>%
    ggplot(aes(x = nt, y = mean_cpm))+
    geom_line(size = 1, colour = line_col)+
    ylim(ylims)+
    xlab("nt\n(relative to stop codon)")+
    publication_theme()+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) -> UTR3_start_plot
  
  png(filename = file.path(parent_dir, "plots/offset_aligned_single_nt_plots", paste(sample, "offset aligned single nt plot.png")), width = 5500, height = 900, res = 300)
    grid.arrange(UTR5_end_plot,
                            CDS_start_plot,
                            CDS_end_plot,
                            UTR3_start_plot,
                            nrow = 1,
                            widths = c(1,2,2,1))
  dev.off()
}

#all mean of all samples----
do.call("rbind", spliced_list) %>%
  group_by(region, nt) %>%
  summarise(mean_cpm = mean(mean_cpm)) -> all_data

ylims <- c(0,max(all_data$mean_cpm))

#5'UTR
all_data[all_data$region == "UTR5" & all_data$nt < 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to start codon)")+
  ylab("mean CPM")+
  my_theme -> UTR5_end_plot

#CDS
all_data[all_data$region == "CDS" & all_data$nt > 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to start codon)")+
  my_theme+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) -> CDS_start_plot

all_data[all_data$region == "CDS" & all_data$nt < 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to stop codon)")+
  my_theme+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) -> CDS_end_plot

#3'UTR
all_data[all_data$region == "UTR3" & all_data$nt > 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to stop codon)")+
  my_theme+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) -> UTR3_start_plot

png(filename = file.path(parent_dir, "plots/offset_aligned_single_nt_plots/all samples offset aligned single nt plot.png"), width = 1300, height = 300)
grid.arrange(UTR5_end_plot, CDS_start_plot, CDS_end_plot, UTR3_start_plot, nrow = 1, widths = c(1,2,2,1))
dev.off()
