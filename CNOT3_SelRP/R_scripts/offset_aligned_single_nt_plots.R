#load packages----
library(tidyverse)
library(grid)
library(gridExtra)

#read in common variables----
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/R_scripts/common_variables.R")

#read in data----
region_lengths <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_region_lengths.csv", col_names = c("transcript", "UTR5_len", "CDS_len", "UTR3_len"))

#create themes----
my_theme <- theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 16),
        legend.position="none")

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
            axis.title = element_text(face = "bold",size = 20),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(size = 25, color = "black"),
            axis.text.y = element_text(size = 25, color = "black"),
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
            plot.margin=unit(c(4,2,2,2),"mm"),
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

sample = RPF_sample_names[[1]]

#plot samples individually----
for (sample in RPF_sample_names) {
  df <- spliced_list[[sample]]
  
  cond <- unlist(str_split(sample, "_"))[3]
  
  if(cond == "Tot"){
    line_col = "#FFC20A"
  } else if(cond == "CNOT3"){
    line_col = "#0C7BDC"
  } 
  
  ylims <- c(0,max(df$mean_cpm))
  
  REP <- as.numeric(str_remove(unlist(str_split(sample, "_"))[1], "REP")) -1
  tit <- paste0("REP", REP, " ", cond)
  
  tytle <- 
    ggplot() +
    ggtitle(paste0(tit)) +
    publication_theme()
  
  #5'UTR
  df[df$region == "UTR5" & df$nt < 0,] %>%
    ggplot(aes(x = nt, y = mean_cpm))+
    geom_line(size = 1, colour = line_col)+
    ylim(ylims)+
    xlab("nt\n(relative to start codon)")+
    ylab("mean CPM")+
    publication_theme() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank()) -> UTR5_end_plot
  
  #CDS
  df[df$region == "CDS" & df$nt > 0,] %>%
    ggplot(aes(x = nt, y = mean_cpm))+
    geom_line(size = 1, colour = line_col)+
    ylim(ylims)+
    xlab("nt\n(relative to start codon)")+
    publication_theme()+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
  axis.title.x = element_blank()) -> CDS_start_plot
  
  df[df$region == "CDS" & df$nt < 0,] %>%
    ggplot(aes(x = nt, y = mean_cpm))+
    geom_line(size = 1, colour = line_col)+
    ylim(ylims)+
    xlab("nt\n(relative to stop codon)")+
    publication_theme()+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
  axis.title.x = element_blank()) -> CDS_end_plot
  
  #3'UTR
  df[df$region == "UTR3" & df$nt > 0,] %>%
    ggplot(aes(x = nt, y = mean_cpm))+
    geom_line(size = 1, colour = line_col)+
    ylim(ylims)+
    xlab("nt\n(relative to stop codon)")+
    publication_theme()+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
  axis.title.x = element_blank()) -> UTR3_start_plot
  
  png(filename = file.path(parent_dir, "plots/offset_aligned_single_nt_plots", paste(tit, "offset aligned single nt plot.png")), res = 300, width = 5500, height = 1250)

    grid.arrange(UTR5_end_plot,
                 CDS_start_plot,
                 CDS_end_plot,
                 UTR3_start_plot,
                 nrow = 1,
                 widths = c(1.5,2,2,1.5))
    
  dev.off()
}

#all mean of all samples----
avdat <- do.call("rbind", spliced_list)

avdat$MD30 <- factor(str_remove_all(str_extract(avdat$sample, "_.*_"), "_"),
                     levels = c("M", "D30"))

avdat %>%
  group_by(MD30, region, nt) %>%
  summarise(mean_cpm = mean(mean_cpm)) -> all_data


all_data_M <- all_data %>% filter(MD30 == "M")
all_data_D30 <- all_data %>% filter(MD30 == "D30")

ylims <- c(0,max(all_data$mean_cpm))

#5'UTR
all_data_M[all_data_M$region == "UTR5" & all_data_M$nt < 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to start codon)")+
  ylab("mean CPM")+
  my_theme -> UTR5_end_plot

#CDS
all_data_M[all_data_M$region == "CDS" & all_data_M$nt > 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to start codon)")+
  my_theme+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) -> CDS_start_plot

all_data_M[all_data_M$region == "CDS" & all_data_M$nt < 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to stop codon)")+
  my_theme+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) -> CDS_end_plot

#3'UTR
all_data_M[all_data_M$region == "UTR3" & all_data_M$nt > 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to stop codon)")+
  my_theme+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) -> UTR3_start_plot

png(filename = file.path(parent_dir, "plots/offset_aligned_single_nt_plots/all samples offset aligned single nt plot M.png"), width = 1300, height = 300)
grid.arrange(UTR5_end_plot, CDS_start_plot, CDS_end_plot, UTR3_start_plot, nrow = 1, widths = c(1,2,2,1),
             top = textGrob("Monosomes",gp=gpar(fontsize=20,font=3)))
dev.off()













#5'UTR
all_data_D30[all_data_D30$region == "UTR5" & all_data_D30$nt < 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to start codon)")+
  ylab("mean CPM")+
  my_theme -> UTR5_end_plot

#CDS
all_data_D30[all_data_D30$region == "CDS" & all_data_D30$nt > 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to start codon)")+
  my_theme+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) -> CDS_start_plot

all_data_D30[all_data_D30$region == "CDS" & all_data_D30$nt < 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to stop codon)")+
  my_theme+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) -> CDS_end_plot

#3'UTR
all_data_D30[all_data_D30$region == "UTR3" & all_data_D30$nt > 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to stop codon)")+
  my_theme+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) -> UTR3_start_plot

png(filename = file.path(parent_dir, "plots/offset_aligned_single_nt_plots/all samples offset aligned single nt plot D30.png"), width = 1300, height = 300)
grid.arrange(UTR5_end_plot, CDS_start_plot, CDS_end_plot, UTR3_start_plot, nrow = 1, widths = c(1,2,2,1),
             top = textGrob("D30s",gp=gpar(fontsize=20,font=3)))
dev.off()
