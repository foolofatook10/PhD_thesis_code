#load libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(parallel)

#read in common variables
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/R_scripts/common_variables.R")
lengths <- 16:37

#functions
#write a function that will read in a csv file for use with parLapply
read_counts_csv <- function(k){
  df <- read.csv(file = k)
  df$fyle <- rep(k)
  return(df)
}

#write theme
myTheme <- theme_bw()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.title = element_blank())

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
            axis.title = element_text(face = "bold",size = 20),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(size = 18, color = "black"),
            axis.text.y = element_text(size = 18, color = "black"),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = 18),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(1, "cm"),
            legend.margin = unit(0.5, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

#read in data----
#generate a list of file names
fyle_list <- list()
for(sample in RPF_sample_names) {
  for(i in lengths){
    fyle_list[[paste(sample, i, sep = "_")]] <- file.path(parent_dir, "Analysis/periodicity", paste0(sample, "_pc_L", i, "_Off0_periodicity.csv"))
  }
}

#read in the data with parLapply
no_cores <- detectCores() - 1 #sets the number of cores to use (all but one)
cl <- makeCluster(no_cores) #Initiates cluster
data_list <- parLapply(cl, fyle_list, read_counts_csv) #reads in the data
stopCluster(cl) #Stops cluster

#combine data_list into one data frame
#extract sample and read length from fyle
do.call("rbind", data_list) %>%
  mutate(read_length = str_remove(fyle, ".+pc_L"),
         read_length = as.numeric(str_remove(read_length, "_Off0_.+")),
         sample = str_remove(fyle, ".+periodicity/"),
         sample = factor(str_remove(sample, "_pc.+"))) %>%
  select(-fyle) -> all_data

summary(all_data)

#plot data----
#individual samples
for (sample in RPF_sample_names) {
  
  sample_split <- str_split(sample, "_")
  
  if(sample_split[[1]][[3]] == "CNOT3"){
    
    element_two = "CNOT3 IP'd RPFs"
    
  } else{
    
    element_two = "Total RPFs"
    
  }
  
  element_one <- paste0("REP", (as.numeric(str_remove(sample_split[[1]][[1]], "REP")) - 1))
  
  tytle <- paste(element_one, element_two)
  
  
  df <- all_data[all_data$sample == sample,]
  
  df %>%
    group_by(read_length) %>%
    summarise(total_counts = sum(counts)) -> summed_counts
  
  df_temp <- 
  df %>%
    group_by(read_length, frame) %>%
    summarise(frame_counts = sum(counts)) %>%
    inner_join(summed_counts, by = "read_length") %>%
    mutate(frame_perc = (frame_counts / total_counts) * 100,
           frame = factor(frame, levels = c("f2", "f1", "f0"), ordered = T)) 
  
  #print(df_temp)
  
  df_temp %>% 
    ggplot(aes(x = read_length, y = frame_perc, fill = frame))+
    geom_col()+
    scale_x_continuous(breaks = seq(min(lengths), max(lengths),2))+
    ylab("% Counts")+
    xlab("Read length")+
    ggtitle(tytle)+
    publication_theme() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal")-> periodicity_col_plot
  
  leg <- cowplot::get_legend(periodicity_col_plot)
  
  png(filename = file.path(parent_dir, paste0("plots/periodicity/", sample, "_periodicity_leg.png")), res = 300, width = 1000, height = 200)
  print(cowplot::ggdraw(leg))
  dev.off()
  
}
  
  plot_list <- list()
  for (i in lengths) {
    df[df$read_length == i,] %>%
      group_by(read_length, transcript) %>%
      summarise(total_counts = sum(counts)) %>%
      filter(total_counts > 50) -> transcript_counts
    
    df[df$read_length == i,] %>%
      inner_join(transcript_counts, by = c("transcript", "read_length")) %>%
      mutate(frame_perc = (counts / total_counts) * 100) %>%
      ggplot(aes(x = frame, y = frame_perc, fill = frame))+
      geom_boxplot()+
      xlab("Frame")+
      ylab("% Counts")+
      publication_theme() +
      theme(legend.position = "none")+
      ggtitle(paste("read length", i)) -> plot_list[[i]]
  }
  # png(filename = file.path(parent_dir, paste0("plots/periodicity/", sample, "_periodicity_boxplots.png")), width = 1000, height = 750)
  # grid.arrange(plot_list[[28]], plot_list[[29]], plot_list[[30]],
  #              plot_list[[31]], plot_list[[32]], plot_list[[33]],
  #              plot_list[[34]], plot_list[[35]], plot_list[[36]], nrow = 3)
  # dev.off()
}

#all samples
all_data$MD30 <- 
  str_remove_all(
    str_extract(all_data$sample, pattern = "_.*_"), "_")

all_data$MD30 <- factor(all_data$MD30, levels = c("M", "D30"))

all_data %>%
  group_by(MD30, read_length, sample) %>%
  summarise(total_counts = sum(counts)) -> summed_counts



all_data %>%
  group_by(MD30, read_length, frame, sample) %>%
  summarise(frame_counts = sum(counts)) %>%
  inner_join(summed_counts, by = c("MD30", "read_length", "sample")) %>%
  mutate(frame_perc = (frame_counts / total_counts) * 100,
         frame = factor(frame, levels = c("f2", "f1", "f0"), ordered = T)) %>%
  ggplot(aes(x = factor(read_length), y = frame_perc, fill = frame))+
  geom_boxplot(width = 0.5, outlier.shape=NA)+
  facet_wrap(~MD30, nrow = 2) +
  ylab("% counts")+
  xlab("read length")+
  myTheme -> periodicity_col_plot

png(filename = file.path(parent_dir, paste0("plots/periodicity/all_samples_periodicity.png")), width = 600, height = 500)
print(periodicity_col_plot)
dev.off()
