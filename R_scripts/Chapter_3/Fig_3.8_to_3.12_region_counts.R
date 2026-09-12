#load libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(parallel)

#read in common variables
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/R_scripts/common_variables.R")

#set the read lengths you wish to plot
lengths <- 25:35

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
  
#read in data----
#generate a list of file names
fyle_list <- list()
for(sample in RPF_sample_names) {
  for(i in lengths){
    fyle_list[[paste(sample, i, sep = "_")]] <- file.path(parent_dir, "Analysis/region_counts", paste0(sample, "_pc_L", i, "_Off0_region_counts.csv"))
  }
}

#read in the data with parLapply
no_cores <- detectCores() - 1 #sets the number of cores to use (all but one)
cl <- makeCluster(no_cores) #Initiates cluster
data_list <- parLapply(cl, fyle_list, read_counts_csv) #reads in the data
stopCluster(cl) #Stops cluster

#combine data_list into one data frame
#extract sample and read length from fyle and inner_join with sample info to get condition and replicate
do.call("rbind", data_list) %>%
  mutate(read_length = str_remove(fyle, ".+pc_L"),
         read_length = as.numeric(str_remove(read_length, "_Off0_region_counts.csv")),
         sample = str_remove(fyle, ".+region_counts/"),
         sample = factor(str_remove(sample, "_pc.+"))) %>%
  inner_join(RPF_sample_info, by = "sample") %>%
  select(-fyle) -> all_data

summary(all_data)

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
            legend.text = element_text(size = 14),
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

#plot data----
#plot length distribution of all counts
#summed counts
all_data %>%
  mutate(condition = case_when(condition == "siRPS25" ~ "sieS25", .default = condition)) %>% 
  group_by(condition, replicate, read_length) %>%
  summarise(read_length_counts = sum(counts, na.rm = T)) %>%
  mutate(replicate = factor(paste0("REP", replicate), levels = c("REP1", "REP2", "REP3", "REP4", "REP5"))) %>% 
  ggplot(aes(x = read_length, y = read_length_counts, colour = condition))+
  facet_wrap(~replicate, nrow = 1) +
  scale_color_manual(values = c("#FFC20A", "#88CCEE", "#984EA3")) +
  geom_line(size = 1)+
  geom_point(size = 2)+
  publication_theme() +
  ylab("Total counts")+
  xlab("Read length")+
  ggtitle("RPF read length counts") +
  scale_x_continuous(breaks = seq(min(lengths), max(lengths),2)) +
  theme(panel.spacing.x = unit(1.5, "lines"))-> total_counts_plot

png(filename = file.path(parent_dir, "plots/summed_counts/lengths_count_plot.png"), width = 500, height = 300)
print(total_counts_plot)
dev.off()

png(filename = file.path("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/riboseq_read_length_counts.png"), width = 3500, height = 1250, res = 300)
print(total_counts_plot)
dev.off()


#percentage counts
all_data %>%
  group_by(condition, replicate) %>%
  summarise(total_counts = sum(counts, na.rm = T)) -> total_counts

all_data %>%
  group_by(condition, replicate, read_length) %>%
  summarise(read_length_counts = sum(counts, na.rm = T)) %>%
  inner_join(total_counts, by = c("condition", "replicate")) %>%
  mutate(perc_counts = (read_length_counts / total_counts) * 100) %>%
  ggplot(aes(x = read_length, y = perc_counts, colour = condition, lty = replicate, shape = replicate))+
  geom_line(size = 1)+
  geom_point(size = 2)+
  myTheme+
  ylab("% counts")+
  xlab("read length")+
  scale_x_continuous(breaks = seq(min(lengths), max(lengths),2)) -> perc_counts_plot

png(filename = file.path(parent_dir, "plots/summed_counts/lengths_perc_plot.png"), width = 500, height = 300)
print(perc_counts_plot)
dev.off()

#for each sample, plot summed counts within each region across the length distribution
for (sample in RPF_sample_names) {
  all_data[all_data$sample == sample,] %>%
    mutate(region = factor(region, levels = c("UTR3", "UTR5", "CDS"), labels = c("3\'UTR", "5\'UTR", "CDS"), ordered = T)) %>%
    group_by(read_length, region) %>%
    summarise(region_counts = sum(counts, na.rm = T)) %>%
    ggplot(aes(x = read_length, y = region_counts, fill = region))+
    geom_col()+
    publication_theme()+
    xlab("read length")+
    ylab("total counts")+
    scale_x_continuous(breaks = seq(min(lengths), max(lengths),2))+
    ggtitle(str_remove(str_replace_all(sample, "_", " "), "\\sM$")) -> region_counts_plot
  
  png(filename = file.path(parent_dir, paste0("plots/summed_counts/", sample, "_summed_region_counts.png")), width = 1000, height = 1000, res = 300)
  print(region_counts_plot)
  dev.off()
}

#for each sample, plot percentage counts within each region across the length distribution
all_data %>%
  group_by(sample, read_length) %>%
  summarise(total_counts = sum(counts, na.rm = T)) -> read_length_counts

for (sample in RPF_sample_names) {
  
  preprocessed_percentage_data <- 
  all_data[all_data$sample == sample,] %>%
    group_by(sample, region, read_length) %>%
    summarise(region_counts = sum(counts, na.rm = T)) %>%
    ungroup() %>%
    inner_join(read_length_counts, by = c("sample", "read_length")) %>%
    mutate(perc_counts = (region_counts / total_counts * 100),
           region = factor(region, levels = c("UTR3", "UTR5", "CDS"), labels = c("3\'UTR", "5\'UTR", "CDS"), ordered = T))
    #summary()
  
  nom <- str_split(preprocessed_percentage_data$sample, "_")
  first_term <- unique(unlist(lapply(nom, function(x){x[[1]]})))
  second_term <- unique(unlist(lapply(nom, function(x){x[[2]]})))
  
  if(second_term == "siRPS25"){
    
    second_term <- "sieS25"
  }
  
  gg_title <- paste(first_term, second_term)
  
  preprocessed_percentage_data %>% 
    ggplot(aes(x = read_length, y = perc_counts, fill = region))+
    geom_col()+
    publication_theme() +
    xlab("Read length")+
    ylab("Protein coding reads alignment (%)")+
    scale_x_continuous(breaks = seq(min(lengths), max(lengths),2))+
    ggtitle(gg_title) +
    theme(legend.position = "none",
          plot.margin=unit(c(2,2,2,2),"mm"))-> region_perc_plot
  
  png(filename = file.path(parent_dir, paste0("plots/summed_counts/", sample, "_perc_region_counts.png")), width = 1000, height = 1200, res = 300)
  print(region_perc_plot)
  dev.off()
}

all_data[all_data$sample == "REP1_NTC_M",] %>%
  group_by(sample, region, read_length) %>%
  summarise(region_counts = sum(counts, na.rm = T)) %>%
  ungroup() %>%
  inner_join(read_length_counts, by = c("sample", "read_length")) %>%
  mutate(perc_counts = (region_counts / total_counts * 100),
         region = factor(region, levels = c("UTR3", "UTR5", "CDS"), labels = c("3\'UTR", "5\'UTR", "CDS"), ordered = T)) %>%
  #summary()
  ggplot(aes(x = read_length, y = perc_counts, fill = region))+
  geom_col()+
  publication_theme() +
  xlab("read length")+
  ylab("% counts")+
  scale_x_continuous(breaks = seq(min(lengths), max(lengths),2))+
  ggtitle(str_remove(str_replace_all(sample, "_", " "), "\\sM$")) +
  theme(
        plot.margin=unit(c(2,2,2,2),"mm"))-> region_perc_plot_NTC_REP1

get_only_legend <- function(plot) { 
  
  # get tabular interpretation of plot 
  plot_table <- ggplot_gtable(ggplot_build(plot))  
  
  #  Mark only legend in plot 
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")  
  
  # extract legend 
  legend <- plot_table$grobs[[legend_plot]] 
  
  # return legend 
  return(legend)  
}

cowplot::ggdraw(cowplot::get_legend(region_perc_plot_NTC_REP1))

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/region_percentages_legend.png", res = 300, height = 100,width = 1500)
cowplot::ggdraw(cowplot::get_legend(region_perc_plot_NTC_REP1))
dev.off()

#plot region percentages for all samples in one plot
all_data %>%
  group_by(sample, region, read_length) %>%
  summarise(region_counts = sum(counts, na.rm = T)) %>%
  ungroup() %>%
  inner_join(read_length_counts, by = c("sample", "read_length")) %>%
  mutate(perc_counts = (region_counts / total_counts * 100),
         region = factor(region, levels = c("UTR3", "UTR5", "CDS"), labels = c("3\'UTR", "5\'UTR", "CDS"), ordered = T)) %>%
  ggplot(aes(x = read_length, y = perc_counts, colour = region))+
  geom_point()+
  #geom_boxplot(aes(fill = region))+
  myTheme+
  xlab("read length")+
  ylab("% counts")+
  scale_x_continuous(breaks = seq(min(lengths), max(lengths),2)) -> region_perc_plot

png(filename = file.path(parent_dir, "plots/summed_counts/all_samples_perc_region_counts.png"), width = 500, height = 300)
print(region_perc_plot)
dev.off()
