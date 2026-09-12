#load libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(parallel)

#read in common variables
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/R_scripts/common_variables.R")

#set the read lengths you wish to plot
lengths <- 16:37

#functions
#write a function that will read in a csv file for use with parLapply
read_counts_csv <- function(k){
  df <- read.csv(file = k)
  df$fyle <- rep(k)
  return(df)
}

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
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
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

#plot data----
#plot length distribution of all counts
#summed counts
read_lengths <- 
  all_data %>%
  group_by(condition2, condition, replicate, read_length) %>%
  summarise(read_length_counts = sum(counts, na.rm = T))

read_lengths$condition2 <- factor(read_lengths$condition2,
                                  levels = c("M", "D30"))

read_lengths <- 
read_lengths %>% 
  mutate(replicate = factor(case_when(replicate == "2" ~"REP1",
                               replicate == "3" ~"REP2",
                               replicate == "4" ~"REP3"),
                            levels = c("REP1", "REP2", "REP3")))

read_lengths %>% 
  filter(condition2 == "M") %>% 
  ggplot(aes(x = read_length, y = read_length_counts, colour = condition))+
  geom_line(size = 1)+
  geom_point(size = 2)+
  facet_wrap(~replicate, drop = T, nrow = 3, ncol = 1) +
  scale_color_manual(values = c("#0C7BDC", "#FFC20A")) +
  geom_vline(xintercept = 28, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 36, linetype = "dashed", size = 1) +
  publication_theme()+
  ylab("Raw counts")+
  xlab("read length")+
  ggtitle("RPF Length distribution") +
  scale_x_continuous(breaks = seq(min(lengths), max(lengths),2)) -> total_counts_plot

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3SelRP_read_length_distribution_raw.png",
    res = 300, height = 2000, width = 1250)
print(total_counts_plot)
dev.off()

# gplot_mod <- function(gplot){
# g <- ggplotGrob(gplot)
# # get the grobs that must be removed
# rm_grobs <- g$layout$name %in% c("panel-1-2")
# # remove grobs
# g$grobs[rm_grobs] <- NULL
# g$layout <- g$layout[!rm_grobs, ]
# return(g)
# }
# 
# total_counts_plot_final <- gplot_mod(total_counts_plot)
# 
# grid.newpage()
# grid.draw(total_counts_plot_final)


# png(filename = file.path(parent_dir, "plots/summed_counts/lengths_count_plot.png"), width = 700, height = 500)
# grid.newpage()
# grid.draw(total_counts_plot_final)
# dev.off()

#percentage counts
all_data %>%
  group_by(condition2, condition, replicate) %>%
  summarise(total_counts = sum(counts, na.rm = T)) -> total_counts

percentage_counts <-
  all_data %>%
  group_by(condition2, condition, replicate, read_length) %>%
  summarise(read_length_counts = sum(counts, na.rm = T)) %>%
  inner_join(total_counts, by = c("condition2", "condition", "replicate")) %>%
  mutate(perc_counts = (read_length_counts / total_counts) * 100) 

percentage_counts$condition2 <- factor(percentage_counts$condition2, levels = c("M","D30"))

percentage_counts %>% 
  mutate(replicate = factor(case_when(replicate == "2" ~"REP1",
                                      replicate == "3" ~"REP2",
                                      replicate == "4" ~"REP3"),
                            levels = c("REP1", "REP2", "REP3"))) %>% 
  filter(condition2 != "D30") %>% 
  ggplot(aes(x = read_length, y = perc_counts, colour = condition))+
  facet_wrap(~replicate, drop = T, nrow = 3) +
  geom_line(size = 1)+
  geom_point(size = 2)+
  scale_color_manual(values = c("#0C7BDC", "#FFC20A")) +
  geom_vline(xintercept = 28, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 36, linetype = "dashed", size = 1) +
  ggtitle("RPF Length distribution") +
  publication_theme()+
  ylab("% counts")+
  xlab("read length")+
  scale_x_continuous(breaks = seq(min(lengths), max(lengths),2)) -> perc_counts_plot

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3SelRP_read_length_distribution_percentage.png",
    res = 300, height = 2000, width = 1250)
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
    myTheme+
    xlab("read length")+
    ylab("total counts")+
    scale_x_continuous(breaks = seq(min(lengths), max(lengths),2))+
    ggtitle(str_replace_all(sample, "_", " ")) -> region_counts_plot
  
  png(filename = file.path(parent_dir, paste0("plots/summed_counts/", sample, "_summed_region_counts.png")), width = 500, height = 300)
  print(region_counts_plot)
  dev.off()
}

#for each sample, plot percentage counts within each region across the length distribution

all_data <- 
all_data %>%
  mutate(sample = case_when(sample == "REP2_M_Tot" ~ "REP1 Tot",
                             sample == "REP3_M_Tot" ~ "REP2 Tot",
                             sample == "REP4_M_Tot" ~ "REP3 Tot",
                             sample == "REP2_M_CNOT3" ~ "REP1 CNOT3",
                             sample == "REP3_M_CNOT3" ~ "REP2 CNOT3",
                             sample == "REP4_M_CNOT3" ~ "REP3 CNOT3",
                            .default = sample))  

all_data$sample %>% unique()

all_data %>%   
group_by(condition2, sample, read_length) %>%
  summarise(total_counts = sum(counts, na.rm = T)) -> read_length_counts


new_RPF_sample_names <- unique(all_data$sample)

sample = "REP1 CNOT3"

for (sample in new_RPF_sample_names) {
  all_data[all_data$sample == sample,] %>%
    group_by(condition2, sample, region, read_length) %>%
    summarise(region_counts = sum(counts, na.rm = T)) %>%
    ungroup() %>%
    inner_join(read_length_counts, by = c("condition2", "sample", "read_length")) %>%
    mutate(perc_counts = (region_counts / total_counts * 100),
           region = factor(region, levels = c("UTR3", "UTR5", "CDS"), labels = c("3\'UTR", "5\'UTR", "CDS"), ordered = T)) %>%
    #summary()
    ggplot(aes(x = read_length, y = perc_counts, fill = region))+
    geom_col()+
    publication_theme()+
    xlab("read length")+
    ylab("% counts")+
    scale_x_continuous(breaks = seq(min(lengths), max(lengths),2))+
    ggtitle(str_replace_all(sample, "_", " ")) +
    theme(legend.position = "none")-> region_perc_plot
  
  png(filename = file.path(parent_dir, paste0("plots/summed_counts/", sample, "_perc_region_counts.png")), width = 1500, height = 1000, res = 300)
  print(region_perc_plot)
  dev.off()
}

#plot region percentages for all samples in one plot
all_data %>%
  group_by(condition2, sample, region, read_length) %>%
  summarise(region_counts = sum(counts, na.rm = T)) %>%
  ungroup() %>%
  inner_join(read_length_counts, by = c("condition2", "sample", "read_length")) %>%
  mutate(perc_counts = (region_counts / total_counts * 100),
         region = factor(region, levels = c("UTR3", "UTR5", "CDS"), labels = c("3\'UTR", "5\'UTR", "CDS"), ordered = T)) %>%
  ggplot(aes(x = read_length, y = perc_counts, colour = region))+
  geom_point()+
  #geom_line()+
  #geom_boxplot(aes(fill = region))+
  myTheme+
  xlab("read length")+
  ylab("% counts")+
  scale_x_continuous(breaks = seq(min(lengths), max(lengths),2)) -> region_perc_plot

png(filename = file.path(parent_dir, "plots/summed_counts/all_samples_perc_region_counts.png"), width = 500, height = 300)
print(region_perc_plot)
dev.off()
