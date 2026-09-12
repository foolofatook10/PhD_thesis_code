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

#themes
myTheme <- theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5))

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


#read in data----
#generate a list of file names
fyle_list <- list()
for(sample in RPF_sample_names) {
  for(i in lengths){
    fyle_list[[paste(sample, i, "start", sep = "_")]] <- file.path(parent_dir, "Analysis/spliced_counts", paste0(sample, "_pc_L", i, "_Off0_start_site.csv"))
    fyle_list[[paste(sample, i, "stop", sep = "_")]] <- file.path(parent_dir, "Analysis/spliced_counts", paste0(sample, "_pc_L", i, "_Off0_stop_site.csv"))
  }
}

#read in the data with parLapply
no_cores <- detectCores() - 1 #sets the number of cores to use (all but one)
cl <- makeCluster(no_cores) #Initiates cluster
data_list <- parLapply(cl, fyle_list, read_counts_csv) #reads in the data
stopCluster(cl) #Stops cluster

#combine data_list into one data frame
#extract sample, read length and splice position from fyle
do.call("rbind", data_list) %>%
  mutate(read_length = str_remove(fyle, ".+pc_L"),
         read_length = as.numeric(str_remove(read_length, "_Off0_.+")),
         sample = str_remove(fyle, ".+spliced_counts/"),
         sample = factor(str_remove(sample, "_pc.+")),
         splice = str_remove(fyle, ".+Off0_"),
         splice = factor(str_remove(splice, ".csv"))) %>%
  select(-fyle) -> all_data

summary(all_data)

#plot offset plots----
for (sample in RPF_sample_names) {
  
  start_site_plot_list <- list()
  stop_site_plot_list <- list()
  
  for (i in lengths) {
    start_site_data <- all_data[all_data$sample == sample & all_data$splice == "start_site" & all_data$read_length == i,]
    stop_site_data <- all_data[all_data$sample == sample & all_data$splice == "stop_site" & all_data$read_length == i,]
    
    start_site_data %>%
      ggplot(aes(x = position, y = counts)) + 
      geom_col()+
      xlab("Position relative to start codon")+
      ylab("Total counts")+
      geom_vline(xintercept = -12, colour = "red", lty=2)+
      scale_x_continuous(limits = c(-25, 25), breaks = c(-25, -12, 0, 25))+
      publication_theme()+
      ggtitle(paste("read length", i)) -> start_site_plot_list[[i]]
    
    stop_site_data %>%
      ggplot(aes(x = position, y = counts)) + 
      geom_col()+
      xlab("Position relative to stop codon")+
      ylab("Total counts")+
      geom_vline(xintercept = -18, colour = "red", lty=2)+
      scale_x_continuous(limits = c(-25, 25), breaks = c(-25, -18, 0, 25))+
      publication_theme()+
      ggtitle(paste("read length", i)) -> stop_site_plot_list[[i]]
  }
  png(filename = file.path(parent_dir, paste0("plots/offset/", sample, "_start_site_offset_21nts.png")), width = 1000, height = 500)
  grid.arrange(start_site_plot_list[[18]], start_site_plot_list[[19]], start_site_plot_list[[20]],
               start_site_plot_list[[21]], start_site_plot_list[[22]], start_site_plot_list[[23]],
               start_site_plot_list[[24]], start_site_plot_list[[25]], start_site_plot_list[[26]], nrow = 3)
  dev.off()
  
  png(filename = file.path(parent_dir, paste0("plots/offset/", sample, "_stop_site_offset_21nts.png")), width = 1000, height = 500)
  grid.arrange(stop_site_plot_list[[18]], stop_site_plot_list[[19]], stop_site_plot_list[[20]],
               stop_site_plot_list[[21]], stop_site_plot_list[[22]], stop_site_plot_list[[23]],
               stop_site_plot_list[[24]], stop_site_plot_list[[25]], stop_site_plot_list[[26]], nrow = 3)
  dev.off()
}

#all samples
start_site_plot_list_M <- list()
stop_site_plot_list_M <- list()

all_data$MD30 <- str_remove_all(str_extract(all_data$sample, "_.*_"), "_")

for (i in lengths) {
  
  all_data2 <- all_data %>% filter(MD30 == "M")
  
  start_site_data <- all_data2[all_data2$splice == "start_site" & all_data2$read_length == i,]
  stop_site_data <- all_data2[all_data2$splice == "stop_site" & all_data2$read_length == i,]
  
  start_site_data %>%
    ggplot(aes(x = position, y = counts)) + 
    geom_col()+
    xlab("Position relative to start codon")+
    ylab("Total counts")+
    geom_vline(xintercept = -12, colour = "red", lty=2)+
    scale_x_continuous(limits = c(-25, 25), breaks = c(-25, -12, 0, 25))+
    publication_theme()+
    ggtitle(paste("read length", i)) -> start_site_plot_list_M[[i]]
  
  stop_site_data %>%
    ggplot(aes(x = position, y = counts)) + 
    geom_col()+
    xlab("Position relative to stop codon")+
    ylab("Total counts")+
    geom_vline(xintercept = -18, colour = "red", lty=2)+
    scale_x_continuous(limits = c(-25, 25), breaks = c(-25, -18, 0, 25))+
    publication_theme()+
    ggtitle(paste("read length", i)) -> stop_site_plot_list_M[[i]]
}
png(filename = file.path(parent_dir, paste0("plots/offset/all_samples_start_site_offset_M_21nts.png")), width = 1000, height = 500)
grid.arrange(start_site_plot_list_M[[18]], start_site_plot_list_M[[19]], start_site_plot_list_M[[20]],
             start_site_plot_list_M[[21]], start_site_plot_list_M[[22]], start_site_plot_list_M[[23]],
             start_site_plot_list_M[[24]], start_site_plot_list_M[[25]], start_site_plot_list_M[[26]], nrow = 3)
dev.off()

png(filename = file.path(parent_dir, paste0("plots/offset/all_samples_stop_site_offset_M_21nts.png")), width = 1000, height = 500)
grid.arrange(stop_site_plot_list_M[[18]], stop_site_plot_list_M[[19]], stop_site_plot_list_M[[20]],
             stop_site_plot_list_M[[21]], stop_site_plot_list_M[[22]], stop_site_plot_list_M[[23]],
             stop_site_plot_list_M[[24]], stop_site_plot_list_M[[25]], stop_site_plot_list_M[[26]], nrow = 3)
dev.off()



start_site_plot_list_D30 <- list()
stop_site_plot_list_D30 <- list()


for (i in lengths) {
  
  all_data3 <- all_data %>% filter(MD30 == "D30")
  
  start_site_data <- all_data3[all_data3$splice == "start_site" & all_data3$read_length == i,]
  stop_site_data <- all_data3[all_data3$splice == "stop_site" & all_data3$read_length == i,]
  
  start_site_data %>%
    ggplot(aes(x = position, y = counts)) + 
    geom_col()+
    xlab("Position relative to start codon")+
    ylab("Total counts")+
    geom_vline(xintercept = -12, colour = "red", lty=2)+
    scale_x_continuous(limits = c(-25, 25), breaks = c(-25, -12, 0, 25))+
    publication_theme()+
    ggtitle(paste("read length", i)) -> start_site_plot_list_D30[[i]]
  
  stop_site_data %>%
    ggplot(aes(x = position, y = counts)) + 
    geom_col()+
    xlab("Position relative to stop codon")+
    ylab("Total counts")+
    geom_vline(xintercept = -18, colour = "red", lty=2)+
    scale_x_continuous(limits = c(-25, 25), breaks = c(-25, -18, 0, 25))+
    publication_theme()+
    ggtitle(paste("read length", i)) -> stop_site_plot_list_D30[[i]]
}
png(filename = file.path(parent_dir, paste0("plots/offset/all_samples_start_site_offset_D30_21nts.png")), width = 1000, height = 500)
grid.arrange(start_site_plot_list_D30[[18]], start_site_plot_list_D30[[19]], start_site_plot_list_D30[[20]],
             start_site_plot_list_D30[[21]], start_site_plot_list_D30[[22]], start_site_plot_list_D30[[23]],
             start_site_plot_list_D30[[24]], start_site_plot_list_D30[[25]], start_site_plot_list_D30[[26]], nrow = 3)
dev.off()

png(filename = file.path(parent_dir, paste0("plots/offset/all_samples_stop_site_offset_D30_21nts.png")), width = 1000, height = 500)
grid.arrange(stop_site_plot_list_D30[[18]], stop_site_plot_list_D30[[19]], stop_site_plot_list_D30[[20]],
             stop_site_plot_list_D30[[21]], stop_site_plot_list_D30[[22]], stop_site_plot_list_D30[[23]],
             stop_site_plot_list_D30[[24]], stop_site_plot_list_D30[[25]], stop_site_plot_list_D30[[26]], nrow = 3)
dev.off()

