###This script was written by Joe and will make plots of individual traces and all traces overlaid.
###Within the parent directory (set below) the traces need to be in a directory named "traces" and there needs to be an empty directory called "plots" to save these too.

#Imports
library(tidyverse)
library(RColorBrewer)

#set parent_dir as the folder with the csv file in----
parent_dir <- paste0("\\\\data.beatson.gla.ac.uk/data/JETTLES/CNOT3_Riboseq/gradients_final")

#set the lower and upper limit for the y axis
ylims <- c(0, 1.75)

#make a list of files to plot data from
setwd(file.path(parent_dir, "traces"))
fylenames <- dir(pattern = "*.csv") # creates the list of all the csv files in the directory

#use a for loop to read in the data, reformat it, plot the trace a save the reformatted data to a list
data_list <- list() #creates a list to save the data in

replicate <- as.numeric(str_remove_all(unlist(lapply(str_split(fylenames, "_"), function(a){a[[5]]})), "Riboseq"))
condition <- unlist(lapply(str_split(fylenames, "_"), function(a){a[[6]]}))
digestion <- str_remove_all(unlist(lapply(str_split(fylenames, "_"), function(a){a[[7]]})), ".csv")

process_data <- function(x){
  
  df <- read_csv(x, skip = 45)

  #create a continuous scale of cumulative volume to plot against absorbance
  total_vol <- sum(as.numeric(df$`Fraction Volume(ml)`), na.rm = T)
  step <- total_vol / nrow(df)
  df$cum_vol <- seq(step, total_vol, step)
  
  #extract sample name
  sample_name <- str_remove(x, "\\.csv") #removes the .csv
  
  replicate <- unlist(lapply(str_split(sample_name, "_"), function(f){f[[5]]}))
  replicate <- as.numeric(str_remove(replicate, "Riboseq"))
  
  condition <- unlist(lapply(str_split(sample_name, "_"), function(f){f[[6]]}))
  
  digestion <- str_remove(unlist(lapply(str_split(sample_name, "_"), function(f){f[[7]]})), ".csv")
  
  df$replicate <- rep(replicate)
  df$condition <- condition
  df$digestion <- digestion
  
  return(df)
  
  
}


undigested_traces <- 
do.call("rbind", lapply(fylenames, process_data)) %>% 
  filter(digestion == "undigested")

undigested_split <- split(undigested_traces, undigested_traces$replicate)

publication_theme <- function(base_size=14, base_family="Helvetica") {
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
            axis.title.y = element_text(angle=90,vjust =2, size = 14),
            axis.title.x = element_text(vjust = -0.2, size = 14),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


plot_CNOT3_traces <- 
  function(df){
    
    total_vol <- sum(as.numeric(df$`Fraction Volume(ml)`), na.rm = T)
    repli <- unique(df$replicate)
    
    print(df)
    
    df2 <- df %>%
      filter(condition %in% c("CNOT3", "NTC"))
    
    df2$condition <- factor(df2$condition,
                           levels = c("NTC", "CNOT3"),
                           labels = c("NTC", "siCNOT3"))
    
    
    df2 %>%
      #filter(condition %in% c("CNOT3", "NTC")) %>%
      ggplot(aes(x = cum_vol, y = Absorbance, colour = condition))+
      geom_line(size = 1.5)+
      scale_colour_brewer(palette = "Dark2")+
      theme_bw()+
      scale_x_continuous(limits = c(2, 12), breaks = 2:12)+
      scale_color_manual(values = c("#FFC20A", "#88CCEE")) +
      ylim(ylims)+
      ylab("UV Absorbance (254nm)")+
      xlab("Fraction")+
      #ggtitle(paste0("Replicate: ", repli)) +
      publication_theme() +
      theme(axis.title = element_text(size = 20),
            axis.text = element_text(size = 18),
            axis.line = element_line(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.text = element_text(size = 20),
            legend.title = element_blank())
    
  }


plot_RPS25_traces <- 
  function(df){
    
    total_vol <- sum(as.numeric(df$`Fraction Volume(ml)`), na.rm = T)
    repli <- unique(df$replicate)
    
    print(df)
    
    df2 <- df %>%
      filter(condition %in% c("RPS25", "NTC"))
    
    df2$condition <- factor(df2$condition,
                            levels = c("NTC", "RPS25"),
                            labels = c("NTC", "siRPS25"))
    
    df2 %>%
      #filter(condition %in% c("RPS25", "NTC")) %>%
      ggplot(aes(x = cum_vol, y = Absorbance, colour = condition))+
      geom_line(size = 1.5)+
      scale_colour_brewer(palette = "Dark2")+
      theme_bw()+
      scale_x_continuous(limits = c(2, 12), breaks = 2:12)+
      scale_color_manual(values = c("#FFC20A", "#984EA3")) +
      ylim(ylims)+
      ylab("UV Absorbance (254nm)")+
      xlab("Fraction")+
      #ggtitle(paste0("Replicate: ", repli)) +
      publication_theme() +
      theme(axis.title = element_text(size = 20),
            axis.text = element_text(size = 18),
            axis.line = element_line(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.text = element_text(size = 20),
            legend.title = element_blank())
    
  }



CNOT3_undigested_traces <- lapply(undigested_split, plot_CNOT3_traces)
RPS25_undigested_traces <- lapply(undigested_split, plot_RPS25_traces)

fyle_names <- paste0("CNOT3_undigested_REP", 1:5, ".png")
fyle_names_RPS25 <- paste0("RPS25_undigested_REP", 1:5, ".png")

save_plots <- function(x){

  
  png(file.path("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures", fyle_names[[x]]),
      res = 300, height = 1000, width = 1500)
  print(CNOT3_undigested_traces[[x]])
  dev.off()
  
}

save_plots_RPS25 <- function(x){
  
  
  png(file.path("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures", fyle_names_RPS25[[x]]),
      res = 300, height = 1000, width = 1500)
  print(RPS25_undigested_traces[[x]])
  dev.off()
  
}


lapply(1:length(CNOT3_undigested_traces), save_plots)
lapply(1:length(RPS25_undigested_traces), save_plots_RPS25)

all_data <- do.call("rbind", data_list)

all_data$fraction_number <- max(1:max(all_data$cum_vol))

test <- all_data[657:nrow(all_data),]

CNOT3_data <- all_data %>% filter(sample %in% c("NTC undigested", "CNOT3 undigested"))

CNOT3_data$sample <- factor(CNOT3_data$sample, levels = c("NTC undigested", "CNOT3 undigested"))

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
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

trace_plot <- ggplot(CNOT3_data, aes(x = cum_vol, y = Absorbance, colour = sample))+
  geom_line(size = 1)+
  scale_colour_brewer(palette = "Dark2")+ 
  theme_bw()+
  scale_x_continuous(limits = c(2, total_vol), breaks = 2:total_vol)+
  scale_color_manual(values = c()) +
  ylim(ylims)+
  ylab("absorbance")+
  xlab("volume")+
  publication_theme() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        axis.line = element_line(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_blank())#,
#legend.position = "bottom")

RPS25_data <- all_data %>% filter(sample %in% c("NTC undigested", "RPS25 undigested"))

RPS25_data$sample <- factor(RPS25_data$sample, levels = c("NTC undigested", "RPS25 undigested"))



trace_plot <- ggplot(RPS25_data, aes(x = cum_vol, y = Absorbance, colour = sample))+
  geom_line(size = 1)+
  scale_colour_brewer(palette = "Dark2")+ 
  theme_bw()+
  scale_x_continuous(limits = c(2, total_vol), breaks = 2:total_vol)+
  scale_color_manual(values = c("#E41A1C", "#984EA3")) +
  ylim(ylims)+
  ylab("absorbance")+
  xlab("volume")+
  publication_theme() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        axis.line = element_line(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_blank())#,
#legend.position = "bottom")

png(file = file.path(parent_dir, "plots/overlaid_traces.png"), height = 400, width = 1000)
print(trace_plot)
dev.off()

digested_traces <- 
  do.call("rbind", lapply(fylenames, process_data)) %>% 
  filter(digestion == "digested")


digested_split <- split(digested_traces, digested_traces$replicate)

lapply(digested_split, plot_CNOT3_traces)
