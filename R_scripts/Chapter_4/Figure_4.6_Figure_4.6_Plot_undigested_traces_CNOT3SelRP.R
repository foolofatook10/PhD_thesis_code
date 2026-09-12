###This script was written by Joe and will make plots of individual traces and all traces overlaid.
###Within the parent directory (set below) the traces need to be in a directory named "traces" and there needs to be an empty directory called "plots" to save these too.

#Imports
library(tidyverse)
library(RColorBrewer)

# publication theme
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

#set parent_dir as the folder with the csv file in----
parent_dir <- paste0("\\\\data.beatson.gla.ac.uk/data/JETTLES/CNOT3_Sel_RibProf/sucrose_gradients")

#set the lower and upper limit for the y axis
ylims <- c(0, 4)

#make a list of files to plot data from
setwd(file.path(parent_dir, "traces"))
fylenames <- dir(pattern = "*.csv") # creates the list of all the csv files in the directory

#use a for loop to read in the data, reformat it, plot the trace a save the reformatted data to a list
data_list <- list() #creates a list to save the data in
for (fylename in fylenames){
  #read in data
  df <- read_csv(file = fylename, col_names = T, skip = 45) #skip the first 45 rows which contain information about the run settings
  
  #create a continuous scale of cumulative volume to plot against absorbance
  total_vol <- sum(as.numeric(df$`Fraction Volume(ml)`), na.rm = T)
  step <- total_vol / nrow(df)
  df$cum_vol <- seq(step, total_vol, step)
  
  #extract sample name
  sample_name <- str_remove(fylename, "\\.csv") #removes the .csv
  sample_name <- str_replace_all(sample_name, "\\_", " ") #replaces underscores with spaces
  
  df$sample <- rep(sample_name)
  
  #plot data
  trace_plot <- ggplot(data = df, aes(x = cum_vol, y = Absorbance))+
    geom_line(size = 1, color = "red")+
    theme_bw()+
    scale_x_continuous(limits = c(0, total_vol), breaks = 1:total_vol)+
    ylim(ylims)+
    ylab("absorbance")+
    xlab("volume")+
    ggtitle(sample_name)+
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 18),
          axis.line = element_line(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
          legend.text = element_text(size = 20),
          legend.title = element_blank())
  
  png(file = file.path(parent_dir, "plots", str_replace(fylename, ".csv", ".png")), height = 400, width = 800)
  print(trace_plot)
  dev.off()
  
  #save data to list
  data_list[[fylename]] <- df
}

all_data <- do.call("rbind", data_list)

all_data$rep <- unlist(lapply(str_split(all_data$sample, " "), function(x){x[[2]]}))
all_data$condition <- unlist(lapply(str_split(all_data$sample, " "), function(x){x[[6]]}))

undigesteds <- all_data %>% filter(condition == "undigested0")
undigesteds_split <- split(undigesteds, f = undigesteds$rep)

plot_data <- function(x){
  
  x %>% ggplot(aes(x = cum_vol, y = Absorbance))+
    geom_line(size = 1)+
    scale_colour_brewer(palette = "Dark2")+ 
    theme_bw()+
    scale_x_continuous(limits = c(0, total_vol), breaks = 1:total_vol)+
    ylim(ylims)+
    ylab("absorbance")+
    xlab("volume")+
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 18),
          axis.line = element_line(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          legend.text = element_text(size = 20),
          legend.title = element_blank())#,
  #legend.position = "bottom")
}

x_lab = "Volume (ml)"
y_lab = "UV Absorbance (254nm)"

replicate1_undigested <- 
undigesteds_split$REPLICATE2 %>% 
  filter(cum_vol > 1) %>% 
ggplot(aes(x = cum_vol, y = Absorbance))+
  geom_line(size = 1)+
  scale_colour_brewer(palette = "Dark2")+ 
  theme_bw()+
  scale_x_continuous(limits = c(1, total_vol), breaks = 1:total_vol)+
  ylim(0, 4.2)+
  ylab(y_lab)+
  xlab(x_lab)+
  ggtitle("Replicate 1") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        axis.line = element_line(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_blank()) +
  publication_theme()

replicate2_undigested <- 
undigesteds_split$REPLICATE3 %>% 
  filter(cum_vol > 1) %>% 
  ggplot(aes(x = cum_vol, y = Absorbance))+
  geom_line(size = 1)+
  scale_colour_brewer(palette = "Dark2")+ 
  theme_bw()+
  scale_x_continuous(limits = c(1, total_vol), breaks = 1:total_vol)+
  ylim(0, 3.5)+
  ylab(y_lab)+
  xlab(x_lab)+
  ggtitle("Replicate 2") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        axis.line = element_line(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_blank()) +
  publication_theme()

replicate3_undigested <- 
undigesteds_split$REPLICATE4 %>% 
  filter(cum_vol > 1) %>% 
  ggplot(aes(x = cum_vol, y = Absorbance))+
  geom_line(size = 1)+
  scale_colour_brewer(palette = "Dark2")+ 
  theme_bw()+
  scale_x_continuous(limits = c(1, total_vol), breaks = 1:total_vol)+
  ylim(0, 4.1)+
  ylab(y_lab)+
  xlab(x_lab)+
  ggtitle("Replicate 3") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        axis.line = element_line(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_blank()) +
  publication_theme()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/replicate1_undigested.png",
    res = 300, height = 1000, width = 1500)
print(replicate1_undigested)
dev.off()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/replicate2_undigested.png",
    res = 300, height = 1000, width = 1500)
print(replicate2_undigested)
dev.off()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/replicate3_undigested.png",
    res = 300, height = 1000, width = 1500)
print(replicate3_undigested)
dev.off()

trace_plot <- ggplot(data = all_data, aes(x = cum_vol, y = Absorbance, colour = sample))+
  geom_line(size = 1)+
  scale_colour_brewer(palette = "Dark2")+ 
  theme_bw()+
  scale_x_continuous(limits = c(0, total_vol), breaks = 1:total_vol)+
  ylim(ylims)+
  ylab("absorbance")+
  xlab("volume")+
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


