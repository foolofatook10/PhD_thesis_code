###This script was written by Joe and will make plots of individual traces and all traces overlaid.
###Within the parent directory (set below) the traces need to be in a directory named "traces" and there needs to be an empty directory called "plots" to save these too.

#Imports
library(tidyverse)
library(RColorBrewer)




#set parent_dir as the folder with the csv file in----
parent_dir <- "\\\\data.beatson.gla.ac.uk/data/JETTLES/CNOT3_Sel_RibProf/optimisations/crosslinking_gradients"

#set the lower and upper limit for the y axis
ylims <- c(0, 2)

#make a list of files to plot data from
setwd(file.path(parent_dir, "traces"))
fylenames <- dir(pattern = "*.csv") # creates the list of all the csv files in the directory

#use a for loop to read in the data, reformat it, plot the trace a save the reformatted data to a list
data_list <- list() #creates a list to save the data in

fylename = fylenames[[1]]

for (fylename in fylenames){
  #read in data
  df <- read_csv(file = fylename, col_names = T, skip = 45) #skip the first 45 rows which contain information about the run settings
  
  #create a continuous scale of cumulative volume to plot against absorbance
  total_vol <- sum(as.numeric(df$`Fraction Volume(ml)`), na.rm = T)
  step <- total_vol / nrow(df)
  df$cum_vol <- seq(step, total_vol, step)
  
  #extract sample name
  sample_name <- str_remove(fylename, "\\.csv") #removes the .csv
  sample_info <- str_split(sample_name, "_")
  
  replicate_number <- unlist(lapply(sample_info, function(x){x[[1]]}))
  crosslinkning_time <- unlist(lapply(sample_info, function(x){x[[2]]}))
  lysis_buffer <- unlist(lapply(sample_info, function(x){x[[3]]}))
  
  df$sample <- rep(sample_name)
  df$replicate <- rep(replicate_number)
  df$crosslinking <- rep(crosslinkning_time)
  df$lysis <- rep(lysis_buffer)
  
  
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
  
  #png(file = file.path(parent_dir, "plots", str_replace(fylename, ".csv", ".png")), height = 400, width = 800)
  print(trace_plot)
  #dev.off()
  
  #save data to list
  data_list[[fylename]] <- df
}

all_data <- do.call("rbind", data_list)

DDM_data <- all_data %>% filter(lysis == "DDM")

DDM_data$crosslinking <- factor(DDM_data$crosslinking,
                                levels = c("0mins",
                                           "5mins",
                                           "10mins",
                                           "15mins"))

DDM_data_split <- split(DDM_data, DDM_data$replicate)

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

plot_data <- function(x){
  
  tytle <- unique(x$replicate)
  
  lower <- 
    x %>% 
    filter(`Fraction Number` == 3) %>% pull(cum_vol) %>% min() 
  
  max_y <- 
    x %>% 
    filter(cum_vol > lower)%>% 
    pull(Absorbance) %>% max()
  
  max_y <- round(max_y,1) + 0.1
  
  if(tytle %in% c("REP1", "REP2", "REP3")){
    colour2use <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  } else{
    colour2use <- c("#E41A1C", "#4DAF4A")
  }
  
  if(unique(x$lysis) == "DIG"){
    tytle = "DIG"
  } else{
    tytle = tytle
    }
  
  x %>% 
    ggplot(aes(x = cum_vol, 
               y = Absorbance, 
               colour = crosslinking))+
    geom_line(size = 1)+
    scale_colour_manual(values = colour2use)+ 
    theme_bw()+
    scale_x_continuous(limits = c(0, total_vol), 
                       breaks = 1:total_vol)+
    ylim(c(NA, max_y))+
    ylab("absorbance")+
    xlab("volume")+
    ggtitle(paste0(tytle)) +
    publication_theme() +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 18),
          axis.line = element_line(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.direction = "vertical",
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
}

gplots <- lapply(DDM_data_split, plot_data)

legend <- cowplot::get_legend(gplots[[1]])



plot_data <- function(x){
  
  tytle <- unique(x$replicate)
  
  lower <- 
    x %>% 
    filter(`Fraction Number` == 3) %>% pull(cum_vol) %>% min() 
  
  max_y <- 
    x %>% 
    filter(cum_vol > lower)%>% 
    pull(Absorbance) %>% max()
  
  max_y <- round(max_y,1) + 0.1
  
  if(tytle %in% c("REP1", "REP2", "REP3")){
    colour2use <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  } else{
    colour2use <- c("#E41A1C", "#4DAF4A")
  }
  
  if(unique(x$lysis) == "DIG"){
    tytle = "DIG"
  } else{
    tytle = tytle
  }
  
  x %>% 
    filter(cum_vol >= 1) %>% 
    ggplot(aes(x = cum_vol, 
               y = Absorbance, 
               colour = crosslinking))+
    geom_line(size = 1)+
    scale_colour_manual(values = colour2use)+ 
    theme_bw()+
    scale_x_continuous(limits = c(1, total_vol), 
                       breaks = 1:total_vol)+
    ylim(c(NA, max_y))+
    ylab("absorbance")+
    xlab("volume")+
    ggtitle(paste0(tytle)) +
    publication_theme() +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 18),
          axis.line = element_line(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          legend.position = "none",
          legend.direction = "horizontal",
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
}

gplots <- lapply(DDM_data_split, plot_data)

gplot <- cowplot::plot_grid(gplots[[1]],
                   gplots[[2]],
                   gplots[[3]],
                   gplots[[4]],
                   rel_heights = c(0.4,0.4,0.4,0.4),
                   label_x = "x axis",
                   ncol = 1)

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/optimiastions_redone.png",res = 300, height = 2000, width = 1250)
print(gplot)
dev.off()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/optimiastions_redone_legend.png",res = 300, height = 500, width = 500)
print(cowplot::ggdraw(legend))
dev.off()

for(i in 1:length(gplots)){
  
  rep_num <- paste0("optimisation_REP", i, ".png")
  
  png(filename = file.path("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures",
                           rep_num),
      height = 1000,
      width = 1750,
      res = 300)
  print(gplots[i])
  dev.off()
  
}

DIG_plot <- 
all_data %>% 
  filter(lysis == "DIG") %>% 
  plot_data()

png(filename = file.path("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/DIG_plot.png"),
    height = 1000,
    width = 1750,
    res = 300)
print(DIG_plot)
dev.off()
