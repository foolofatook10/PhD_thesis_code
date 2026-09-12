###Imports
library(tidyverse)

#set parent_dir as the folder with the csv file in----
parent_dir <- "\\\\data.beatson.gla.ac.uk/data/JETTLES/CNOT3_Sel_RibProf/sucrose_gradients2"
save_location <- "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures"

#set the lower and upper limit for the y axis
ylims <- c(0, 0.7)

#make a list of files to plot data from
setwd(file.path(parent_dir, "traces"))
fylenames <- dir(pattern = "digested\\d.csv") # creates the list of all the csv files in the directory

#use a for loop to read in the data, reformat it, plot the trace a save the reformatted data to a list
data_list <- list() #creates a list to save the data in

# fraction_lines <- 
#   data.table(sample = c(rep("NTC", 5), rep("CNOT3", 5), rep("RPS25", 5)),
#              replicate = c(rep(1:5, 3)),
#              fraction_collection = c("2-5", "4-7", "3-6", "3-6", "2-5",
#                                      "4-7", "5-8", "3-6", "3-6", "3-6",
#                                      "3-6", "4-7", "3-6", "3-6", "2-5"),
#              axis_y = c(0.5, 0.35, 0.65, 0.5, 0.5,
#                         0.5, 0.35, 0.65, 0.6, 0.6,
#                         0.4, 0.35, 0.65, 0.4, 0.5))
fylename = "CNOT3IP_REPLICATE3_23_10_03_digested1.csv"
for (fylename in fylenames){
  
  print(fylename)
  #fylename = "CNOT3IP_REPLICATE3_23_10_30_digested1.csv"

  df <- read_csv(file = fylename, col_names = T, skip = 45) #skip the first 45 rows which contain information about the run settings
  
  split_name <- str_split(fylename, "_")
  
  REP <- as.numeric(str_remove(unlist(lapply(split_name, function(x){x[[2]]})), "REPLICATE"))
  REP <- REP -1
  samp <- unlist(lapply(split_name, function(x){x[[6]]}))
  
  samp <- paste0("Digested", str_extract(samp, "\\d"))
  
  # which_plot <- fraction_lines %>% filter(sample == samp & replicate == REP) %>% pull(fraction_collection)
  # ylimit <- c(0,fraction_lines %>% filter(sample == samp & replicate == REP) %>% pull(axis_y))
  
  total_vol <- sum(as.numeric(df$`Fraction Volume(ml)`), na.rm = T)
  
  step <- total_vol / nrow(df)
  df$cum_vol <- seq(step, total_vol, step)
  
  #extract sample name
  sample_name <- str_remove(fylename, "\\.csv") #removes the .csv
  sample_name <- str_replace_all(sample_name, "\\_", " ") #replaces underscores with spaces
  
  df$sample <- rep(sample_name)
  
  #plot data
  df$`Fraction Number` <- as.numeric(df$`Fraction Number`)
  df %>%
    filter(!(is.na(`Fraction Number`))) -> fractions
  
  
    if(fylename == "CNOT3IP_REPLICATE3_23_10_03_digested1.csv"){
      
      xint1 = fractions$cum_vol[fractions$`Fraction Number` == 3]
      xint2 = fractions$cum_vol[fractions$`Fraction Number` == 6]
      
    trace_plot <- ggplot(data = df %>% filter(cum_vol > 1), aes(x = cum_vol, y = Absorbance))+
      geom_line(size = 1)+
      theme_bw()+
      scale_x_continuous(limits = c(1, total_vol), breaks = 1:total_vol)+
      #coord_cartesian(ylim = c(0, 1)) +
      ylim(c(0, 2.5))+
      ylab("UV Absorbance (254nm)")+
      xlab("Volume (ml)")+
      ggtitle(paste0("REP", REP," ", samp))+
      annotate("rect", xmin = xint1, xmax = xint2, ymin = -Inf, ymax = Inf,
               fill = "blue", alpha = 0.2) +
      theme(axis.title = element_text(size = 20),
            axis.text = element_text(size = 18),
            axis.line = element_line(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
            legend.text = element_text(size = 20),
            legend.title = element_blank())
    } else{
      
      xint1 = fractions$cum_vol[fractions$`Fraction Number` == 1]
      xint2 = fractions$cum_vol[fractions$`Fraction Number` == 2]
      
      trace_plot <- ggplot(data = df %>% filter(cum_vol > 1), aes(x = cum_vol, y = Absorbance))+
        geom_line(size = 1)+
        theme_bw()+
        scale_x_continuous(limits = c(1, total_vol), breaks = 1:total_vol)+
        #coord_cartesian(ylim = c(0, 1)) +
        ylim(c(0, 2.5))+
        ylab("UV absorbance (254nm)")+
        xlab("Volume (ml)")+
        ggtitle(paste0("REP", REP," ", samp))+
        annotate("rect", xmin = xint1, xmax = xint2, ymin = -Inf, ymax = Inf,
                 fill = "blue", alpha = 0.2) +
        theme(axis.title = element_text(size = 20),
              axis.text = element_text(size = 18),
              axis.line = element_line(),
              panel.border = element_blank(),
              panel.grid = element_blank(),
              plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
              legend.text = element_text(size = 20),
              legend.title = element_blank())
      
    }
    print(trace_plot)
    
    

    png(file = file.path(save_location, str_replace(fylename, ".csv", ".png")), height = 1250, width = 1500, res = 300)
    print(trace_plot)
    dev.off()
    
}
    
    trace_plot <- ggplot(data = df, aes(x = cum_vol, y = Absorbance))+
      geom_line(size = 1)+
      theme_bw()+
      scale_x_continuous(limits = c(0, total_vol), breaks = 1:total_vol)+
      ylim(ylimit)+
      ylab("absorbance")+
      xlab("volume")+
      ggtitle(paste0(samp, " REP", REP))+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 1])+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 2])+
      geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 3])+
      geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 4])+
      geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 5])+
      geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 6])+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 7])+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 8])+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 9])+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 10])+
      theme(axis.title = element_text(size = 20),
            axis.text = element_text(size = 18),
            axis.line = element_line(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
            legend.text = element_text(size = 20),
            legend.title = element_blank())
    
    #png(file = file.path(save_location, str_replace(fylename, ".csv", ".png")), height = 1000, width = 1500, res = 300)
    print(trace_plot)
    #dev.off()
    
  } else if(which_plot == "4-7") {
    
    trace_plot <- ggplot(data = df, aes(x = cum_vol, y = Absorbance))+
      geom_line(size = 1)+
      theme_bw()+
      scale_x_continuous(limits = c(0, total_vol), breaks = 1:total_vol)+
      ylim(ylimit)+
      ylab("absorbance")+
      xlab("volume")+
      ggtitle(paste0(samp, " REP", REP))+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 1])+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 2])+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 3])+
      geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 4])+
      geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 5])+
      geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 6])+
      geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 7])+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 8])+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 9])+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 10])+
      theme(axis.title = element_text(size = 20),
            axis.text = element_text(size = 18),
            axis.line = element_line(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
            legend.text = element_text(size = 20),
            legend.title = element_blank())
    
    #png(file = file.path(save_location, str_replace(fylename, ".csv", ".png")), height = 1000, width = 1500, res = 300)
    print(trace_plot)
    #dev.off()
    
  } else{
    
    trace_plot <- ggplot(data = df, aes(x = cum_vol, y = Absorbance))+
      geom_line(size = 1)+
      theme_bw()+
      scale_x_continuous(limits = c(0, total_vol), breaks = 1:total_vol)+
      ylim(ylimit)+
      ylab("absorbance")+
      xlab("volume")+
      ggtitle(paste0(samp, " REP", REP))+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 1])+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 2])+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 3])+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 4])+
      geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 5])+
      geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 6])+
      geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 7])+
      geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 8])+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 9])+
      #geom_vline(xintercept = fractions$cum_vol[fractions$`Fraction Number` == 10])+
      theme(axis.title = element_text(size = 20),
            axis.text = element_text(size = 18),
            axis.line = element_line(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
            legend.text = element_text(size = 20),
            legend.title = element_blank())
    
    #png(file = file.path(save_location, str_replace(fylename, ".csv", ".png")), height = 1000, width = 1500, res = 300)
    print(trace_plot)
    #dev.off()
    
    
  }
  
  #save data to list
  data_list[[fylename]] <- df
}
all_data <- do.call("rbind", data_list)

all_data %>%
  filter(sample == "2025 04 29 CNOT3RPS25 Riboseq1 CNOT3 digested") %>%
  ggplot(aes(x = cum_vol, y = Absorbance))+
  geom_line(size = 1)+
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
        legend.title = element_blank()) -> trace_plot

png(file = file.path(parent_dir, "plots/A2_overlaid_dig_1_traces.png"), height = 400, width = 1000)
print(trace_plot)
dev.off()


