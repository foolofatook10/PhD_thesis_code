###Imports
library(tidyverse)

#set parent_dir as the folder with the csv file in----
parent_dir <- "\\\\data.beatson.gla.ac.uk/data/JETTLES/CNOT3_Riboseq/gradients_final"
save_location <- "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures"

#set the lower and upper limit for the y axis
ylims <- c(0, 0.7)

#make a list of files to plot data from
setwd(file.path(parent_dir, "traces"))
fylenames <- dir(pattern = "_digested.csv") # creates the list of all the csv files in the directory

#use a for loop to read in the data, reformat it, plot the trace a save the reformatted data to a list
data_list <- list() #creates a list to save the data in

fraction_lines <-
  data.frame(sample = c(rep("NTC", 5), rep("CNOT3", 5), rep("eS25", 5)),
             replicate = c(rep(1:5, 3)),
             fraction_collection = c("2-5", "4-7", "3-6", "3-6", "2-5",
                                     "4-7", "5-8", "3-6", "3-6", "3-6",
                                     "3-6", "4-7", "3-6", "3-6", "2-5"),
             axis_y = c(0.5, 0.35, 0.65, 0.5, 0.5,
                        0.5, 0.35, 0.65, 0.6, 0.6,
                        0.4, 0.35, 0.65, 0.4, 0.5))

fraction_lines$fraction_start <-  as.numeric(unlist(lapply(str_split(fraction_lines$fraction_collection, "-"), function(x){x[[1]]})))
fraction_lines$fraction_end <- as.numeric(unlist(lapply(str_split(fraction_lines$fraction_collection, "-"), function(x){x[[2]]})))


for (fylename in fylenames){
  
  print(fylename)
  #fylename = "CNOT3IP_REPLICATE3_23_10_30_digested1.csv"
  
  df <- read_csv(file = fylename, col_names = T, skip = 45) #skip the first 45 rows which contain information about the run settings
  
  split_name <- str_split(fylename, "_")
  
  REP <- as.numeric(str_remove(unlist(lapply(split_name, function(x){x[[5]]})), "Riboseq"))
  
  samp <- unlist(lapply(split_name, function(x){x[[6]]}))
  
  if(samp == "RPS25"){
    
    samp = str_replace(samp, pattern = "RPS", "eS")
    
  }
  
  fraction_start <- fraction_lines %>% filter(sample == samp & replicate == REP) %>% pull(fraction_start)
  fraction_end <- fraction_lines %>% filter(sample == samp & replicate == REP) %>% pull(fraction_end)
  ylimit <- c(0,fraction_lines %>% filter(sample == samp & replicate == REP) %>% pull(axis_y))
  
  total_vol <- sum(as.numeric(df$`Fraction Volume(ml)`), na.rm = T)
  step <- total_vol / nrow(df)
  df$cum_vol <- seq(step, total_vol, step)
  
  #extract sample name
  sample_name <- paste0(samp, "_REP", REP)
  
  df$sample <- rep(sample_name)
  
  #plot data
  df$`Fraction Number` <- as.numeric(df$`Fraction Number`)
  df %>%
    filter(!(is.na(`Fraction Number`))) -> fractions
  
  xint1 <- fractions$cum_vol[fractions$`Fraction Number` == fraction_start]
  xint2 <- fractions$cum_vol[fractions$`Fraction Number` == (fraction_start + 1)]
  xint3 <- fractions$cum_vol[fractions$`Fraction Number` == (fraction_end - 1)]
  xint4 <- fractions$cum_vol[fractions$`Fraction Number` == fraction_end]
  
  x_lab <- "Volume (ml)"
  y_lab <- "UV Absorbance (254nm)"
  
  
  
  trace_plot <- ggplot(data = df, aes(x = cum_vol, y = Absorbance))+
    geom_line(size = 1)+
    theme_bw()+
    scale_x_continuous(limits = c(0, total_vol), breaks = 1:total_vol)+
    coord_cartesian(xlim = c(3,NA)) +
    ylim(ylimit)+
    ylab(y_lab)+
    xlab(x_lab)+
    #geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = labs), 
     #         alpha = 0.2, data = rect_df, inherit.aes = FALSE) +
    annotate("rect", xmin = xint1, xmax = xint2, ymin = -Inf, ymax = Inf,
             fill = "blue", alpha = 0.2) +
    annotate("rect", xmin = xint3, xmax = xint4, ymin = -Inf, ymax = Inf,
             fill = "red", alpha = 0.2) +
  ggtitle(paste0(samp, " REP", REP))+
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 18),
          axis.line = element_line(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
          legend.text = element_text(size = 20),
          legend.title = element_blank())
  
  print(trace_plot)
  
  png(file = file.path(save_location, str_replace(fylename, ".csv", ".png")), height = 1000, width = 1250, res = 300)
  print(trace_plot)
  dev.off()
  
}


rect_df <- data.frame(xmin = c(xint1, xint3), xmax = c(xint2, xint4),
                      labs = c("Monosome", "Disome"))

rect_df$labs <- factor(rect_df$labs, levels = c("Monosome", "Disome"))

trace_plot <- ggplot(data = df, aes(x = cum_vol, y = Absorbance))+
  geom_line(size = 1)+
  theme_bw()+
  scale_x_continuous(limits = c(0, total_vol), breaks = 1:total_vol)+
  coord_cartesian(xlim = c(3,NA)) +
  ylim(ylimit)+
  ylab(y_lab)+
  xlab(x_lab)+
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = labs), 
            alpha = 0.2, data = rect_df, inherit.aes = FALSE) +
  scale_fill_manual(values = c("blue", "red")) +
  ggtitle(paste0(samp, " REP", REP))+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 18),
        axis.line = element_line(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")

legend <- cowplot::get_legend(trace_plot)

png(file = file.path(save_location, "digested_traces_legend.png"), height = 100, width = 1000, res = 300)
cowplot::ggdraw(legend)
dev.off()

