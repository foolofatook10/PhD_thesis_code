###This script was written by Joe and will calculate subs/polys ratios for each trace 
###Within the parent directory (set below) the traces need to be in a directory named "traces" and there needs to be an empty directory called "plots" to save these too.

publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = 28, hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.title.x = element_text(size = 16, vjust = -0.2),
            axis.text.y = element_text(size = 16), 
            axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.text = element_text(size = 12),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.75, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

#Imports
library(tidyverse)

#set parent_dir as the folder with the csv files in----
parent_dir <- "\\\\data.beatson.gla.ac.uk/data/JETTLES/CNOT3_Riboseq/gradients_final"

#set the what you believe the baseline to be (this will effect the magnitude of subs:polys ratios but not the direction). If a blank was carried out use this as a guide
baseline <- 0.12

#set end of polys (this can remove weird things in the traces which sometimes occurs towards the end, 12 is sensible unless there is a good reason overwise)
max_fraction <- 12

#set the upper ylim for plots
upper_ylim <- 0.5

#make a list of files to plot data from
setwd(file.path(parent_dir, "Undigested_traces"))
fylenames <- dir(pattern = ".csv") # creates the list of all the csv files in the directory

#use a for loop to read in the data, reformat it, calculate where the subs and polys start and plot the trace with these values included for visual inspection
#It then calculates the area under the curve (minus the baseline) for subs and polys and saves as a data frame to a list
data_list <- list() #creates a list to save the data to

fylename = fylenames[[1]]

for (fylename in fylenames){
  
  print(fylename)
  
  #read in data
  df <- read_csv(file = fylename, col_names = T, skip = 45) #skip the first 45 rows which contain information about the run settings
  
  #replic <- as.numeric(str_remove(unlist(lapply(str_split(fylename, "_"),function(x){x[[5]]})), "Riboseq"))
  conditioner <- unlist(lapply(str_split(fylename, "_"),function(x){x[[1]]}))
  replycate <- str_remove_all(unlist(lapply(str_split(fylename, "_"),function(x){x[[2]]})), ".csv")
  
  samp <- str_remove(fylename, ".csv")
  
  if(samp %in% c("NTC_REP1", "RPS25_REP1", "RPS25_REP2", "CNOT3_REP4", "NTC_REP5")){
    
    trisome_start_sum_col = 8.3
    trisome_end_cum_vol = 10
    
  } else{
    
    trisome_start_sum_col = 8
    trisome_end_cum_vol = 9
  }
  
  if(samp == "RPS25_REP5"){
    
    fortyS_start = 3.5
  } else{
    
    fortyS_start = 4
  }
  
  #create a continuous scale of cumulative volume to plot against absorbance
  total_vol <- sum(as.numeric(df$`Fraction Volume(ml)`), na.rm = T)
  step <- total_vol / nrow(df)
  df$cum_vol <- seq(step, total_vol, step)
  
  #extract sample name
  sample_name <- str_remove(fylename, "\\.csv") #removes the .csv
  sample_name <- str_replace_all(sample_name, "\\_", " ") #replaces underscores with spaces
  df$sample <- rep(sample_name)
  
  #calculate start of 40S peak
  #The following pipe finds the lowest absorbance between fractions 0.5 and 2, which should be the start of the 40S peak but needs manual inspection to be sure
  df %>%
    filter(cum_vol < fortyS_start & cum_vol > 3) %>% 
    top_n(n = -1, wt = Absorbance) %>% #selects the row with minimum absorbance
    pull(cum_vol) -> start_40s
  
  df %>%
    filter(cum_vol < 5 & cum_vol > 4) %>% 
    top_n(n = -1, wt = Absorbance) %>% #selects the row with minimum absorbance
    pull(cum_vol) -> start_60s
  
  #calculate end of 80S peak
  #The following pipe finds the lowest absorbance between fractions 4 and 5, which typically should be the end of the 80S peak but needs manual inspection to be sure
  #adjust the values below if required
  df %>%
    filter(cum_vol > 5 & cum_vol < 6) %>%
    top_n(n = -1, wt = Absorbance) %>% #selects the row with minimum absorbance
    pull(cum_vol) -> start_80s
  
  df %>%
    filter(cum_vol > 6 & cum_vol < 7) %>%
    top_n(n = -1, wt = Absorbance) %>% #selects the row with minimum absorbance
    pull(cum_vol) -> start_disome
  
  df %>%
    filter(cum_vol > 7 & cum_vol < 8) %>%
    top_n(n = -1, wt = Absorbance) %>% #selects the row with minimum absorbance
    pull(cum_vol) -> start_trisome
  
  df %>%
    filter(cum_vol > trisome_start_sum_col & cum_vol < trisome_end_cum_vol) %>%
    top_n(n = -1, wt = Absorbance) %>% #selects the row with minimum absorbance
    pull(cum_vol) -> end_trisome
  
  #plot data to check the thresholds calculated above are appropriate
  trace_plot <- ggplot(data = df, aes(x = cum_vol, y = Absorbance))+
    geom_line(size = 1)+
    theme_bw()+
    scale_x_continuous(limits = c(0, total_vol), breaks = 1:total_vol)+
    ylab("absorbance")+
    xlab("volume")+
    ylim(c(0, upper_ylim))+
    ggtitle(sample_name)+
    geom_vline(xintercept = start_40s, linetype="dashed")+
    geom_vline(xintercept = start_60s, linetype="dashed")+
    geom_vline(xintercept = start_80s, linetype="dashed")+
    geom_vline(xintercept = start_disome, linetype="dashed")+
    geom_vline(xintercept = start_trisome, linetype="dashed")+
    geom_vline(xintercept = end_trisome, linetype="dashed")+
    geom_vline(xintercept = max_fraction, linetype="dashed")+
    geom_hline(yintercept = baseline, linetype="dashed")+
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 18),
          axis.line = element_line(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
          legend.text = element_text(size = 20),
          legend.title = element_blank())
  
  #png(file = file.path(parent_dir, "subs_polys", str_replace(fylename, ".csv", ".png")), height = 400, width = 800)
  print(trace_plot)
  #dev.off()
  
  #quantify subs/polys area under the curve
  df %>%
    filter(cum_vol > start_40s & cum_vol <= start_60s) %>%
    mutate(corrected_absorbance = Absorbance - baseline) %>%
    summarise(subs = sum(corrected_absorbance)) %>%
    pull(subs) -> R40S
  
  df %>%
    filter(cum_vol > start_60s & cum_vol <= start_80s) %>%
    mutate(corrected_absorbance = Absorbance - baseline) %>%
    summarise(subs = sum(corrected_absorbance)) %>%
    pull(subs) -> R60S
  
  df %>%
    filter(cum_vol > start_80s & cum_vol <= start_disome) %>%
    mutate(corrected_absorbance = Absorbance - baseline) %>%
    summarise(subs = sum(corrected_absorbance)) %>%
    pull(subs) -> R80S
  
  df %>%
    filter(cum_vol > start_disome & cum_vol <= start_trisome) %>%
    mutate(corrected_absorbance = Absorbance - baseline) %>%
    summarise(subs = sum(corrected_absorbance)) %>%
    pull(subs) -> disome
  
  df %>%
    filter(cum_vol > start_trisome & cum_vol <= end_trisome) %>%
    mutate(corrected_absorbance = Absorbance - baseline) %>%
    summarise(subs = sum(corrected_absorbance)) %>%
    pull(subs) -> trisome
  
  df %>%
    filter(cum_vol > end_trisome & cum_vol <= max_fraction) %>%
    mutate(corrected_absorbance = Absorbance - baseline) %>%
    summarise(subs = sum(corrected_absorbance)) %>%
    pull(subs) -> hpolys
  
  df %>%
    filter(cum_vol > start_40s & cum_vol <= start_disome) %>%
    mutate(corrected_absorbance = Absorbance - baseline) %>%
    summarise(polys = sum(corrected_absorbance)) %>%
    pull(polys) -> subs
  
  df %>%
    filter(cum_vol > start_disome & cum_vol <= max_fraction) %>%
    mutate(corrected_absorbance = Absorbance - baseline) %>%
    summarise(polys = sum(corrected_absorbance)) %>%
    pull(polys) -> polys
  
  #save subs/polys area under curves to a list as a data frame
  data_list[[fylename]] <- data.frame(sample = sample_name, 
                                      R40S = R40S,
                                      R60S = R60S,
                                      R80S = R80S,
                                      disome = disome,
                                      trisome = trisome,
                                      heavy_polys = hpolys,
                                      subs = subs,
                                      polys = polys,
                                      lpolys = (disome + trisome),
                                      all = subs + polys)
}

#quantify subs vs polys
subs_vs_polys <- do.call("rbind", data_list) #row binds all data frames in the list

subs_vs_polys <- 
  subs_vs_polys %>% 
  rownames_to_column("sample_id") 

names_split <- str_split(subs_vs_polys$sample, " ")

samples_cleaned <- unlist(lapply(names_split, function(x){x[[1]]}))
replicates_cleaned <- as.numeric(str_remove_all(unlist(lapply(names_split, function(x){x[[2]]})), pattern = "REP"))

subs_vs_polys$sample <- samples_cleaned
subs_vs_polys$replicate <- replicates_cleaned

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
            axis.title = element_text(size = 20),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(size = 18, color = "black"),
            axis.text.y = element_text(size = 18, color = "black"),
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


CNOT3_data <- 
  subs_vs_polys %>% 
  filter(sample %in% c("NTC", "CNOT3")) %>% 
  mutate(R40S_ratio = R40S/all * 100,
         R60S_ratio = R60S/all * 100,
         R80S_ratio = R80S/all * 100,
         disome_ratio = disome/all * 100,
         trisome_ratio = trisome/all * 100,
         light_polys = lpolys/all * 100,
         h_polys = heavy_polys/all * 100) %>% 
  select(sample, replicate, R40S_ratio:h_polys) %>% 
  gather(key = "population", value = "auc_percent", R40S_ratio:h_polys)  

CNOT3_data$population <- factor(CNOT3_data$population, 
                                levels = c("R40S_ratio", 
                                           "R60S_ratio",
                                           "R80S_ratio",
                                           "disome_ratio",
                                           "trisome_ratio",
                                           "light_polys",
                                           "h_polys"),
                                labels = c("40S",
                                           "60S",
                                           "80S",
                                           "disomes\n(D)",
                                           "trisomes\n(T)",
                                           "L. polys",
                                           "H. polys"))

CNOT3_data$sample <- factor(CNOT3_data$sample,
                            levels = c("NTC", "CNOT3"),
                            labels = c("NTC", "siCNOT3")) 


CNOT3_data$replicate <- as.character(CNOT3_data$replicate)

CNOT3_gplot <- 
  CNOT3_data %>%  
  filter(!population %in% c("disomes\n(D)", "trisomes\n(T)")) %>% 
  ggplot(aes(x = population, y = auc_percent, colour = sample, fill = sample, shape = replicate, group = sample)) +
  geom_boxplot(aes(x = population, y = auc_percent, colour = sample), 
               position = position_dodge(width = 0.9), 
               outliers = F,
               inherit.aes = F, show.legend = F) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.75), size = 2) +
  scale_color_manual(values = c("#FFC20A", "#88CCEE"), name = "Condition") +
  scale_fill_manual(values = c("#FFC20A", "#88CCEE"), name = "Condition") +
  scale_shape_manual(values = c(21,22,23,24,25), name = "Replicate") +
  ylab("Normalised AUC") +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "none",
        legend.direction = "vertical",
        legend.margin = unit(0.5, "cm"),
        legend.box = "vertical",
        legend.box.spacing = unit(1, "cm"))

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_traces_quantification.png",
    res = 300, height = 1250, width = 1500)
print(CNOT3_gplot)
dev.off()

library(broom)

CNOT3_data %>%
  filter(!population %in% c("disomes\n(D)", "trisomes\n(T)")) %>%
  select(population, replicate, sample, auc_percent) %>%
  pivot_wider(
    names_from = sample,
    values_from = auc_percent
  ) %>%
  group_by(population) %>%
  summarise(
    p_siCNOT3 = t.test(siCNOT3, NTC, paired = TRUE)$p.value,
    #p_sieS25 = t.test(sieS25, NTC, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_siCNOT3 = p.adjust(p_siCNOT3, method = "BH"),
    #p_adj_sieS25 = p.adjust(p_sieS25, method = "BH")
  )



RPS25_data <- 
  subs_vs_polys %>% 
  filter(sample %in% c("NTC", "RPS25")) %>% 
  mutate(R40S_ratio = R40S/all * 100,
         R60S_ratio = R60S/all * 100,
         R80S_ratio = R80S/all * 100,
         disome_ratio = disome/all * 100,
         trisome_ratio = trisome/all * 100,
         light_polys = lpolys/all * 100,
         h_polys = heavy_polys/all * 100) %>% 
  select(sample, replicate, R40S_ratio:h_polys) %>% 
  gather(key = "population", value = "auc_percent", R40S_ratio:h_polys)  

RPS25_data$population <- factor(RPS25_data$population, 
                                levels = c("R40S_ratio", 
                                           "R60S_ratio",
                                           "R80S_ratio",
                                           "disome_ratio",
                                           "trisome_ratio",
                                           "light_polys",
                                           "h_polys"),
                                labels = c("40S",
                                           "60S",
                                           "80S",
                                           "disomes\n(D)",
                                           "trisomes\n(T)",
                                           "L. polys",
                                           "H. polys"))

RPS25_data$sample <- factor(RPS25_data$sample,
                            levels = c("NTC", "RPS25"),
                            labels = c("NTC", "sieS25")) 


RPS25_data$replicate <- as.character(RPS25_data$replicate)

RPS25_gplot <- 
  RPS25_data %>%  
  filter(!population %in% c("disomes\n(D)", "trisomes\n(T)")) %>% 
  ggplot(aes(x = population, y = auc_percent, colour = sample, fill = sample, shape = replicate, group = sample)) +
  geom_boxplot(aes(x = population, y = auc_percent, colour = sample), 
               position = position_dodge(width = 0.9), 
               outliers = F,
               inherit.aes = F, show.legend = F) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.75), size = 2) +
  scale_color_manual(values = c("#FFC20A", "#984EA3"), name = "Condition") +
  scale_fill_manual(values = c("#FFC20A", "#984EA3"), name = "Condition") +
  scale_shape_manual(values = c(21,22,23,24,25), name = "Replicate") +
  ylab("Normalised AUC") +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "none",
        legend.direction = "vertical",
        legend.margin = unit(0.5, "cm"),
        legend.box = "vertical",
        legend.box.spacing = unit(1, "cm")) 



png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/RPS25_traces_quantification.png",
    res = 300, height = 1250, width = 1500)
print(RPS25_gplot)
dev.off()


RPS25_data %>%
  filter(!population %in% c("disomes\n(D)", "trisomes\n(T)")) %>%
  select(population, replicate, sample, auc_percent) %>%
  pivot_wider(
    names_from = sample,
    values_from = auc_percent
  ) %>%
  group_by(population) %>%
  summarise(
    #p_siCNOT3 = t.test(siCNOT3, NTC, paired = TRUE)$p.value,
    p_sieS25 = t.test(sieS25, NTC, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    #p_adj_siCNOT3 = p.adjust(p_siCNOT3, method = "BH"),
    p_adj_sieS25 = p.adjust(p_sieS25, method = "BH")
  )

