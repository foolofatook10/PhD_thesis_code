###This script was written by Joe and will calculate subs/polys ratios for each trace 
###Within the parent directory (set below) the traces need to be in a directory named "traces" and there needs to be an empty directory called "plots" to save these too.

#Imports
library(tidyverse)

#set parent_dir as the folder with the csv files in----
parent_dir <- "\\\\data.beatson.gla.ac.uk/data/JETTLES/CNOT3_Sel_RibProf/optimisations/crosslinking_gradients"

#set the what you believe the baseline to be (this will effect the magnitude of subs:polys ratios but not the direction). If a blank was carried out use this as a guide
baseline <- 0.12

#set end of polys (this can remove weird things in the traces which sometimes occurs towards the end, 12 is sensible unless there is a good reason overwise)
max_fraction <- 12

#set the upper ylim for plots
upper_ylim <- 2

#make a list of files to plot data from
setwd(file.path(parent_dir, "traces"))
fylenames <- dir(pattern = "DDM.csv") # creates the list of all the csv files in the directory

#use a for loop to read in the data, reformat it, calculate where the subs and polys start and plot the trace with these values included for visual inspection
#It then calculates the area under the curve (minus the baseline) for subs and polys and saves as a data frame to a list
data_list <- list() #creates a list to save the data to

fylename = fylenames[[2]]

gplot_lyst <- list()

for (fylename in fylenames){
  #read in data
  df <- read_csv(file = fylename, col_names = T, skip = 45) #skip the first 45 rows which contain information about the run settings
  
  df$Absorbance <- df$Absorbance + 0.5 
  
  replic <- as.numeric(str_remove(unlist(lapply(str_split(fylename, "_"),function(x){x[[1]]})), "REP"))
  conditioner <- unlist(lapply(str_split(fylename, "_"),function(x){x[[2]]}))
  lysis <- str_remove(unlist(lapply(str_split(fylename, "_"),function(x){x[[3]]})), ".csv")
  
  samp <- paste0(conditioner,replic)
  
    
  if(replic ==4){
    
    S80 = 3.2
  } else{
    
    S80 = 3
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
    filter(cum_vol < 2 & cum_vol > 1) %>% 
    top_n(n = -1, wt = Absorbance) %>% #selects the row with minimum absorbance
    pull(cum_vol) -> start_40s
  
  df %>%
    filter(cum_vol < 3 & cum_vol > 2) %>% 
    top_n(n = -1, wt = Absorbance) %>% #selects the row with minimum absorbance
    pull(cum_vol) -> start_60s
  
  #calculate end of 80S peak
  #The following pipe finds the lowest absorbance between fractions 4 and 5, which typically should be the end of the 80S peak but needs manual inspection to be sure
  #adjust the values below if required
  df %>%
    filter(cum_vol > S80 & cum_vol < 4) %>%
    top_n(n = -1, wt = Absorbance) %>% #selects the row with minimum absorbance
    pull(cum_vol) -> start_80s
  
  df %>%
    filter(cum_vol > 4 & cum_vol < 5) %>%
    top_n(n = -1, wt = Absorbance) %>% #selects the row with minimum absorbance
    pull(cum_vol) -> start_disome
  
 

  
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
    #geom_vline(xintercept = start_trisome, linetype="dashed")+
    #geom_vline(xintercept = end_trisome, linetype="dashed")+
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
  
  gplot_lyst[[paste0(samp)]] <- trace_plot
  
  #png(file = file.path(parent_dir, "subs_polys", str_replace(fylename, ".csv", ".png")), height = 400, width = 800)
  #print(trace_plot)
  #dev.off()
  



  #quantify subs/polys area under the curve
  df %>%
    filter(cum_vol >= start_40s & cum_vol <= start_60s) %>%
    mutate(corrected_absorbance = Absorbance - baseline) %>%
    summarise(subs = sum(corrected_absorbance)) %>%
    pull(subs) -> R40S
  
  df %>%
    filter(cum_vol >= start_60s & cum_vol <= start_80s) %>%
    mutate(corrected_absorbance = Absorbance - baseline) %>%
    summarise(subs = sum(corrected_absorbance)) %>%
    pull(subs) -> R60S
  
  df %>%
    filter(cum_vol >= start_80s & cum_vol <= start_disome) %>%
    mutate(corrected_absorbance = Absorbance - baseline) %>%
    summarise(subs = sum(corrected_absorbance)) %>%
    pull(subs) -> R80S
  
  df %>%
    filter(cum_vol >= start_disome & cum_vol <= max_fraction) %>%
    mutate(corrected_absorbance = Absorbance - baseline) %>%
    summarise(subs = sum(corrected_absorbance)) %>%
    pull(subs) -> polys
  
  subs = (R40S + R60S + R80S)
  
  
  #save subs/polys area under curves to a list as a data frame
  data_list[[fylename]] <- data.frame(sample = sample_name, 
                                      R40S = R40S,
                                      R60S = R60S,
                                      R80S = R80S,
                                      #disome = disome,
                                      #trisome = trisome,
                                      #polys = polys,
                                      subs = subs,
                                      polys = polys,
                                      all = subs + polys)
}


gplot_lyst[[1]]
gplot_lyst[[2]]
gplot_lyst[[3]]
gplot_lyst[[4]]

gplot_lyst[[5]]
gplot_lyst[[6]]
gplot_lyst[[7]]
gplot_lyst[[8]]

gplot_lyst[[9]]
gplot_lyst[[10]]
gplot_lyst[[11]]

gplot_lyst[[12]]
gplot_lyst[[13]]

#quantify subs vs polys
subs_vs_polys <- do.call("rbind", data_list) #row binds all data frames in the list

subs_vs_polys <- 
  subs_vs_polys %>% 
  rownames_to_column("sample_id") 

names_split <- str_split(subs_vs_polys$sample_id, "_")

samples_cleaned <- unlist(lapply(names_split, function(x){x[[2]]}))
replicates_cleaned <- as.numeric(str_remove(unlist(lapply(names_split, function(x){x[[1]]})), "REP"))

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
            axis.title = element_text(face = "bold",size = rel(1)),
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


AUC <- 
  subs_vs_polys %>% 
  #filter(sample %in% c("NTC", "CNOT3")) %>% 
  mutate(R40S_ratio = R40S/all * 100,
         R60S_ratio = R60S/all * 100,
         R80S_ratio = R80S/all * 100,
         subs_ratio = subs/all * 100,
         polys_ratio = polys/all * 100) %>% 
         #disome_ratio = disome/all * 100,
         #trisome_ratio = trisome/all * 100,
         #light_polys = lpolys/all * 100,
         #h_polys = heavy_polys/all * 100) %>% 
  select(sample, replicate, R40S_ratio:polys_ratio) %>% 
  gather(key = "population", value = "auc_percent", R40S_ratio:polys_ratio)  

AUC$population <- factor(AUC$population, 
                                levels = c("R40S_ratio", 
                                           "R60S_ratio",
                                           "R80S_ratio",
                                           "subs_ratio",
                                           "polys_ratio"),
                                labels = c("40S",
                                           "60S",
                                           "80S",
                                           "Subs",
                                           "Polys"))

AUC$sample <- factor(AUC$sample,
                            levels = c("0mins", "5mins", "10mins", "15mins"))


#CNOT3_data$replicate <- as.character(CNOT3_data$replicate)

gplot <- 
  AUC %>% 
  filter(population %in% c("Subs", "Polys")) %>% 
  ggplot(aes(x = population, y = auc_percent, colour = sample)) +
  geom_boxplot(position = position_dodge(width = 0.9), outliers = F) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9), size = 2) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0,NA)) +
  ylab("Normalised AUC") +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "none",
        #legend.direction = "vertical",
        legend.title = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_traces_quantification.png",
    res = 300, height = 1500, width = 1000)
print(gplot)
dev.off()

aov_data <- 
AUC %>% 
  filter(population %in% c("Subs", "Polys")) %>% 
  mutate(condition = paste0(population, "_", sample))

aov_object <- aov(formula = auc_percent ~ condition, data = aov_data)

summary(aov_object)

tukey_data <- as.data.frame(TukeyHSD(aov_object)$condition)

tukey_data %>% 
  rownames_to_column("comparison") %>% 
  dplyr::rename(padj = `p adj`) %>% 
  filter(`padj` < 0.05)
