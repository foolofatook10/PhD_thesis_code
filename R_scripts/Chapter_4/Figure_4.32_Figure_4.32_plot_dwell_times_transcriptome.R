library(tidyverse)
library(data.table)

dts <- read_csv("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/dwell_times_CNOT3SelRP_entire_transcriptome.csv") %>% 
  select(-`...1`)

polarised <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/polarised_RPFs") %>% 
  pull("x")

# step 1 filter out poolarised transcripts
dts <- dts %>% filter(!transcript %in% polarised)

dts_summarised <- 
  dts %>% 
  group_by(sample_id, codon) %>% 
  summarise(`-4` = sum(minus4),
            `-3` = sum(minus3),
            `E` = sum(E_site),
            `P` = sum(P_site),
            `A` = sum(A_site),
            `+1` = sum(plus1),
            `+2` = sum(plus2),
            total = sum(total))

dts_summarised_split <- split(dts_summarised, dts_summarised$sample_id)

downstream_processing <- function(df){

  df <-
    df %>%
    filter(codon != "tga" & codon != "tag" & codon != "taa") %>%
    mutate(minus_4 = (log2(`-4`)) - log2(total / 7),
           minus_3 = (log2(`-3`)) - log2(total / 7),
           E_site = (log2(E)) - log2(total / 7),
           P_site = (log2(P)) - log2(total / 7),
           A_site = (log2(A)) - log2(total / 7),
           plus_1 = (log2(`+1`)) - log2(total / 7),
           plus_2 = (log2(`+2`)) - log2(total / 7)) %>%
    select(codon, minus_4, minus_3, E_site, P_site, A_site, plus_1, plus_2, sample_id) %>%
    dplyr::rename(sample = sample_id)
  
  return(df)
}



data <- 
  do.call("rbind", lapply(dts_summarised_split, downstream_processing)) %>% 
  ungroup()

#gather the data into tidy format----
data %>%
  pull(codon) %>%
  unique() -> codons

RPF_sample_names <- unique(data$sample)

gathered_list <- list()
for (sample in RPF_sample_names) {
  for (codon in codons) {
    data[data$codon == codon & data$sample == sample,] %>%
      select(-sample) %>%
      gather(key = codon, value = freq) %>%
      rename(position = codon) %>%
      mutate(codon = rep(codon),
             sample = factor(rep(sample))) -> gathered_list[[paste(codon, sample, sep = "_")]]
  }
}
gathered_data <- do.call("rbind", gathered_list)
# 
gathered_data %>%
  mutate(position = as.numeric(case_when(position == "minus_4" ~ -4,
                                         position == "minus_3" ~ -3,
                                         position == "E_site" ~ -2,
                                         position == "P_site" ~ -1,
                                         position == "A_site" ~ 0,
                                         position == "plus_1" ~ 1,
                                         position == "plus_2" ~ 2)),
         wobble = factor(str_sub(codon, 3,3))) -> plot_data

for (sample in RPF_sample_names) {
  plot_title <- str_replace_all(sample, "_", " ")
  plot_title <- str_remove(plot_title, " RPFs")
  
  plot_data[plot_data$sample == sample,] %>%
    ggplot(aes(x = position, y = freq, colour = wobble))+
    geom_point(size = 2)+
    stat_summary(fun=mean, geom="line", size = 1)+
    ylab("normalised codon frequency")+
    scale_x_continuous(limits = c(-4,2), breaks = -4:2)+
    xlab("Codon position (0 = A-site)")+
    ggtitle(plot_title) -> codon_plot
  
  print(codon_plot)
  
}
# 
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "-3") %>% arrange(-freq)})
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "-2") %>% arrange(-freq)})
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "-1") %>% arrange(-freq)})
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "0") %>% arrange(-freq)})
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "1") %>% arrange(-freq)})
# 
# 
# lapply(split(plot_data, plot_data$sample), function(x){x %>% filter(position == "-1") %>% arrange(-freq)})
# 
#        
mean_data <- 
  plot_data %>% 
  mutate(condition = unlist(lapply(str_split(sample, "_"),function(x){x[[3]]}))) %>% 
  group_by(condition, position, codon, wobble) %>% 
  summarise(mean_freq = mean(freq)) 

cca <- mean_data %>% filter(codon == "cca" & position == -2)
cct <- mean_data %>% filter(codon == "cct"& position == -2)
gaa <- mean_data %>% filter(codon == "gaa"& position == 0)
tac <- mean_data %>% filter(codon == "tac"& position == 0)
cgg <- mean_data %>% filter(codon == "cgg"& position == -1)
cga <- mean_data %>% filter(codon == "cga"& position == -1)
agg <- mean_data %>% filter(codon == "agg"& position == -1)
gat <- mean_data %>% filter(codon == "gat"& position == -3)
cgt <- mean_data %>% filter(codon == "cgt"& position == 1)
# 
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
            plot.margin=unit(c(2,6,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

mean_data_final <- 
  mean_data %>% 
  mutate(to_col = case_when(codon == "agg" & position == -1 ~ "Arg-AGG",
                            codon == "cga" & position == -1 ~ "Arg-CGA",
                            codon == "cgg" & position == -1 ~ "Arg-CGG",
                            codon == "gaa" & position == 0 ~ "Glu-GAA",
                            codon == "tac" & position == 0 ~ "Tyr-TAC",
                            codon == "cca" & position == -2 ~ "Pro-CCA",
                            codon == "cct" & position == -2 ~ "Pro-CCT",
                            codon == "ccg" & position == -2 ~ "Pro-CCG",
                            codon == "ccc" & position == -2 ~ "Pro-CCC",.default = "other")) 

mean_data_final$to_col <- 
  factor(mean_data_final$to_col, levels = c("Arg-AGG", "Arg-CGA", "Arg-CGG",
                                            "Glu-GAA","Tyr-TAC",
                                            "Pro-CCA", "Pro-CCT", "Pro-CCG", "Pro-CCC", "other"))


mean_data_final$condition <- factor(mean_data_final$condition,
                                    levels = c("Tot", "CNOT3"))
set.seed(123456788)

final_gplot <- 
  mean_data_final %>% 
  ggplot(aes(x = position, y = mean_freq, colour = to_col)) +
  geom_point(size = 2, position = position_jitter(width = 0.1)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#984EA3", "#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F", "#40E0D0", "#E5C494", "#d3d3d3")) +
  facet_wrap(~condition, ncol = 2) +
  scale_x_continuous(limits = c(-4,2), 
                     breaks = -4:2, 
                     labels = c("-2", "-1", "E", "P", "A", "+1", "+2"))+
  #xlab("Codon position (0 = A-site)") +
  ylab("Log2 enrichment score") +
  publication_theme() +
  ylim(-1,1) +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank())


png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/dwell_times_transcriptome.png",
    res = 300, height = 1000, width = 2250)
print(final_gplot)
dev.off()

diffs <- 
  mean_data_final %>% 
  spread(key = condition, value = mean_freq) %>% 
  mutate(differential = CNOT3-Tot) 

diffs$position %>% unique()

diffs_coloured <- 
  diffs %>% 
  mutate(to_col2 = case_when(codon == "agg" & position == -1 ~ "Arg-AGG",
                             codon == "cga" & position == -1 ~ "Arg-CGA",
                             codon == "cgg" & position == -1 ~ "Arg-CGG",
                             codon == "tcg" & position == -1 ~ "Ser-TCG",
                             codon == "tca" & position == -1 ~ "Ser-TCA",
                             codon == "tgt" & position == -1 ~ "Cys-TGT",
                             codon == "ata" & position == 0 ~ "Ile-ATA",
                             codon =="gaa" & position == 0 ~ "Glu-GAA",
                             .default = "other"))

diffs_coloured$to_col2 <- factor(diffs_coloured$to_col2,
                                 levels = c("Arg-AGG", "Arg-CGA", "Arg-CGG",
                                            "Ser-TCG", "Ser-TCA", "Cys-TGT",
                                            "Ile-ATA", "Glu-GAA", "other")) 

diff_gplot <-  
  ggplot(diffs_coloured, aes(x = position, y = differential, colour = to_col2)) +
  geom_point(size = 2, position = position_jitter(width = 0.1)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#984EA3", "#4DAF4A", "#BF5B17", "#FF7F00", "#F0027F", "#40E0D0", "#d3d3d3")) +
  #facet_wrap(~condition, ncol = 2) +
  scale_x_continuous(limits = c(-4,2), 
                     breaks = -4:2, 
                     labels = c("-2", "-1", "E", "P", "A", "+1", "+2"))+
  #xlab("Codon position (0 = A-site)") +
  ylab("Differential enrichment") +
  ggtitle("Transcriptome") +
  publication_theme() +
  ylim(-0.5,0.5) +
  #coord_fixed() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank())


png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/differential_gplot.png",
    res = 300, height = 1250, width = 1500)
print(diff_gplot)
dev.off()

mean_data_final2 <- 
  mean_data_final %>% 
  mutate(AA_colour = case_when(codon %in% c("ccc", "ccg", "cca", "cct") ~ "Pro",
                               codon %in% c("cga", "cgg", "cgt", "cgc", "aga", "agg") ~"Arg",
                               codon %in% c("gaa", "gag") ~ "Glu",
                               #codon %in% c("ggg", "ggc", "gga", "ggt") ~ "Gly",
                               .default = "other"
  )) %>% 
  spread(key = condition, value = mean_freq)

mean_data_final2 

data_split <- split(mean_data_final2, mean_data_final2$position)

library(ggsci)

# pal_npg("nrc")(9)
# "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF"
# [8] "#DC0000FF" "#7E6148FF"


create_gplot <- 
  function(x){
    
    tytle <- case_when(unique(x$position) == "-1" ~ "P",
                       unique(x$position) == "-2" ~ "E",
                       unique(x$position) == "0" ~ "A",
                       unique(x$position) == "-3" ~ "-3",
                       unique(x$position) == "-4" ~ "-4",
                       unique(x$position) == "1" ~ "+1",
                       unique(x$position) == "2" ~ "+2",
                       .default = "Not P")
    
    axis_start <- case_when(tytle == "E"~ -0.6,
                            tytle == "P" ~-0.75,
                            tytle == "A" ~ -0.75)
    
    axis_end <- case_when(tytle == "E"~ 0.6,
                          tytle == "P" ~0.75,
                          tytle == "A" ~ 0.75)
    
    dat_lab <- x %>% filter(AA_colour != "other")
    
    x$AA_colour <- factor(x$AA_colour,
                          levels = c("Pro", "Arg", "Glu", "other"))
    
    my_colours <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF","#A3A3A3")
    
    x %>% 
      ggplot(aes(x = Tot, y = CNOT3, colour = AA_colour)) +
      geom_point() +
      ggrepel::geom_text_repel(data = dat_lab, 
                               aes(x = Tot, y = CNOT3, label = codon), 
                             show.legend = F) +
    ggtitle(paste0(tytle)) +
    scale_color_manual(values = my_colours) +
    geom_abline(lty = "dashed") +
    xlim(c(axis_start, axis_end)) +
    ylim(c(axis_start, axis_end)) +
    publication_theme() +
    theme(legend.title = element_blank())
  
}

gplots_bothdts <- lapply(data_split, create_gplot)

E <- gplots_bothdts[[3]]
P <- gplots_bothdts[[4]]
A <- gplots_bothdts[[5]]

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/RDTs_transcriptome_Esite.png", 
    res = 300,
    height = 1250,
    width = 1250)
print(E)
dev.off()

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/RDTs_transcriptome_Psite.png", 
    res = 300,
    height = 1250,
    width = 1250)
print(P)
dev.off()

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/RDTs_transcriptome_Asite.png", 
    res = 300,
    height = 1250,
    width = 1250)
print(A)
dev.off()
