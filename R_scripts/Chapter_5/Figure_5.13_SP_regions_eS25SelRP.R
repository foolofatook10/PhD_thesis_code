library(tidyverse)
library(data.table)

setwd("\\\\data.beatson.gla.ac.uk/data/R11/external_sequencing_data/rps25_sel-RiboSeq/CPMs")
save_path = "\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/plots/signal_sequences"
# read in CPMs
CDS <- as_tibble(fread("CDS_CPMs.csv"))

# read in features
features <- as_tibble(fread("\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/features_mouse.csv", 
                            header = T, 
                            drop = "V1"))

#CDS$MD30 <- factor(CDS$MD30, levels = c("M", "D30"))

#summary(CDS)

#CDS2 <- CDS %>% filter(MD30 != "D30")

#CDS2$MD30 <- droplevels(CDS2$MD30)

#summary(CDS2)

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
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

signal_sequences <- 
  features %>% 
  select(ENSMUST, signal_sequence)

# filter transcripts that are over 100 codons long
over100codons <- 
  features %>% 
  mutate(cds_length = cds_length/3) %>% 
  filter(cds_length >100) %>% 
  pull(ENSMUST)

CDS_sigseqs <- 
  CDS %>% filter(transcript %in% over100codons) %>% 
  inner_join(signal_sequences, 
             by = c("transcript" = "ENSMUST"))

CDS_sigseqs %>% 
  select(transcript) %>% 
  unique() %>% 
  inner_join(features %>% select(ENSMUST, signal_sequence), by = c("transcript" = "ENSMUST")) %>% 
  pull(signal_sequence) %>% table()

CPMS <-
  CDS_sigseqs %>% 
  select(replicate, IP, transcript, signal_sequence, codon, mean_CPM_codon)


# considering just CPMs

# consider only first 30 codons
CPMS_processed <- 
  CPMS %>%
  ungroup() %>% 
  group_by(replicate, IP, signal_sequence, transcript) %>%
  slice_head(n = -1) %>% 
  slice_head(n = 500) %>% 
  ungroup()

CPMs_replicates <- 
  CPMS_processed %>% 
  group_by(replicate, IP, signal_sequence, codon) %>% 
  summarise(mean_CPM = mean(mean_CPM_codon)) %>% 
  ggplot(aes(x = codon, y = mean_CPM, colour = IP))+
  geom_line() +
  facet_grid(replicate~signal_sequence)

CPMs_average <- 
  CPMS_processed %>% 
  group_by(IP, signal_sequence, codon) %>% 
  summarise(mean_CPM = mean(mean_CPM_codon)) %>% 
  ggplot(aes(x = codon, y = mean_CPM, colour = IP))+
  geom_line() +
  facet_wrap(~signal_sequence)

deltas_replicates <- 
  CPMS_processed %>% 
  spread(key = IP, value= mean_CPM_codon) %>% 
  mutate(delta = rps25 - Total) %>% 
  group_by(replicate, signal_sequence, codon) %>% 
  summarise(mean_delta = mean(delta)) %>% 
  ggplot(aes(x = codon, y = mean_delta, colour = signal_sequence))+
  geom_line() +
  facet_wrap(~replicate)

final <- 
CPMS_processed %>% 
  spread(key = IP, value= mean_CPM_codon) %>% 
  mutate(delta = rps25 - Total) %>% 
  group_by(signal_sequence, codon) %>% 
  summarise(
    n = n(),
    mean_delta = mean(delta),
    sd_delta = sd(delta, na.rm = TRUE),
    SE_delta = sd_delta / sqrt(n))


final$se_upper <- final$mean_delta + final$SE_delta
final$se_lower <- final$mean_delta - final$SE_delta

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
            axis.text.x = element_text(size = 20, color = "black"),
            axis.text.y = element_text(size = 20, color = "black"),
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
            plot.margin=unit(c(10,2,10,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

deltas_average <- 
  final  %>% 
  ggplot(aes(x = codon, y = mean_delta, colour = signal_sequence, fill = signal_sequence))+
  geom_rect(aes(xmin=1, xmax=30, ymin=-Inf, ymax=Inf), fill = alpha("#DECBE4", 0.05), colour = alpha("#DECBE4", 0.05)) +
  geom_rect(aes(xmin=30, xmax=60, ymin=-Inf, ymax=Inf), fill = alpha("#98FBCB", 0.05), colour = alpha("#98FBCB", 0.05)) +
  geom_rect(aes(xmin=60, xmax=150, ymin=-Inf, ymax=Inf), fill = alpha("#B3CDE3", 0.05), colour = alpha("#B3CDE3", 0.05)) +
  geom_line(size = 1) +
  geom_ribbon(aes(x = codon, y = mean_delta, ymin = se_lower, ymax = se_upper),
              alpha = 0.5) +
  #geom_smooth() +
  #geom_point() +
  #geom_function() +
  ylab("Delta CPM (eS25 - Tot)") +
  xlab("Codon") +
  #ggtitle("CNOT3 SelRP") +
  coord_cartesian(xlim = c(0,500))  +
  scale_color_manual(values = c("#FF7F0EFF", "#1F77B4FF"), labels = c("No signal sequence", "signal sequence")) +
  scale_fill_manual(values = c("#FF7F0EFF", "#1F77B4FF"), labels = c("No signal sequence", "signal sequence")) +
  scale_x_continuous(breaks = c(1,100,200,300,400,500)) +
  publication_theme() +
  theme(legend.position = "none",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())



#RColorBrewer::brewer.pal(8,"Pastel1")
#"#CCEBC5" "#DECBE4" "#FED9A6" "#FFFFCC" "#E5D8BD" "#FDDAEC"
#alpha()

deltas_average_0_30 <- 
  final %>% 
  ggplot(aes(x = codon, y = mean_delta, colour = signal_sequence, fill = signal_sequence)) +
  geom_rect(aes(xmin=1, xmax=30, ymin=-Inf, ymax=Inf), fill = alpha("#DECBE4", 0.05), colour = alpha("#DECBE4", 0.05)) +
  geom_rect(aes(xmin=30, xmax=60, ymin=-Inf, ymax=Inf), fill = alpha("#98FBCB", 0.05), colour = alpha("#98FBCB", 0.05)) +
  geom_rect(aes(xmin=60, xmax=150, ymin=-Inf, ymax=Inf), fill = alpha("#B3CDE3", 0.05), colour = alpha("#B3CDE3", 0.05)) +
  geom_line(size = 1) +
  geom_ribbon(aes(x = codon, y = mean_delta, ymin = se_lower, ymax = se_upper),
              alpha = 0.5) +
  #geom_smooth() +
  #geom_point() +
  #geom_function() +
  ylab("Delta CPM (CNOT3 - Tot)") +
  xlab("Codon") +
  #ggtitle("CNOT3 SelRP") +
  coord_cartesian(xlim = c(0,150))  +
  scale_color_manual(values = c("#FF7F0EFF", "#1F77B4FF"), labels = c("No signal sequence", "signal sequence")) +
  scale_fill_manual(values = c("#FF7F0EFF", "#1F77B4FF"), labels = c("No signal sequence", "signal sequence")) +
  scale_x_continuous(breaks = c(1,30,60,100,150)) +
  publication_theme() +
  theme(legend.position = "none",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/eS25SelRP_SPs.png",
    res = 300, height = 2000, width = 1250)
print(cowplot::plot_grid(deltas_average, deltas_average_0_30, ncol = 1, align = "y"))
dev.off()


mean_deltas <- 
  CPMS_processed %>% 
  spread(key = IP, value= mean_CPM_codon) %>% 
  mutate(delta = rps25 - Total) %>% 
  group_by(transcript, signal_sequence, codon) %>% 
  summarise(mean_delta = mean(delta))


signal_peptides <-  features %>% select(ENSMUST, signal_sequence)



mean_deltas_1_30 <- mean_deltas %>% filter(codon <= 30) 
mean_deltas_30_60 <- mean_deltas %>% filter(codon >30 & codon <= 60)
mean_deltas_50_150 <- mean_deltas %>% filter(codon >60 & codon <= 150) 



mean_deltas_1_30$mean_delta <- (mean_deltas_1_30$mean_delta + abs(mean_deltas_1_30$mean_delta %>% min()) + 0.1)
mean_deltas_30_60$mean_delta <- (mean_deltas_30_60$mean_delta + abs(mean_deltas_30_60$mean_delta %>% min()) + 0.1)
mean_deltas_50_150$mean_delta <- (mean_deltas_50_150$mean_delta + abs(mean_deltas_50_150$mean_delta %>% min()) + 0.1)


mean_deltas_split <- split(mean_deltas_50_150, mean_deltas_50_150$transcript)

#mean_deltas_split[[1]] %>% pull(mean_delta)

library(smplot2)

aucs <- 
  do.call("rbind", lapply(mean_deltas_split, function(x){
    
    transcript <- unique(x$transcript)
    codons <- x$codon
    deltas <- x$mean_delta
    #SP <-unique(x$signal_sequence)
    
    auc <- sm_auc(codons, deltas)
    return(data.frame(transcript = transcript, auc = auc))
  }))

aucs <- 
  aucs %>% inner_join(signal_peptides, by = c("transcript" = "ENSMUST"))

data <- 
  aucs %>% 
  mutate(signal_sequence = case_when(signal_sequence == F ~ "Transcripts\nwithout SP",
                                     signal_sequence == T ~ "Transcripts\nwith SP"))

data$signal_sequence <- factor(data$signal_sequence,
                               levels = c("Transcripts\nwithout SP",
                                          "Transcripts\nwith SP"))

gplot <- 
  ggplot(data, aes(x = signal_sequence, y = auc, fill = signal_sequence)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  #scale_y_continuous(limits = c(0,NA))# +
  #
  #coord_cartesian(ylim = c(125,127.5)) +
  ylab("Area under curve") +
  stat_summary(fun = mean, geom = "point", size = 4, color = "red") +
  scale_fill_manual(values = c("#FF7F0EFF", "#1F77B4FF")) +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

wilcox.test(auc~signal_sequence,data = data)
