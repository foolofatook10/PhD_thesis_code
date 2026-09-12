library(tidyverse)
library(data.table)

features <- fread("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/features.csv")

signal_peptides <-  features %>% select(ENST, signal_sequence)

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


normalised_CPMs <- fread("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/riboseq_datasets/siDDX6_CDS_CPMs.csv")

#table(is.na(normalised_CPMs))

normalised_CPMs <- 
normalised_CPMs %>% 
  dplyr::rename(KD = condition) %>% 
  mutate(KD = case_when(KD == "DDX6si" ~ "siDDX6", .default = KD)) %>% 
  dplyr::rename(normalised_CPM = CPM)

normalised_CPMs <- 
  normalised_CPMs %>% 
  mutate(id = paste0(replicate, "_", KD))

normalised_CPMs_DDX6 <- normalised_CPMs %>% filter(KD %in% c("siDDX6", "NTC"))

transcripts_per_condition <- lapply(split(normalised_CPMs_DDX6, normalised_CPMs_DDX6$id), function(x){x %>% pull(transcript)})

final_transcripts <- Reduce(intersect, transcripts_per_condition)

normalised_CPMs_DDX6 <- 
  normalised_CPMs_DDX6 %>% 
  filter(transcript %in% final_transcripts)

filt_transcripts <- intersect(unique(normalised_CPMs_DDX6$transcript), features %>% filter(cds_length > 300) %>% pull(ENST))

normalised_CPMs_DDX6 <- 
  normalised_CPMs_DDX6 %>% filter(transcript %in% filt_transcripts)

#number of SPs and non SPs in this analysis
inner_join(features %>% select(ENST, signal_sequence),
           unique(normalised_CPMs %>% select(transcript)),
           by = c("ENST"  = "transcript")) %>%
  pull(signal_sequence) %>% table()


normalised_CPMs_DDX6$KD <- factor(normalised_CPMs_DDX6$KD)
#summary(normalised_CPMs)

deltas <- 
  normalised_CPMs_DDX6 %>%
  select(replicate, KD, transcript, codon, normalised_CPM) %>% 
  spread(key = "KD", value = normalised_CPM) %>%
  mutate(siDDX6_delta = siDDX6 - NTC) %>% 
  select(replicate, transcript, codon, siDDX6_delta)

table(is.na(deltas$siDDX6_delta))

final <- 
  deltas %>% 
  group_by(replicate, transcript) %>% 
  slice_head(n = -1) %>% 
  slice_head(n = 500) %>% 
  inner_join(signal_peptides, by = c("transcript" = "ENST")) %>% 
  group_by(signal_sequence, codon) %>% 
  summarise(
    n = n(),
    mean_delta = mean(siDDX6_delta),
    sd_delta = sd(siDDX6_delta, na.rm = TRUE),
    SE_delta = sd_delta / sqrt(n))

final$se_upper <- final$mean_delta + final$SE_delta
final$se_lower <- final$mean_delta - final$SE_delta

mean_deltas <- 
  deltas %>% 
  ungroup() %>% 
  group_by(transcript, codon) %>% 
  summarise(mean_delta = mean(siDDX6_delta)) %>% 
  ungroup() %>% 
  group_by(transcript) %>% 
  slice_head(n = -1)

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

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/siDDX6_SPs.png",
    res = 300, height = 2000, width = 1250)
print(cowplot::plot_grid(deltas_average, deltas_average_0_30, ncol = 1, align = "y"))
dev.off()


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
  aucs %>% inner_join(signal_peptides, by = c("transcript" = "ENST"))

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
