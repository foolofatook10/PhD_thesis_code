library(data.table)
library(tidyverse)
library(ggsci)


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

setwd("\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/reporter_assays/GAUSSIA+REDFF_TRANSFECTIONS/Tranfections_from_mini_preps")

data <- 
  read_csv("CUG_siCNOT3RPS25.csv") %>% 
  select(signal_peptide:condition) %>% na.omit()

raw <- 
  data %>%
  group_by(condition, signal_peptide, codon, sample, channel, biorep) %>%
  summarise(mean_value = mean(value))

raw_lyst <- split(raw, f = raw$signal_peptide)

for(i in 1:length(raw_lyst)){
  
  dat <- raw_lyst[[i]]
  
  #dat <- dat %>% filter(codon %in% c("CUG", "CUC", "UUA"))
  
  dat$codon <- factor(dat$condition, levels = c("NTC", "siCNOT3", "siRPS25"))
  
  dat$biorep <- as.character(dat$biorep)
  
  media <- 
    dat %>%
    filter(sample == "media") %>%
    ggplot(aes(x = condition, y = mean_value, colour = biorep)) +
    geom_point() +
    ggtitle(paste(dat$signal_peptide[1]," Media")) +
    facet_wrap(~channel) + 
    expand_limits(y = 0) +
    publication_theme()
  
  print(media)
  
  lysate <- 
    dat %>%
    filter(sample == "lysate") %>%
    ggplot(aes(x = condition, y = mean_value, colour = biorep)) +
    geom_point() +
    ggtitle(paste(dat$signal_peptide[1]," Lysate")) +
    facet_wrap(~channel) + 
    expand_limits(y = 0) +
    publication_theme()
  
  
  print(lysate)
  
}

test <- 
data %>%
  mutate(ID = paste0(channel, "_", sample)) %>% 
  select(-sample, -channel) %>% 
  spread(key = ID, value = value) %>%
  mutate(media_normalised = gaussia_media/red_firefly_lysate,
         lysate_normalised = gaussia_lysate/red_firefly_lysate) %>% 
  na.omit()

normalised_data <- 
  data %>%
  mutate(ID = paste0(channel, "_", sample)) %>% 
  select(-sample, -channel) %>% 
  spread(key = ID, value = value) %>%
  mutate(media_normalised = gaussia_media/red_firefly_lysate,
         lysate_normalised = gaussia_lysate/red_firefly_lysate) %>% 
  na.omit() %>% 
  select(condition, signal_peptide, codon, techrep, biorep, media_normalised, lysate_normalised) %>% 
  gather(key = sample, value = normalised_value, media_normalised:lysate_normalised) %>% 
  mutate(sample = str_remove(sample, "_normalised")) %>% 
  filter(!codon %in% c("CTT", "TTG"))


normalised_data$condition <- factor(normalised_data$condition, levels = c("NTC", "siCNOT3", "siRPS25"))

#summary(normalised_data)

normalised_data$biorep <- factor(normalised_data$biorep, 
                                 levels = c("1", "2", "3"))

normalised_data %>% 
  filter(signal_peptide == "SHH") %>% 
  #filter(sample == "media") %>% 
  ggplot(aes(x = condition, y = normalised_value)) +
  geom_point(alpha = 0.3) + 
  #position = position_jitter()) +
  ylab("Normalised luciferase") +
  facet_grid(biorep~sample) +
  #ylim(0,25) +
  publication_theme() +
  theme(axis.title.x = element_blank())

test <- 
  normalised_data %>% 
  group_by(condition, signal_peptide, codon, sample, biorep) %>% 
  summarise(mean_ratio = mean(normalised_value))

test_split <- split(test, f = test$signal_peptide)

percentage_normalisation <- function(x){
  
  y <- x  %>% spread(key = "condition", value = "mean_ratio") 
  
  z <- 
    y %>% 
    mutate(NTC_normalised = NTC/NTC * 100,
           siCNOT3_normalised = siCNOT3/NTC * 100,
           siRPS25_normalised = siRPS25/NTC * 100) %>% 
    select(c(-NTC, -siCNOT3, -siRPS25)) %>% 
    gather(key = "condition", value = "percent", NTC_normalised: siRPS25_normalised) %>% 
    mutate(condition = str_remove(condition, pattern = "_normalised"))
  
  return(z)
  
}

final_data <- do.call("rbind", lapply(test_split, percentage_normalisation)) 

final_data$condition <- factor(final_data$condition, levels = c("NTC", "siCNOT3", "siRPS25"))

final_data$biorep <- factor(final_data$biorep, 
                            levels = c("1", "2", "3"))

final_data_bars <- final_data %>% 
  group_by(condition, signal_peptide, sample) %>% 
  summarise(mean_percent = mean(percent))

media <- 
ggplot(data = final_data %>% 
         filter(sample == "media"), 
       aes(x = condition, y = percent)) +
  #geom_boxplot() +
  geom_bar(stat = "summary", aes(group = "biorep"), alpha = 0.5, width = 0.5, show.legend = F) +
  geom_errorbar(stat='summary', aes(group = "biorep"), width=.2) +
  geom_point(position = position_jitter(width = 0.25)) +
  facet_wrap(~signal_peptide, nrow = 3) +
  ggtitle("Media") +
  scale_y_continuous(limits = c(0,250), breaks = seq(0, 220, by = 20)) +
  publication_theme() +
  scale_color_npg() +
  theme(legend.direction = "vertical",
        legend.position = "right")

lysate <- 
ggplot(data = final_data %>% 
         filter(sample == "lysate"), 
       aes(x = condition, y = percent)) +
  #geom_boxplot() +
  geom_bar(stat = "summary", aes(group = "biorep"), alpha = 0.5, width = 0.5, show.legend = F) +
  geom_errorbar(stat='summary', aes(group = "biorep"), width=.2) +
  geom_point(position = position_jitter(width = 0.25)) +
  facet_wrap(~signal_peptide, nrow = 3) +
  ggtitle("Lysate") +
  scale_y_continuous(limits = c(0,350), breaks = seq(0, 300, by = 50)) +
  publication_theme() +
  scale_color_npg() +
  theme(legend.direction = "vertical",
        legend.position = "right")

png("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/constructs_media_KD.png", res = 300, height = 2200, width = 1250)
print(media)
dev.off()

png("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/constructs_lysate_KD.png", res = 300, height = 2200, width = 1250)
print(lysate)
dev.off()

final_gplot <- 
  ggplot(data = final_data %>% 
           filter(sample == "media"), 
         aes(x = condition, y = percent)) +
  #geom_boxplot() +
  geom_bar(stat = "summary", aes(group = "biorep"), alpha = 0.5, width = 0.5, show.legend = F) +
  geom_errorbar(stat='summary', aes(group = "biorep"), width=.2) +
  geom_point() +
  facet_wrap(~signal_peptide, nrow = 3) +
  #scale_y_continuous(breaks = seq(0, 250, by = 20)) +
  publication_theme() +
  scale_color_npg() +
  theme(legend.direction = "vertical",
        legend.position = "right")

TukeyHSD(aov(percent~condition,
final_data %>% 
  filter(sample == "media") %>% 
  filter(signal_peptide == "SHH")))

TukeyHSD(aov(percent~condition,
             final_data %>% 
               filter(sample == "media") %>% 
               filter(signal_peptide == "COL")))

TukeyHSD(aov(percent~condition,
             final_data %>% 
               filter(sample == "lysate") %>% 
               filter(signal_peptide == "SHH")))

TukeyHSD(aov(percent~condition,
             final_data %>% 
               filter(sample == "lysate") %>% 
               filter(signal_peptide == "COL")))

tiff(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/reporter_assays/GAUSSIA+REDFF_TRANSFECTIONS/siCNOT3siRPS25.tiff", res = 300, height = 2000, width = 2000)
print(final_gplot)
dev.off()


ggplot(data = final_data %>% 
         filter(sample == "lysate",
                codon %in% c("CTG", "CTC", "TTA")), 
       aes(x = codon, y = percent, color = biorep)) +
  #geom_boxplot() +
  geom_bar(stat = "summary", aes(group = "biorep"), alpha = 0.5, width = 0.5, show.legend = F) +
  geom_errorbar(stat='summary', aes(group = "biorep"), width=.2) +
  geom_point(position = position_jitter(width = 0.25)) +
  facet_wrap(~signal_peptide, nrow = 3) +
  scale_y_continuous(breaks = seq(0, 120, by = 20)) +
  publication_theme() +
  scale_color_npg() +
  theme(legend.direction = "vertical",
        legend.position = "none")
