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
  read_csv("humanised_constructs_final2.csv") %>% 
  select(signal_peptide:biorep) %>% na.omit()

split(data, f = data$signal_peptide)


raw <- 
  data %>%
  group_by(signal_peptide, codon, sample, channel, biorep) %>%
  summarise(mean_value = mean(value))


raw_lyst <- split(raw, f = raw$signal_peptide)

for(i in 1:length(raw_lyst)){
  
  dat <- raw_lyst[[i]]
  
  dat <- dat %>% filter(codon %in% c("CTG", "CTC", "TTA"))
  
  dat$codon <- factor(dat$codon, levels = c("CTG", "CTC", "TTA"))
  
  dat$biorep <- as.character(dat$biorep)
  
  media <- 
    dat %>%
    filter(sample == "media") %>%
    ggplot(aes(x = codon, y = mean_value, colour = biorep)) +
    geom_point() +
    ggtitle(paste(dat$signal_peptide[1]," Media")) +
    facet_wrap(~channel, scales = "free") + 
    publication_theme()
  
  print(media)
  
  lysate <- 
    dat %>%
    filter(sample == "lysate") %>%
    ggplot(aes(x = codon, y = mean_value, colour = biorep)) +
    geom_point() +
    ggtitle(paste(dat$signal_peptide[1]," Lysate")) +
    facet_wrap(~channel, scales = "free") + 
    publication_theme()
  
  
  print(lysate)
  
}

normalised_data <- 
  data %>%
  mutate(ID = paste0(channel, "_", sample)) %>% 
  select(-sample, -channel) %>% 
  spread(key = ID, value = value) %>%
  mutate(media_normalised = gaussia_media/red_firefly_lysate,
         lysate_normalised = gaussia_lysate/red_firefly_lysate) %>% 
  select(signal_peptide, codon, techrep, biorep, media_normalised, lysate_normalised) %>% 
  gather(key = sample, value = normalised_value, media_normalised:lysate_normalised) %>% 
  mutate(sample = str_remove(sample, "_normalised")) %>% 
  filter(!codon %in% c("CTT", "TTG"))


normalised_data$codon <- factor(normalised_data$codon, levels = c("CTG", "CTC", "TTA"))

#summary(normalised_data)

normalised_data$biorep <- factor(normalised_data$biorep, 
                                 levels = c("1", "2", "3", "4", "5"))

normalised_data %>% 
  filter(signal_peptide == "COL8A2") %>% 
  #filter(sample == "media") %>% 
  ggplot(aes(x = codon, y = normalised_value)) +
  geom_point(alpha = 0.3) + 
  #position = position_jitter()) +
  ylab("Normalised luciferase") +
  facet_grid(biorep~sample) +
  #ylim(0,25) +
  publication_theme() +
  theme(axis.title.x = element_blank())

test <- 
  normalised_data %>% 
  group_by(signal_peptide, codon, sample, biorep) %>% 
  summarise(mean_ratio = mean(normalised_value, na.rm = T))

test_split <- split(test, f = test$signal_peptide)

percentage_normalisation <- function(x){
  
  y <- x  %>% spread(key = "codon", value = "mean_ratio") 
  
  z <- 
    y %>% 
    mutate(CTG_normalised = CTG/CTG * 100,
           CTC_normalised = CTC/CTG * 100,
           #CTA_normalised = CTA/CTG * 100,
           TTA_normalised = TTA/CTG * 100) %>% 
    select(c(-CTG, -CTC, -TTA)) %>% 
    gather(key = "codon", value = "percent", CTG_normalised: TTA_normalised) %>% 
    mutate(codon = str_remove(codon, pattern = "_normalised"))
  
  return(z)
  
}

final_data <- do.call("rbind", lapply(test_split, percentage_normalisation)) 

final_data$codon <- factor(final_data$codon, levels = c("CTG", "CTC", "TTA"))

final_data$biorep <- factor(final_data$biorep, 
                            levels = c("1", "2", "3", "4", "5"))

final_data_bars <- final_data %>% 
  group_by(signal_peptide, sample, codon) %>% 
  summarise(mean_percent = mean(percent))

ggplot(data = final_data %>% 
         filter(sample == "media"), 
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
        legend.position = "right")

final_gplot <- 
  ggplot(data = final_data %>% 
           filter(sample == "media",
                  codon %in% c("CTG", "CTC", "TTA")), 
         aes(x = codon, y = percent, color = biorep)) +
  #geom_boxplot() +
  geom_bar(stat = "summary", aes(group = "biorep"), alpha = 0.5, width = 0.5, show.legend = F) +
  geom_errorbar(stat='summary', aes(group = "biorep"), width=.2) +
  geom_point(position = position_jitter(width = 0.25)) +
  facet_wrap(~signal_peptide, nrow = 3) +
  #ylim(0,120) +
  scale_y_continuous(limits = c(0,150), breaks = seq(0, 150, by = 20)) +
  ggtitle("Media") +
  ylab("Gaussia activity (%)") +
  publication_theme() +
  scale_color_npg() +
  theme(legend.direction = "vertical",
        legend.position = "none")

COL8A2_lysate <- final_data %>% filter(signal_peptide == "COL8A2" & sample == "lysate")
COL8A2_media <- final_data %>% filter(signal_peptide == "COL8A2" & sample == "media")

SHH_lysate <- final_data %>% filter(signal_peptide == "SHH" & sample == "lysate")
SHH_media <- final_data %>% filter(signal_peptide == "SHH" & sample == "media")

TukeyHSD(aov(percent~codon ,data = SHH_media))
TukeyHSD(aov(percent~codon ,data = COL8A2_media))

TukeyHSD(aov(percent~codon ,data = SHH_lysate))
TukeyHSD(aov(percent~codon ,data = COL8A2_lysate))

tiff(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/reporter_assays/GAUSSIA+REDFF_TRANSFECTIONS/CTG_CTC_TTA.tiff", res = 300, height = 1500, width = 1000)
print(final_gplot)
dev.off()

lysate_gplot <- 
ggplot(data = final_data %>% 
         filter(sample == "lysate",
                codon %in% c("CTG", "CTC", "TTA")), 
       aes(x = codon, y = percent, colour = biorep)) +
  #geom_boxplot() +
  geom_bar(stat = "summary", aes(group = "biorep"), alpha = 0.5, width = 0.5, show.legend = F) +
  geom_errorbar(stat='summary', aes(group = "biorep"), width=.2) +
  geom_point(position = position_jitter(width = 0.25)) +
  facet_wrap(~signal_peptide, nrow = 3) +
  scale_y_continuous(limits = c(0,160),
                     breaks = seq(0, 160, by = 20)) +
  ggtitle("Lysate") +
  ylab("Gaussia activity (%)") +
  publication_theme() +
  scale_color_npg() +
  theme(legend.direction = "vertical",
        legend.position = "none")

tiff(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/reporter_assays/GAUSSIA+REDFF_TRANSFECTIONS/CTG_CTC_TTA_lysate.tiff", res = 300, height = 1500, width = 1000)
print(lysate_gplot)
dev.off()

