library(tidyverse)
library(data.table)

RColorBrewer::brewer.pal(5, "Set1")
all_replicates_colours <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
replicates_1_4_5 <- c("#E41A1C", "#984EA3", "#FF7F00")

path2data <- "C:/Users/Jettles/OneDrive - University of Glasgow/RiboseqKD_quantification.csv"

data <- read_csv(path2data)

data %>% filter(Sample == "TUB")

data_summarised <- 
data %>% 
  group_by(Replicate, Sample, Band, Lane) %>% 
  summarise(Intensity_summed = sum(Intensity)) %>% 
  select(-Lane)

data_percentages <- 
data_summarised %>% 
  spread(key = Band, value = Intensity_summed) %>% 
  mutate(eS25 = if_else(is.na(eS25), 0, eS25)) %>% 
  mutate(CNOT1_normalised = CNOT1/TUB,
         CNOT3_normalised = CNOT3/TUB,
         CNOT7_normalised = CNOT7/TUB,
         eS25_normalised = eS25/TUB) %>% 
  select(Replicate, Sample, CNOT1_normalised, CNOT3_normalised, CNOT7_normalised, eS25_normalised) %>% 
  gather(key = band_value, value = normalised_intensity, CNOT1_normalised:eS25_normalised) %>% 
  spread(key = Sample,, value = normalised_intensity) %>% 
  mutate(id_tag = paste0(Replicate, "_", band_value)) %>% 
  filter(id_tag != "REP4_eS25_normalised" & id_tag != "REP5_eS25_normalised") %>% 
  select(-id_tag) %>% 
  mutate(NTC_normalised_perc = NTC/NTC * 100,
         CNOT3_normalised_perc = siCNOT3/NTC * 100,
         sieS25_normalised_perc = sieS25/NTC * 100)

long_format <- 
data_percentages %>% 
  select(-NTC, -siCNOT3, -sieS25) %>% 
  dplyr::rename(NTC = NTC_normalised_perc,
                siCNOT3 = CNOT3_normalised_perc,
                sieS25 = sieS25_normalised_perc) %>% 
  mutate(band_value = str_remove_all(band_value, "_normalised")) %>% 
  gather(key = condition, value = Intensity, NTC:sieS25)

long_format %>% 
  ggplot(aes(x = band_value, y = Intensity, fill = condition)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~Replicate)


###
unormalised_data <- 
data_summarised %>% 
  spread(key = Sample, value = Intensity_summed) %>% 
  mutate(NTC_perc = NTC/NTC * 100,
         siCNOT3_perc = siCNOT3/NTC * 100,
         sieS25_perc = sieS25/NTC * 100) %>% 
  select(-NTC, -siCNOT3, -sieS25) %>% 
  dplyr::rename(NTC = NTC_perc,
                siCNOT3 = siCNOT3_perc,
                sieS25 = sieS25_perc) %>% 
  gather(key = condition, value = intensity, NTC:sieS25) %>% 
  mutate(ID = paste0(Replicate, "_", Band)) %>% 
  filter(ID != "REP2_TUB" & ID != "REP3_TUB") %>% 
  select(-ID) %>% 
  mutate(Band = case_when(Band == "TUB" ~ "Tubulin", .default = Band)) %>% 
  ggplot(aes(x = condition, y = intensity, colour = Replicate)) +
  #geom_boxplot() +
  #geom_point(position = position_jitter(width = 0.25)) +
  facet_wrap(~Band, nrow = 1, scale = "free_y") +
  scale_color_manual(values = all_replicates_colours) +
  geom_point(position = position_jitter(width = 0.25, height = 0, seed = 117)) +
  stat_summary(
    aes(group = condition),
    fun = mean,
    geom = "bar",
    fill = "grey70",
    colour = "black",
    alpha = 0.5,
    width = 0.6
  ) +
  ylab("Unnormalised intensity (%)") +
  #scale_y_break(c(175, 350)) +
  scale_y_continuous(limits = c(0,NA), breaks = c(0,25,50,75,100,150,200,300)) +
  publication_theme() +
  theme(axis.title.x = element_blank())

png(file.path(save, "Riboseq_KDs_unnormalised_data.png"), res = 300, width = 2500, height = 1250)
print(unormalised_data)
dev.off()

library(cowplot)
png(file.path(save, "Riboseq_KDs_unnormalised_data_legend.png"), res = 300, width = 1750, height = 100)
ggdraw(get_legend(unormalised_data))
dev.off()

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

normalised_data_for_gplot <- 
data_summarised %>% 
  spread(key = Band, value = Intensity_summed) %>% 
  mutate(eS25 = if_else(is.na(eS25), 0, eS25)) %>% 
  mutate(CNOT1_normalised = CNOT1/TUB,
         CNOT3_normalised = CNOT3/TUB,
         CNOT7_normalised = CNOT7/TUB,
         eS25_normalised = eS25/TUB) %>% 
  select(Replicate, Sample, CNOT1_normalised, CNOT3_normalised, CNOT7_normalised, eS25_normalised) %>% 
  gather(key = band_value, value = normalised_intensity, CNOT1_normalised:eS25_normalised) %>% 
  filter(Replicate %in% c("REP1", "REP4", "REP5")) %>% 
  mutate(band_value = str_remove_all(band_value, "_normalised")) %>% 
  spread(Sample, normalised_intensity) %>% 
  filter(!c(Replicate == "REP4" & band_value == "eS25")) %>% 
  filter(!c(Replicate == "REP5" & band_value == "eS25")) %>% 
  mutate(NTC_perc = NTC/NTC * 100,
         siCNOT3_perc = siCNOT3/NTC * 100,
         sieS25_perc = sieS25/NTC * 100) %>% 
  select(-NTC, -siCNOT3, -sieS25) %>% 
  gather(key = condition, value = intensity_normalised, NTC_perc:sieS25_perc) %>% 
  mutate(condition = str_remove_all(condition, "_perc")) 
  

gplot_normalised_data <- 
normalised_data_for_gplot %>% 
ggplot(aes(x = condition, y = intensity_normalised, colour = Replicate)) +
  geom_point(position = position_jitter(width = 0.25, height = 0, seed = 117)) +
  stat_summary(
    aes(group = condition),
    fun = mean,
    geom = "bar",
    fill = "grey70",
    colour = "black",
    alpha = 0.5,
    width = 0.6
  ) +
  facet_wrap(~band_value, scale = "free_y", nrow = 1) +
  scale_y_continuous(limits = c(0,NA)) +
  scale_color_manual(values = replicates_1_4_5) +
  ylab("Normalised intensity (%)")  +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

normalised_data_for_gplot %>%
  filter(condition %in% c("NTC", "siCNOT3", "sieS25")) %>%
  select(band_value, Replicate, condition, intensity_normalised) %>%
  pivot_wider(
    names_from = condition,
    values_from = intensity_normalised
  ) %>%
  filter(band_value != "eS25") %>% 
  group_by(band_value) %>%
  summarise(
    p_siCNOT3 = t.test(siCNOT3, NTC, paired = T)$p.value,
    p_sieS25  = t.test(sieS25, NTC, paired = T)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_siCNOT3 = p.adjust(p_siCNOT3, method = "BH"),
    p_adj_sieS25  = p.adjust(p_sieS25, method = "BH")
  )

#save = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures"
png(file.path(save, "riboseq_KDs_quantification_normalised.png"), res = 300, width = 2500, height = 1250)
print(gplot_normalised_data)
dev.off()
