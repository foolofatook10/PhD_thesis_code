library(data.table)
library(tidyverse)

CPMs <- fread("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/CDS_CPMs.csv")

CPMs <- 
CPMs %>% 
  filter(MD30 == "M") %>% 
  select(replicate, IP, transcript, codon, summed_CPM_codon) %>% 
  rename(CPM = summed_CPM_codon)

CNOT3_trans <- CPMs$transcript %>% unique()

deltas <- 
  CPMs %>% 
  spread(key = IP, value = CPM) %>% 
  mutate(delta = CNOT3 - Tot) %>% 
  group_by(transcript, codon) %>% 
  summarise(mean_delta = mean(delta))


features <- read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/features.csv")
tmhmm <- read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/extracted_TMs.csv")

SignalP_transcripts <- features %>% filter(signal_sequence == T) %>% pull(ENST)

type1 <- 
tmhmm %>% 
  filter(class == "type1") %>% 
  select(ENST, class) %>% 
  unique() %>% 
  pull(ENST)

type1s_corroborated <- intersect(SignalP_transcripts, type1)

sponly <-   
  tmhmm %>% 
  filter(class == "sponly") %>% 
  select(ENST, class) %>% 
  unique() %>% 
  pull(ENST)

sponlys_corroborated <- intersect(SignalP_transcripts, sponly)

unknown <- 
tmhmm %>% 
  filter(class == "unknown") %>% 
  select(ENST, class) %>% 
  unique() %>% 
  pull(ENST)

intersect(unknown[unknown %in% CNOT3_trans], SignalP_transcripts)

unknowns_corroborated <- intersect(SignalP_transcripts, unknown)

transcript_comparison <- 
ggvenn::ggvenn(
list(SignalP = SignalP_transcripts[SignalP_transcripts %in% CNOT3_trans],
     `Tmhmm Type1` = type1[type1 %in% CNOT3_trans],
     `Tmhmm SP` = sponly[sponly %in% CNOT3_trans]),
show_percentage = F)

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/SiPTmhmm_comparison_C3transcripts.png", res = 300, height = 1250, width = 1250)
print(transcript_comparison)
dev.off()

data_prep <- 
deltas %>% 
  mutate(classification = case_when(transcript %in% sponlys_corroborated ~ "SP",
                                    transcript %in% type1s_corroborated ~ "T1",
                                    transcript %in% unknowns_corroborated ~ "unknown",
                                    transcript %in% SignalP_transcripts & !transcript %in% sponlys_corroborated & !transcript %in% type1s_corroborated & !transcript %in% unknowns_corroborated ~"SigP",
                                    .default = "noSP")) %>% 
  filter(codon <= 300) 

boundaries <- 
data_prep %>%   
group_by(classification, codon) %>% 
  summarise(delta = mean(mean_delta), 
            sd_delta = sd(mean_delta)) %>% 
  mutate(upper = delta + (2* sd_delta),
         lower = delta - (2 * sd_delta)) %>% 
  select(classification, codon, upper, lower)

data <- 
data_prep %>% 
  inner_join(boundaries, by = c("classification", "codon")) %>% 
  filter(mean_delta > lower & mean_delta < upper) %>% 
  #filter(codon == 50 & transcript == "ENST00000265643.4")
  group_by(classification, codon) %>% 
  summarise(mean_delta = mean(mean_delta)) 

# data <- 
# deltas %>% 
#   mutate(classification = case_when(transcript %in% sponlys_corroborated ~ "SP",
#                                     transcript %in% type1s_corroborated ~ "T1",
#                                     transcript %in% unknowns_corroborated ~ "unknown",
#                                     transcript %in% SignalP_transcripts & !transcript %in% sponlys_corroborated & !transcript %in% type1s_corroborated & !transcript %in% unknowns_corroborated ~"SigP",
#                                     .default = "noSP")) %>% 
#   filter(codon <= 300) %>% 
#   filter(!transcript %in% c("ENST00000316448.10", "ENST00000300289.10", "ENST00000265643.4","ENST00000381192.10", "ENST00000303575.9", "ENST00000457354.7", "ENST00000300026.4",  "ENST00000316448.10" )) %>% 
#   group_by(classification, codon) %>% 
#   summarise(mean_delta = mean(mean_delta)) 

deltas_extreme_50 <- 
deltas %>% 
  filter(codon == 50) %>% 
  filter(transcript %in% type1s_corroborated) %>%  
  arrange(mean_delta) %>%
  head(n = 4) %>% 
  pull(transcript)

deltas %>% 
  filter(transcript %in% deltas_extreme_50) %>%
  ggplot(aes(x = codon, y = mean_delta)) +
  geom_line() +
  facet_wrap(~transcript)

deltas %>% 
  filter(codon == 50) %>% 
  filter(transcript %in% sponlys_corroborated) %>%  
  arrange(mean_delta) %>%
  head(n = 2) %>% 
  pull(transcript)

  


data %>% filter(classification == "SigP")

# read in publication theme
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
            legend.text = element_text(size = 12),
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
  

subtype_breakdown <- 
ggplot(data %>% filter(classification %in% c("noSP", "SigP")), aes(x = codon, y = mean_delta, colour = classification)) +
  geom_line() +
  geom_point() +
  #geom_smooth() +
  scale_color_brewer(palette = "Set1") +
  ylab("Delta CPM (CNOT3 - Tot)") +
  coord_cartesian(xlim = c(0,300)) +
  publication_theme() +
  theme(legend.title = element_blank())
#, labels = c("Cytosolic", "N-term", "T1"))

ggplot(data %>% filter(classification %in% c("noSP", "T1")), aes(x = codon, y = mean_delta, colour = classification)) +
  geom_line() +
  geom_point() +
  #geom_smooth() +
  scale_color_brewer(palette = "Set1") +
  ylab("Delta CPM (CNOT3 - Tot)") +
  coord_cartesian(xlim = c(0,300)) +
  publication_theme() +
  theme(legend.title = element_blank())
#, labels = c("Cytosolic", "N-term", "T1"))

ggplot(data %>% filter(classification %in% c("noSP", "SP")), aes(x = codon, y = mean_delta, colour = classification)) +
  geom_line() +
  geom_point() +
  #geom_smooth() +
  scale_color_brewer(palette = "Set1") +
  ylab("Delta CPM (CNOT3 - Tot)") +
  coord_cartesian(xlim = c(0,300)) +
  publication_theme() +
  theme(legend.title = element_blank())
#, labels = c("Cytosolic", "N-term", "T1"))

ggplot(data %>% filter(classification %in% c("noSP", "SigP")), aes(x = codon, y = mean_delta, colour = classification)) +
  geom_line() +
  geom_point() +
  #geom_smooth() +
  scale_color_brewer(palette = "Set1") +
  ylab("Delta CPM (CNOT3 - Tot)") +
  coord_cartesian(xlim = c(0,300)) +
  publication_theme() +
  theme(legend.title = element_blank())
#, labels = c("Cytosolic", "N-term", "T1"))

CNOT3_selective_200 <- 
ggplot(data %>% filter(classification %in% c("noSP", "SP", "T1")), aes(x = codon, y = mean_delta, colour = classification)) +
  geom_line() +
  #geom_point() +
  #geom_smooth() +
  scale_color_brewer(palette = "Set1", labels = c("Cytosolic", "N-term", "T1")) +
  ylab("Delta CPM (CNOT3 - Tot)") +
  xlab("Codon position") +
  coord_cartesian(xlim = c(0,200)) +
  ggtitle("CNOT3 Selective") +
  publication_theme() +
  theme(legend.title = element_blank())
#, labels = c("Cytosolic", "N-term", "T1"))

CNOT3_selective_50 <- 
ggplot(data %>% filter(classification %in% c("noSP", "SP", "T1")), aes(x = codon, y = mean_delta, colour = classification)) +
  geom_line() +
  geom_point() +
  #geom_smooth() +
  scale_color_brewer(palette = "Set1", labels = c("Cytosolic", "N-term", "T1")) +
  ylab("Delta CPM (CNOT3 - Tot)") +
  xlab("Codon position") +
  ggtitle("CNOT Selective") +
  coord_cartesian(xlim = c(0,50)) +
  publication_theme() +
  theme(legend.title = element_blank())
#, labels = c("Cytosolic", "N-term", "T1"))

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_deltas_SP_subtypes_200.png",
    height = 1000, width = 1000, res = 300)
print(CNOT3_selective_200)
dev.off()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_deltas_SP_subtypes_50.png",
    height = 1000, width = 1000, res = 300)
print(CNOT3_selective_50)
dev.off()


DDX6_CPMs <- fread("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/riboseq_datasets/DDX6_CDS_CPMs.csv")

DDX6_mean_deltas <- 
DDX6_CPMs %>% 
  spread(key = condition, value = CPM) %>% 
  mutate(delta = DDX6- Total) %>% 
  group_by(transcript, codon) %>% 
  summarise(mean_delta = mean(delta)) %>% 
  filter(!transcript %in% c("ENST00000316448.10", "ENST00000300289.10", "ENST00000265643.4","ENST00000381192.10", "ENST00000303575.9", "ENST00000457354.7", "ENST00000300026.4",  "ENST00000316448.10" )) 

DDX6_transcripts <- unique(DDX6_mean_deltas$transcript)

DDX6_prep <- 
DDX6_mean_deltas %>% 
  mutate(classification = case_when(transcript %in% sponlys_corroborated ~ "SP",
                                                     transcript %in% type1s_corroborated ~ "T1",
                                                     transcript %in% unknowns_corroborated ~ "unknown",
                                    transcript %in% SignalP_transcripts & !transcript %in% sponlys_corroborated & !transcript %in% type1s_corroborated & !transcript %in% unknowns_corroborated ~"SigP",
                                                     .default = "noSP")) %>% 
  filter(codon <= 300)

DDX6_boundaries <- 
DDX6_prep %>% 
group_by(classification, codon) %>% 
  summarise(delta = mean(mean_delta),
            sd_delta = sd(mean_delta)) %>% 
  mutate(upper = delta + (1* sd_delta),
         lower = delta - (1*sd_delta))

DDX6_final <- 
DDX6_prep %>% 
  inner_join(DDX6_boundaries, by = c("classification", "codon")) %>% 
  filter(mean_delta > lower & mean_delta < upper)%>% 
  group_by(classification, codon) %>% 
  summarise(mean_delta = mean(delta))



DDX6_selective_50 <- 
ggplot(DDX6_final %>% filter(!classification %in% c("SigP", "unknown")), aes(x = codon, y = mean_delta, colour = classification)) +
  geom_line() +
  geom_point() +
  #geom_smooth() +
  scale_color_brewer(palette = "Set1", labels = c("Cytosolic", "N-term", "T1")) +
  ylab("Delta CPM (DDX6-Tot)") +
  xlab("Codon Position") +
  ggtitle("DDX6 Selective") +
  coord_cartesian(xlim = c(0,50)) +
  publication_theme() +
  theme(legend.title = element_blank())

DDX6_selective_200 <- 
ggplot(DDX6_final %>% filter(!classification %in% c("SigP", "unknown")), aes(x = codon, y = mean_delta, colour = classification)) +
  geom_line() +
  #geom_point() +
  #geom_smooth() +
  scale_color_brewer(palette = "Set1", labels = c("Cytosolic", "N-term", "T1")) +
  ylab("Delta CPM (DDX6-Tot)") +
  xlab("Codon Position") +
  ggtitle("DDX6 Selective") +
  coord_cartesian(xlim = c(0,200)) +
  publication_theme() +
  theme(legend.title = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/DDX6_deltas_SP_subtypes_50.png",
    height = 1000, width = 1000, res = 300)
print(DDX6_selective_50)
dev.off()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/DDX6_deltas_SP_subtypes_200.png",
    height = 1000, width = 1000, res = 300)
print(DDX6_selective_200)
dev.off()

fread("")
