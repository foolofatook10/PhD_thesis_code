library(tidyverse)
library(data.table)

path = "\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/reporter_assays/GAUSSIA+REDFF_TRANSFECTIONS/Tranfections_from_mini_preps"

setwd(path)

data <- read_csv("Knockdown_RNA.csv")
colnames(data) <- c("Sample", "Target", "CT")

signal_peptide <- unlist(lapply(str_split(data$Sample, " "), function(x){x[1]}))
siRNA <- unlist(lapply(str_split(data$Sample, " "), function(x){x[2]}))
replicate <- unlist(lapply(str_split(data$Sample, " "), function(x){x[3]}))
reverse_transcription <- unlist(lapply(str_split(data$Sample, " "), function(x){x[4]}))
tech_rep <- rep(1:3,nrow(data)/3)

data <- 
data %>% 
  mutate(signal_peptide = signal_peptide,
         siRNA = siRNA,
         replicate = replicate,
         tech_rep = tech_rep,
         reverse_transcription = reverse_transcription) %>% 
  dplyr::select(signal_peptide, siRNA, replicate, tech_rep, reverse_transcription, Target, CT) %>% 
  filter(CT != "Undetermined")

data$CT <- as.numeric(data$CT) 

final <- 
data %>% 
  group_by(signal_peptide, siRNA, replicate, reverse_transcription, Target) %>% 
  summarise(mean_CT = mean(CT)) %>% 
  spread(key = Target, value = mean_CT) %>% 
  mutate(GAU_REDFF = GAUSSIA - `RED FIREFLY Primer2`) %>% 
  mutate(expression_GAUREDFF = 2^-GAU_REDFF)%>% 
  filter(reverse_transcription == "RT")

final$siRNA <- factor(final$siRNA, levels = c("NTC", "siCNOT3", "siRPS25"))

#Load theme
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
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


final %>% 
  ggplot(aes(x = siRNA, y = expression_GAUREDFF)) +
  geom_point() +
  facet_wrap(~replicate) +
  ylim(c(0,NA)) +
  ylab(expression("Gaussia/Red Firefly")) +
  publication_theme()

RNA <- split(final, final$replicate)

CT_normalisation <- function(x){
  
  y <- 
    x %>% 
    filter(siRNA == "NTC") %>% 
    pull(expression_GAUREDFF)
  
 z <- x %>% 
   mutate(normalised_expression = expression_GAUREDFF/y * 100)
 
 return(z)
}

final_normalised <- do.call("rbind",lapply(RNA, CT_normalisation))

means <- 
final_normalised %>% 
  group_by(signal_peptide, siRNA, reverse_transcription) %>% 
  summarise(mean_expression = mean(normalised_expression))

mRNA_level <- 
final_normalised %>% 
  ggplot(aes(x = siRNA, y = normalised_expression)) +
  geom_bar(stat = "summary", aes(group = "rep"), alpha = 0.5, width = 0.5, show.legend = F)+
  geom_errorbar(stat='summary', aes(group = "rep"), width=.2) +
  geom_point(position = position_jitter(width = 0.1)) +
  ylab("mRNA level") +
  scale_x_discrete(labels = c("NTC", "siCNOT3", "sieS25")) +
  ggtitle("SHH-CTG") +
  ylim(c(0, 150)) +
  publication_theme()

TukeyHSD(aov(normalised_expression~siRNA, final_normalised))

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/SHHCTG_mRNA_KD.png",
    res = 300, height = 1250, width = 1000)
print(mRNA_level)
dev.off()
