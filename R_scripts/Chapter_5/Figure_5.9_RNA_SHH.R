library(data.table)
library(tidyverse)


data <- 
  read_csv("\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/SHHREP2REP3RNA.csv") 

samples <- str_split(data$`Sample Name`, " ")

rev_trans <- function(x, element){
  
  y <- x[element]
  return(y)
}



data$reverse_transcribed <- unlist(lapply(samples, rev_trans, element = 1))
data$construct <- unlist(lapply(samples, rev_trans, element = 2))
data$codon <- unlist(lapply(samples, rev_trans, element = 3))
data$rep <- unlist(lapply(samples, rev_trans, element = 4))

data <- data %>% mutate(codon = case_when(codon == "CTG" ~"CUG",
                                          codon == "CTC" ~ "CUC",
                                          codon == "TTA" ~ "UUA",
                                          .default = codon)) %>% 
  mutate(techrep = rep(1:3, nrow(data)/3)) %>% 
  filter(`Target Name` != "28S") %>% 
  dplyr::select(-`Sample Name`) %>% 
  filter(CT != "Undetermined") %>% 
  rename(target = `Target Name`)

data$CT <- as.numeric(data$CT)

data %>% 
  group_by(target, reverse_transcribed, codon, rep) %>% 
  summarise(mean_CT = mean(CT)) %>% 
  ggplot(aes(x = reverse_transcribed, y = mean_CT, colour = rep)) +
  geom_point() +
  facet_grid(target~codon)






final <- 
data %>% 
  dplyr::select(reverse_transcribed, construct, codon, rep, target, CT) %>% 
  filter(target %in% c("GAUSSIA", "RED FIRFLY PRIMER 2")) %>% 
  group_by(reverse_transcribed, construct, codon, rep, target) %>% 
  summarise(mean_CT = mean(CT)) %>%  
  filter(reverse_transcribed == "RT") %>% 
  spread(key = target, value = mean_CT) %>% 
  #mutate(GAU_28S = GAUSSIA - `28S`) %>%
  #mutate(REDFF_28S = `RED FIRFLY PRIMER 2` - `28S`) %>% 
  mutate(GAU_REDFF = GAUSSIA - `RED FIRFLY PRIMER 2`) %>% 
  mutate(expression_GAUREDFF = 2^-GAU_REDFF)

final$codon <- factor(final$codon, levels = c("CUG", "CUC", "UUA"))

#Load theme
publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA, fill = NA),
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

my_expression <- paste0("\u0394", "\u0394", "CT")

final %>% 
  ggplot(aes(x = codon, y = expression_GAUREDFF, colour = rep)) +
  geom_point() +
  ggtitle(unique(final$construct)) +
  ylim(c(0,NA)) +
  ylab(expression("Gaussia/Red Firefly")) +
  publication_theme()



normalisation <- function(x){

  y <- x %>% filter(codon == "CUG") %>% pull(expression_GAUREDFF)
  
  z <- x %>% mutate(normalised_expression = expression_GAUREDFF/y * 100)
  
  return(z)
  
}  

final_normalised <- do.call("rbind", lapply(split(final, f = final$rep), normalisation))

final_plot <- 
final_normalised %>% 
  ggplot(aes(x = codon, y = normalised_expression)) +
  geom_bar(stat = "summary", aes(group = "rep"), alpha = 0.5, width = 0.5, show.legend = F)+
  geom_errorbar(stat='summary', aes(group = "rep"), width=.2) +
  geom_point(position = position_jitter(width = 0.1)) +
  ylab("mRNA level") +
  scale_x_discrete(labels = c("CTG", "CTC", "TTA")) +
  ggtitle("SHH") +
  ylim(c(0, NA)) +
  publication_theme()

TukeyHSD(aov(normalised_expression~codon, data = final_normalised))

png(filename = "SHH_RNA.png", res = 300, height = 1000, width = 1000)
print(final_plot)
dev.off()

one.way <- aov(normalised_expression ~ codon ,data = final_normalised)

TukeyHSD(one.way)
