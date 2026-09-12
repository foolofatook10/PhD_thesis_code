library(data.table)
library(tidyverse)


data <- 
  read_csv("\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/COL8A2_RNA_data.csv") 

samples <- str_split(data$sample, " ")

rev_trans <- function(x, element){
  
  y <- x[element]
  return(y)
}



data$reverse_transcribed <- unlist(lapply(samples, rev_trans, element = 4))
data$construct <- unlist(lapply(samples, rev_trans, element = 1))
data$codon <- unlist(lapply(samples, rev_trans, element = 2))
data$rep <- unlist(lapply(samples, rev_trans, element = 3))

data <- data %>% rename(target = `Target Name`)

data <- 
data %>% mutate(codon = case_when(codon == "CTG" ~"CUG",
                          codon == "CTC" ~ "CUC",
                          codon == "TTA" ~ "UUA",
                          .default = codon))

data <- 
data %>% 
  filter(reverse_transcribed != "NRT") %>% 
  filter(target != "18S") %>% 
  select(-sample)

data$CT <- as.numeric(data$CT)

final <- 
  data %>% 
  select(reverse_transcribed, construct, codon, rep, target, CT) %>% 
  group_by(reverse_transcribed, construct, codon, rep, target) %>% 
  summarise(mean_CT = mean(CT)) %>% 
  filter(reverse_transcribed == "RT") %>% 
  spread(key = target, value = mean_CT) %>% 
  #mutate(GAU_18S = GAUSSIA - `18S`) %>%
  #mutate(REDFF_18S = `RED FIREFLY PRIMER2` - `18S`) %>% 
  mutate(GAU_REDFF = GAUSSIA - `RED FIREFLY PRIMER2`) %>% 
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
            plot.background = element_rect(colour = NA),
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
  filter(codon != "CUC") %>% 
  ggplot(aes(x = codon, y = normalised_expression)) +
  geom_bar(stat = "summary", aes(group = "rep"), alpha = 0.5, width = 0.5, show.legend = F)+
  geom_errorbar(stat='summary', aes(group = "rep"), width=.2) +
  geom_point(position = position_jitter(width = 0.1)) +
  ylab("mRNA level") +
  scale_x_discrete(labels = c("CTG", "TTA")) +
  ggtitle("COL8A2") +
  ylim(c(0, NA)) +
  publication_theme()

wilcox.test(normalised_expression ~ codon, data = final_normalised %>% filter(codon != "CUC"))

png(filename = "COL8A2_RNA.png", res = 300, height = 1000, width = 1000)
print(final_plot)
dev.off()

test <- final_normalised %>% filter(codon != "CUC")

one.way <- aov(normalised_expression~ codon, data = test)
summary(one.way)
