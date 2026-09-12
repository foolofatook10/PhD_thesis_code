library(tidyverse)
library(data.table)

CPMs <- 
  fread("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/CPMs/CDS_CPMs.csv",
      drop = "V1",
      header = T)


CPMs <- CPMs %>% filter(MD30 != "D30")

CPMs$condition <- paste0(CPMs$MD30, "_", CPMs$IP)

CPMs <- 
  CPMs %>% 
  select(replicate, condition, transcript, codon, summed_CPM_codon) %>% 
  dplyr::rename(CPM = summed_CPM_codon) %>% 
  mutate(CPM = (CPM/3))


# read in features
features <- fread("\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/features.csv",
                  drop = "V1",
                  header = T)

CNOT3_transcripts <- features %>% 
  filter(Gene == "CNOT3") %>% 
  pull(var = ENST)

CNOT3 <- CPMs %>% filter(transcript %in% CNOT3_transcripts)

#create themes----
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
            legend.key.size= unit(1, "cm"),
            legend.margin = unit(0.5, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

# RColorBrewer::brewer.pal(8,"Set2")
# 
# "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF"
# "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494" "#B3B3B3"

data <- 
  CNOT3 %>% 
  group_by(condition, transcript, codon) %>% 
  summarise(mean_CPM = mean(CPM))

data$condition <- factor(data$condition, 
                         levels = c("M_Tot", "M_CNOT3"),
                         labels = c("Tot", "CNOT3"))

CNOT3_reads <- 
  data %>% 
  ggplot(aes(x = codon, y = mean_CPM, colour = condition)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#FFC20A", "#0C7BDC")) +
  ggtitle("CNOT3") +
  xlab("Codon") +
  ylab("CPM") +
  geom_vline(xintercept = 450, linetype="dashed", size = 0.75) +
  geom_vline(xintercept = 753, linetype="dashed", size = 0.75) +
  publication_theme()


tiff(filename = "\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/plots/CNOT3_RPFs.tiff",
     height = 1500,
     width = 2000,
     res = 300)
print(CNOT3_reads)
dev.off()
