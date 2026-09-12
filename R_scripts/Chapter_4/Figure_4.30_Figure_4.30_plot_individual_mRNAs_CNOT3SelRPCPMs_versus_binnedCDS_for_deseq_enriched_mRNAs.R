library(tidyverse)

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

#read in common variables----
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/R_scripts/common_variables.R")

load(file.path(parent_dir, "Counts_files/R_objects/M_binned_list.Rdata"))

binned_data <- do.call("rbind", binned_list)

deseq <- 
read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/Analysis/DESeq2_output/RPFs_CNOT3_DEseq2_apeglm_LFC_shrinkage_M_28nts.csv")

enriched_transcripts <- 
deseq%>% 
  filter(padj < 0.05 & log2FoldChange > 1) %>% 
  filter(gene_sym != "CNOT3") %>% 
  pull(transcript)

enriched_transcripts_summarised <- 
binned_data %>% 
  filter(transcript %in% enriched_transcripts) %>% 
  inner_join(deseq %>% select(gene_sym, transcript), by = "transcript") %>% 
  group_by(gene_sym, transcript, bin, condition, region) %>% 
  summarise(mean_CPM = mean(binned_cpm),
            median_CPM = median(binned_cpm)) %>% 
  mutate(condition = case_when(condition == "M_CNOT3" ~ "CNOT3",
                               condition == "M_Tot" ~ "Tot")) %>% 
  mutate(bin = bin * 2)

split_data <- 
split(enriched_transcripts_summarised, enriched_transcripts_summarised$transcript)

plot_data <- 
function(x){
  
  x %>%
    filter(region == "CDS") %>% 
    ggplot(aes(x = bin, y = mean_CPM, colour = condition)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("#0C7BDC", "#FFC20A")) +
    ggtitle(unique(x$gene_sym)) +
    ylab("CPM") +
    xlab("CDS%") +
    publication_theme()
  
}

gplots <- lapply(split_data, plot_data)
gplots$ENST00000200453.6

for (i in 1:length(gplots)) {
  
  genesymbol <- unique(split_data[[i]]$gene_sym)
 
  
  png(paste0("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/binned_cpms/", genesymbol, ".png"),
      height = 1000, width = 1200, res = 300)
  print(gplots[[i]])
  dev.off()
}
