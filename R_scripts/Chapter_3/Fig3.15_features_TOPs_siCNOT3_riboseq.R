library(tidyverse)
library(data.table)

siCNOT3 <- 
  read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/Analysis/DESeq2_output/FINAL_siCNOT3_merged_DESeq2.csv")

features <- 
  read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/features.csv") %>% 
  select(-`...1`)

most_abundant_transcript <- read_csv("\\\\data.beatson.gla.ac.uk/data/JETTLES/CNOT3_Riboseq/bioinformatics/Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv")

siCNOT3 <- 
  siCNOT3 %>% inner_join(most_abundant_transcript, by = c("gene", "gene_sym"))

joined <- 
  inner_join(siCNOT3,
             features,
             by = c("gene" = "ENSG",
                    "gene_sym" = "Gene",
                    "transcript" = "ENST"))

joined$RPFs_group <- factor(joined$RPFs_group,
                            levels = c("both down",
                                       "both up",
                                       "Totals down",
                                       "Totals up",
                                       "RPFs down",
                                       "RPFs up", 
                                       "no change", "NS"),
                            labels = c("Both\ndown",
                                       "Both\nup",
                                       "Totals\ndown",
                                       "Totals\nup",
                                       "RPFs\ndown",
                                       "RPFs\nup", 
                                       "No\nchange", "NS"))

my_colours <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", "#C8C8C8")
my_colours2 <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", "#C8C8C8")

# read in publication theme
publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            #axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2, size = 20),
            axis.title.x = element_text(vjust = -0.2, size = 20),
            axis.text.x = element_text(size = 20, color = "black"),
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
            legend.margin = unit(1, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

joined %>% 
  ggplot(aes(x = totals_log2FC,y = RPFs_log2FC, colour= RPFs_group)) +
  geom_point() +
  scale_color_manual(values = my_colours) +
  ggpubr::stat_cor(data = joined, aes(x = totals_log2FC,y = RPFs_log2FC), inherit.aes = F)

gplot <- 
  ggplot(joined, aes(x = RPFs_group, y = GC3_cont * 100, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = 0.4) +
  scale_y_continuous(limits = c(0,NA)) +
  ylab("GC3%") +
  scale_fill_manual(values = my_colours) +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        
        legend.position = "none")

TukeyHSD(aov(data = joined, formula = GC3_cont~RPFs_group))



png(filename = "C:/Users/Jettles/OneDrive - University of Glasgow/Documents/GC3_percentage.png",
    res = 300,
    width = 2250,
    height = 1000)
print(gplot)
dev.off()

TOPs <- 
  read_csv("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/Gene_lists/TOPs/pnas.1912864117.sd01.csv", col_names = F) %>% pull(X1)

joined %>% filter(RPFs_group == "RPFs\nup")
joined %>% filter(RPFs_group == "RPFs\nup" & gene_sym %in% TOPs)

nottops <- 
joined %>%
  filter(!gene_sym %in% TOPs)

tops <- 
  joined %>%
  filter(gene_sym %in% TOPs)

TOP_gplot_final <- 
ggplot() +
  geom_point(data = nottops, 
             aes(x = totals_log2FC,y = RPFs_log2FC, colour= RPFs_group),
             colour = "#C8C8C8") +
  geom_point(data = tops,
             aes(x = totals_log2FC,y = RPFs_log2FC, colour= RPFs_group),
             colour = "red")+
  xlab("Cytoplasmic RNA log2FC") +
  ylab("RPFs log2FC") +
  ggtitle("TOP motif\ncontaining mRNAs") +
  publication_theme()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/siCNOT3_riboseq_TOP_mRNAs.png", res = 300, height = 1250, width = 1250)
print(TOP_gplot_final)
dev.off()

gplot <- 
joined %>% 
  mutate(TOP_motif = case_when(gene_sym %in% TOPs ~ "TOP mRNA", .default = "TOP absent mRNA")) %>% 
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = TOP_motif)) +
  geom_point() +
  scale_color_manual(values = c("#C8C8C8", "red")) +
  publication_theme() +
  theme(legend.position = "right",
        legend.direction = "vertical") +
  labs(color = "TOP presence")

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/TOP_legend.png", res = 300, height = 500, width = 750)
ggdraw(cowplot::get_legend(gplot))
dev.off()

ggplot() +
  geom_point(data = nottops, 
             aes(x = totals_log2FC,y = RPFs_log2FC, colour= RPFs_group),
             colour = "#C8C8C8") +
  geom_point(data = tops,
             aes(x = totals_log2FC,y = RPFs_log2FC, colour= RPFs_group),
             colour = "red")+
  xlab("Cytoplasmic RNA log2FC") +
  ylab("RPFs log2FC") +
  ggtitle("TOP motif\ncontaining mRNAs") +
  publication_theme()

ggvenn::ggvenn(
list(TE_up=
joined %>% 
  #filter(TE_group == "TE up") %>% 
  pull(gene_sym),
TOPs = TOPs))

joined %>% 
  filter(TE_group != "TE up") %>% 
  pull(gene_sym)
