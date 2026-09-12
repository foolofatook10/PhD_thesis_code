library(tidyverse)
# read in publication theme
publication_theme <- function(base_size=18, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = 22, hjust = 0.5),
            #text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = 20),
            axis.title.y = element_text(angle=90,vjust =2, size = 20),
            axis.title.x = element_text(vjust = -0.2, size = 20),
            axis.text.x = element_text(size = 18, color = "black"),
            axis.text.y = element_text(size = 18, color = "black"),
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
            legend.title = element_text(face="italic", size = 14),
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}
#Set directories
parent_dir = "\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL"
path2deseq = file.path(parent_dir,"Analysis/DESeq2_output/RPFs_CNOT3_DEseq2_apeglm_LFC_shrinkage_M.csv")
DEseq2 <- read_csv(file = path2deseq)
#Define 3 groups
DEseq2 <- DEseq2 %>% 
  mutate(group = factor(case_when(log2FoldChange <= 0 & padj < 0.05 ~ "downregulated",
                                  log2FoldChange >= 0 & padj < 0.05 ~ "upregulated",
                                  TRUE ~ "unchanged")))

tmhmm <- read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/TM_domains_FINAL.csv")
er_targeted_transcripts <- tmhmm %>% select(ENST, class) %>% unique()

tidied_data <- 
DEseq2 %>% 
  mutate(secretory_mRNAs = case_when(transcript %in% unique(er_targeted_transcripts$ENST) ~ "ER targeted",
                   .default = "Cytosolic")) %>% 
  select(transcript, group, secretory_mRNAs) 

data_for_plotting <- 
inner_join(
tidied_data %>% 
  group_by(group, secretory_mRNAs) %>% 
  summarise(n_mRNAs = n()),
tidied_data %>% 
  group_by(group) %>% 
  summarise(total_mRNAs = n()), 
by = "group") %>% 
  mutate(percentage_mRNAs = n_mRNAs/total_mRNAs * 100) %>% 
  filter(secretory_mRNAs == "ER targeted") 

data_for_plotting$group <- factor(data_for_plotting$group, 
                                  levels = c("downregulated", "unchanged", "upregulated"),
                                  labels = c("Dep", "Unch.", "Enr."))

SP_or_TMs <- 
data_for_plotting %>% 
  ggplot(aes(x = group, fill = group, y = percentage_mRNAs, label = n_mRNAs)) +
  geom_col() +
  geom_text(vjust = -0.2, fontface = "bold") +
  scale_fill_brewer(palette = "Dark2") +
  ylab("Proportion mRNAs\nwith SP/TM domain (%)") +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3SelRP_SPTM_enrichment.png",
    res = 300,
    height = 1500,
    width = 1250)  
print(SP_or_TMs)
dev.off()

inner_join(
    tidied_data %>% 
      group_by(group, secretory_mRNAs) %>% 
      summarise(n_mRNAs = n()),
    tidied_data %>% 
      group_by(group) %>% 
      summarise(total_mRNAs = n()), 
    by = "group") %>% 
  mutate(percentage_mRNAs = n_mRNAs/total_mRNAs * 100)

fisher_table <- matrix( c( 402, 890, 2494, 8277 ), nrow = 2, byrow = TRUE )
  
rownames(fisher_table) <- c("Upregulated", "Not_upregulated") 
colnames(fisher_table) <- c("ER_targeted", "Cytosolic")
  
fisher_result <- fisher.test( fisher_table, alternative = "greater" ) 
fisher_result
  