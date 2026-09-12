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

# Load in TOP mRNAs
TOP_genes <- read_csv("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/Gene_lists/TOPs/pnas.1912864117.sd01.csv", col_names = F) %>% pull(X1)

DEseq2 <- DEseq2 %>% 
  mutate(TOP_mRNA = factor(case_when(gene_sym %in% TOP_genes~ "TOP",.default = "normal"),
                           levels = c("normal", "TOP")))

normal_genes <- DEseq2 %>% filter(TOP_mRNA == "normal")
TOP_genes_DE <- DEseq2 %>% filter(TOP_mRNA == "TOP")

TOP_volcano <- 
ggplot() +
  geom_point(data = normal_genes, aes(x = log2FoldChange, y = -log10(padj)), colour = "black") +
  geom_point(data = TOP_genes_DE, aes(x = log2FoldChange, y = -log10(padj)), colour = "red")  +
  xlab("Log2 Fold Change\n(CNOT3 IP'd - Total RPFs)") +
  ylab("-Log10(padj)") +
  publication_theme()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/TOP_volcano.png",
    height = 1250, width = 1250, res = 300)
print(TOP_volcano)
dev.off()
