library(tidyverse)

my_colours <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", "#D3D3D3")

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
            axis.text.x = element_text(size = 17, color = "black"),
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
            legend.key.size= unit(1, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

get_stats <- function(formula_term){
  
  if(formula_term %in% c("fputr_length", "cds_length", "tputr_length", "miRNA_binding_sites")){
    
    temp_data <- data_joined %>% mutate(fputr_length = log10(fputr_length + 1),
                                        cds_length = log10(cds_length + 1),
                                        tputr_length = log10(tputr_length + 1),
                                        miRNA_binding_sites = log10((miRNA_binding_sites + 1)))
    
  } else{
    
    temp_data = data_joined
    
  }
  
  stats1 <- 
    (TukeyHSD(aov(data = temp_data, formula = as.formula(eval(parse(text = formula_term)) ~ RPFs_group))))$RPFs_group %>% 
    as.data.frame() %>% 
    #filter(`p adj` < 0.05) %>% 
    rownames_to_column("comparison")
  
  terms_to_filter <- c("Both\nup-Both\ndown", 
                       "Totals\nup-Totals\ndown",
                       "No\nchange-Both\ndown",
                       "No\nchange-Both\nup",
                       "RPFs\nup-RPFs\ndown",
                       "No\nchange-Totals\ndown",
                       "No\nchange-Totals\nup",
                       "No\nchange-RPFs\ndown",
                       "No\nchange-RPFs\nup")
  
  stats1 <- stats1 %>% filter(comparison %in% terms_to_filter)
  
  terms1 <- str_replace(unlist(lapply(str_split(stats1$comparison, "-"), function(x){x[[1]]})), "\n", " ")
  terms2 <- str_replace(unlist(lapply(str_split(stats1$comparison, "-"), function(x){x[[2]]})), "\n", " ")
  
  stats1$group1 <-  terms1
  stats1$group2 <- terms2
  
  final_stats <- 
    stats1 %>% select(-comparison) %>% 
    mutate(comparison = formula_term)
  
  final_stats$group1 <- factor(final_stats$group1, 
                               levels =c("Both up",
                                         "Totals up", 
                                         "RPFs up",
                                         "No change"))
  
  return(final_stats)
  
}

# define parent directory
parent_dir <- "\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq"

# read in merged data
data <- read_csv(file.path(parent_dir, "Analysis/DESeq2_output/FINAL_siCNOT3_merged_DESeq2.csv"))

most_abundant_transcript <- read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs_siCNOT3.csv")

data <- 
data %>% inner_join(most_abundant_transcript,
                    by = c("gene", "gene_sym"))

features <- read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/features.csv")

data_joined <- 
  data %>% 
  select(gene_sym, gene, transcript, RPFs_group, TE_group) %>% 
  inner_join(features,
             by = c("gene_sym" = "Gene",
                    "gene" = "ENSG",
                    "transcript" = "ENST")) %>% 
  select(-`...1`) %>% 
  filter(RPFs_group != "NS")

data_joined$RPFs_group %>% unique()

data_joined$RPFs_group <- 
  factor(data_joined$RPFs_group, 
         levels = c("both down", "both up",
                    "Totals down", "Totals up",
                    "RPFs down", "RPFs up",
                    "no change"),
         labels = c("Both\ndown", "Both\nup",
                    "Totals\ndown", "Totals\nup",
                    "RPFs\ndown", "RPFs\nup",
                    "No\nchange"))

summary(data_joined$RPFs_group)

RColorBrewer::brewer.pal.info

UTR_5P <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = GC_cont_fp * 100, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = 0.25) +
  ylab("GC%") +
  #ggtitle("5'UTR") +
  scale_y_continuous(limits = c(0,NA)) +
  publication_theme() +
  scale_fill_manual(values = my_colours) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

get_stats("GC3_cont")
get_stats("GC_cont_fp")
get_stats("GC_cont_cds")
get_stats("GC_cont_tp")
get_stats("AG_content")
get_stats("fputr_length")
get_stats("cds_length")
get_stats("tputr_length")
get_stats("miRNA_binding_sites")

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_Riboseq_UTR5.png",
    res = 300, height = 1000, width = 1750)
print(UTR_5P)
dev.off()

GC_CDS <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = GC_cont_cds * 100, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = 0.25) +
  ylab("GC%") +
  #ggtitle("CDS") +
  scale_y_continuous(limits = c(0,100)) +
  publication_theme() +
  #scale_y_continuous(limits = c(0,100)) +
  scale_fill_manual(values = my_colours) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_Riboseq_CDS.png",
    res = 300, height = 1000, width = 1750)
print(GC_CDS)
dev.off()

GC_3PUTR <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = GC_cont_tp * 100, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = 0.25) +
  ylab("GC%") +
  #ggtitle("3'UTR") +
  scale_fill_manual(values = my_colours) +
  publication_theme() +
  scale_y_continuous(limits = c(0,100)) +
  scale_colour_manual(values = colours) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_Riboseq_TPUTR.png",
    res = 300, height = 1000, width = 1750)
print(GC_3PUTR)
dev.off()

GC3_CDS <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = GC3_cont * 100, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = 0.25) +
  ylab("GC3%") +
  #gtitle("CDS") +
  publication_theme() +
  scale_y_continuous(limits = c(0,100)) +
  scale_fill_manual(values = my_colours) +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_Riboseq_GC3.png",
    res = 300, height = 1000, width = 1750)
print(GC3_CDS)
dev.off()

UTR5P_length <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = fputr_length, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = 0.25) +
  #ylab("Length (nts)") +
  publication_theme() +
  #ggtitle("5'UTR") +
  scale_y_continuous(trans = 'log10') +
  scale_fill_manual(values = my_colours) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_Riboseq_FPUTR_length.png",
    res = 300, height = 1000, width = 1750)
print(UTR5P_length)
dev.off()

CDS_length <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = cds_length, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = 0.25) +
  #ylab("Length (nts)") +
  publication_theme() +
  #ggtitle("CDS") +
  scale_y_continuous(trans = 'log10') +
  scale_fill_manual(values = my_colours) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_Riboseq_CDS_length.png",
    res = 300, height = 1000, width = 1750)
print(CDS_length)
dev.off()

TPUTR_length <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = tputr_length, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = 0.25) +
  #ylab("Length (nts)") +
  publication_theme() +
  #ggtitle("CDS") +
  scale_y_continuous(trans = 'log10') +
  scale_fill_manual(values = my_colours) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_Riboseq_TPUTR_length.png",
    res = 300, height = 1000, width = 1750)
print(TPUTR_length)
dev.off()


AG <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = AG_content * 100, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = 0.25) +
  ylab("AG%") +
  publication_theme() +
  #ggtitle("CDS") +
  scale_y_continuous(limits = c(0,100)) +
  scale_fill_manual(values = my_colours) +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_Riboseq_AGcontent.png",
    res = 300, height = 1000, width = 1750)
print(AG)
dev.off()


miRNA_sites <- 
  data_joined %>%
  mutate(miRNA_binding_sites = (miRNA_binding_sites + 1)) %>% 
  ggplot(aes(x = RPFs_group, y = miRNA_binding_sites, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = 0.25) +
  ylab("miRNA binding sites") +
  publication_theme() +
  ggtitle("miRNAs") +
  scale_y_continuous(trans = 'log10') +
  scale_fill_manual(values = my_colours) +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_Riboseq_microRNAs.png",
    res = 300, height = 1000, width = 1750)
print(miRNA_sites)
dev.off()

signal_peptides <- 
  inner_join(
    data_joined %>% 
      group_by(RPFs_group, signal_sequence) %>% 
      tally(),
    data_joined %>% 
      group_by(RPFs_group) %>% 
      tally() %>% 
      dplyr::rename(sum_n = n),
    by = c("RPFs_group")) %>% 
  mutate(perc = n/sum_n * 100) %>% 
  filter(signal_sequence == T) %>% 
  ggplot(aes(x = RPFs_group, y = perc)) +
  geom_bar(stat = "identity")

tmhmms <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/TM_domains_FINAL.csv") %>% 
  select(ENST, class) %>% unique()

secretory_info <- 
  data_joined %>% 
  left_join(tmhmms, by = c("transcript" = "ENST")) %>% 
  select(RPFs_group, transcript, class) %>% 
  mutate(class = case_when(is.na(class) ~ "cytosolic",.default = "ER")) 

tmhmm_data <- 
  inner_join(
    secretory_info %>% 
      group_by(RPFs_group, class) %>% 
      tally(),
    secretory_info %>% 
      group_by(RPFs_group) %>% 
      tally() %>% 
      dplyr::rename(sum_n = n),
    by = "RPFs_group") %>% 
  mutate(per = n/sum_n * 100) %>% 
  filter(class == "ER") %>% 
  ggplot(aes(x = RPFs_group, y = per, fill = RPFs_group, label = n)) +
  geom_bar(stat = "identity") +
  geom_text(vjust = -0.1, fontface = "bold") +
  scale_fill_manual(values = my_colours) +
  ylab("Percentage") +
  ggtitle("SP/TM presence") +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

TOPs <- read_csv("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/Gene_lists/TOPs/pnas.1912864117.sd01.csv",
                 col_names = F) %>% pull("X1")

# read in merged data
data_full <- read_csv(file.path(parent_dir, "Analysis/DESeq2_output/FINAL_siCNOT3_merged_DESeq2.csv"))

data_full %>% 
  inner_join(most_abundant_transcript, by = c("gene", "gene_sym")) %>% 
  left_join(secretory_info %>% select(-RPFs_group), by = "transcript") -> test
  
test$class <- factor(test$class, levels = c("cytosolic", "ER")) 

gplt <- 
ggplot() +
  geom_point(data = test %>% filter(class == "cytosolic"),
             aes(x = totals_log2FC, RPFs_log2FC), colour = "black", alpha = 0.5) +
  geom_point(data = test %>% filter(class == "ER"),
             aes(x = totals_log2FC, RPFs_log2FC), colour = "blue", alpha = 0.5)+
  xlab("Log2FoldChange mRNA") +
  ylab("Log2FoldChange RPFs") +
  publication_theme()

