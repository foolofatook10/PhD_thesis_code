library(tidyverse)

CPMS <- 
read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/CDS_CPMs.csv")

mean_CPMs <- 
CPMS %>% 
  filter(MD30 == "M") %>% 
  mutate(CPM = (summed_CPM_codon/3)) %>% 
  select(replicate, IP, transcript, codon, CPM) %>% 
  group_by(IP, transcript, codon) %>% 
  summarise(CPM = mean(CPM))

master <- read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/master.csv")

RASs <- 
master %>% 
  filter(Gene %in% c("KRAS", "HRAS", "NRAS"))

CPMS_RASs <- 
inner_join(mean_CPMs,
           RASs %>% select(Gene, ENST),
           by = c("transcript" = "ENST"))

CPMS_RASs_split <- 
split(CPMS_RASs, CPMS_RASs$Gene)

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
            legend.title = element_text(face="italic", size = 14),
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

plot_RAS <- function(data_obj){
  
  ggplot(data = data_obj,
         aes(x = codon, y = CPM, colour = IP)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("#0C7BDC", "#FFC20A")) +
    ggtitle(paste0(unique(data_obj$Gene)))+
    geom_vline(xintercept = 142, lty = "dotted", colour = "black", size = 1) +
    publication_theme()
}


split_into_codons <- function(string) {
  
  # total length of string
  num.chars <- nchar(string)
  
  # the indices where each substr will start
  starts <- seq(1,num.chars, by=3)
  
  # chop it up
  sapply(starts, function(ii) {
    substr(string, ii, ii+2)
  })
}

RAS_ENST <- CPMS_RASs$transcript %>% unique()

KRAS_sequence <- 
master %>% 
  filter(ENST %in% RAS_ENST) %>% 
  filter(Gene == "KRAS") %>% 
  pull(nucleotide_sequence_CDS) %>% 
  split_into_codons()

HRAS_sequence <- 
  master %>% 
  filter(ENST %in% RAS_ENST) %>% 
  filter(Gene == "HRAS") %>% 
  pull(nucleotide_sequence_CDS) %>% 
  split_into_codons()

NRAS_sequence <- 
  master %>% 
  filter(ENST %in% RAS_ENST) %>% 
  filter(Gene == "NRAS") %>% 
  pull(nucleotide_sequence_CDS) %>% 
  split_into_codons()

paste(
HRAS_sequence[141],
HRAS_sequence[142],
HRAS_sequence[143])

paste(
KRAS_sequence[141],
KRAS_sequence[142],
KRAS_sequence[143])

paste(
NRAS_sequence[141],
NRAS_sequence[142],
NRAS_sequence[143])

gplots_RAS <- lapply(CPMS_RASs_split, plot_RAS)

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/NRAS.png",
    res = 300, height = 1000, width = 1250)
print(gplots_RAS$NRAS)
dev.off()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/KRAS.png",
    res = 300, height = 1000, width = 1250)
print(gplots_RAS$KRAS)
dev.off()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/HRAS.png",
    res = 300, height = 1000, width = 1250)
print(gplots_RAS$HRAS)
dev.off()



half_lives <- read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/CNOTdata_mRNAhalflives.csv")

features <- read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/features.csv")

deseq <- 
  read_csv(file.path("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/Analysis/DESeq2_output",
                     "RPFs_CNOT3_DEseq2_apeglm_LFC_shrinkage_M_28nts.csv"))

deseq_filt <- 
master %>% 
  filter(ENST %in% unique(deseq$transcript)) %>% 
  select(Gene, ENST, nucleotide_sequence_CDS)

tidy_up_data <- function(x){
  
  gene <- deseq_filt[x,]$Gene
  transcript <- deseq_filt[x,]$ENST
  sequence_split <- split_into_codons(deseq_filt[x,]$nucleotide_sequence_CDS)
  minus_atg <- sequence_split %>% head(-1) %>% tail(-1)
  
  df <- 
  data.frame(
    gene = gene,
    transcript = transcript,
             codon_id = minus_atg,
             codon_position = 2:(length(minus_atg) +1))
  
  return(df)
}

data_prepared <- 
do.call("rbind",
lapply(1:length(deseq_filt$nucleotide_sequence_CDS), 
       tidy_up_data))

CSCs <- 
read_csv("\\\\data.beatson.gla.ac.uk/data/R11/external_sequencing_data/Gillen_2021_halflife/Analysis/CSCs/CSCs.csv")%>% 
  mutate(codon = gsub("u", "t", tolower(codon)))

optimality_data <- 
data_prepared %>% 
  inner_join(CSCs %>% select(codon,
                             ctrl_cor_estimate, 
                             delta_cor_estimate), 
             by = c("codon_id" = "codon")) %>% 
  group_by(gene, transcript) %>% 
  summarise(ctrl_optimality = mean(ctrl_cor_estimate),
            siCNOT1_optimality = mean(delta_cor_estimate)) %>% 
  inner_join(deseq %>% select(transcript, log2FoldChange), 
             by = "transcript")

ggplot(optimality_data, aes(x = log2FoldChange, y = ctrl_optimality)) +
  geom_point()+
  ggpubr::stat_cor()

ggplot(optimality_data, aes(x = log2FoldChange, y = siCNOT1_optimality)) +
  geom_point() +
  ggpubr::stat_cor()



longest <- 
do.call("rbind",
lapply(split(features, features$Gene), 
       function(x){x %>%
           mutate(the_rank  = rank(-cds_length, ties.method = "first")) %>%
           filter(the_rank == 1) %>% 
           dplyr::select(-the_rank)}))

longest %>% 
  select(Gene, ENST, tct:tga) %>% 
  gather(codon, percent, tct:tga) %>% 
  inner_join(half_lives)
  group_by(Gene, ENST) %>% 
  summarise(sum_percent = sum(percent))

half_lives %>% 
  gather(key = condition, value = hl, siControl_half_life: log2FC_halflife) %>% 
  inner_join(longest %>% select(Gene, ENST, GC3_cont), 
                                                                                           by = c("Gene_ID" = "Gene")) %>% 
  ggplot(aes(x = GC3_cont, y = log10(hl))) +
  geom_point() +
  facet_wrap(~condition)



ip_enrichment <- 
deseq %>% 
  mutate(sig = case_when(padj <0.05 & log2FoldChange > 0 ~ "enriched",
                         padj < 0.05 & log2FoldChange < 0 ~ "depleted",
                         .default = "unchanged")) %>% 
  select(transcript, sig)

tidy_data <- 
ip_enrichment %>% 
  inner_join(master %>% select(Gene, ENST), 
             by = c("transcript" = "ENST")) %>% 
  inner_join(half_lives, by = c("Gene" = "Gene_ID"))

tidy_data %>% 
  gather(key = "condition", value = "half_life", siControl_half_life:log2FC_halflife) %>% 
  #filter(condition != "siCNOT1_half_life") %>% 
  ggplot(aes(x = sig, y = half_life)) +
  geom_violin() +
  geom_boxplot() +
  facet_wrap(~condition, scales = "free_y") 

tidy_data %>% 
  gather(key = "condition", value = "half_life", siControl_half_life:log2FC_halflife) %>% 
  #filter(condition == "siControl_half_life") %>% 
  ggplot(aes(x = sig, y = log10(half_life))) +
  geom_violin() +
  geom_boxplot() 

half_lives %>% 
  gather(key = condition, value = half_life, siControl_half_life:log2FC_halflife) %>% 
  inner_join(deseq, by = c("Gene_ID" = "gene_sym")) %>% 
  ggplot(aes(x = log2FoldChange, y = log10(half_life))) +
  geom_point() +
  facet_wrap(~condition) +
  ggpubr::stat_cor() +
  ylim(c(-3,4))
  
ctrl <- 
half_lives %>% 
  gather(key = condition, value = half_life, siControl_half_life:log2FC_halflife) %>% 
  inner_join(deseq, by = c("Gene_ID" = "gene_sym")) %>% 
  filter(condition == "siControl_half_life")

cor(ctrl$log2FoldChange, ctrl$half_life)

hl_genes <- half_lives$Gene_ID %>% unique()
dseq_gene <- deseq$gene_sym %>% unique()

intersect(hl_genes, dseq_gene)

half_lives %>% 
  gather(key = condition, value = half_life, siControl_half_life:log2FC_halflife) %>% 
  filter(condition == "log2FC_halflife") %>% 
  mutate(stabilisation = case_when(half_life < 0 ~"destabilised",
                                   half_life > 0 ~ "stabilised",
                                   .default = "unchanged")) %>% 
  inner_join(longest %>% select(Gene, GC3_cont),
             by = c("Gene_ID"  = "Gene")) %>% 
  ggplot(aes(x = stabilisation, y = GC3_cont)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  scale_y_continuous(limits = c(0, NA))

half_lives %>% 
  gather(key = condition, 
         value = half_life, 
         siControl_half_life:log2FC_halflife) %>% 
  filter(condition == "log2FC_halflife") %>% 
  # mutate(stabilisation = case_when(half_life < 0 ~"destabilised",
  #                                  half_life > 0 ~ "stabilised",
  #                                  .default = "unchanged")) %>% 
  inner_join(longest %>% select(Gene, GC3_cont),
             by = c("Gene_ID"  = "Gene")) %>% 
  ggplot(aes(x = GC3_cont, y = half_life)) +
  ggtitle("siCNOT1 Log2FoldChange") +
  geom_point() +
  ggpubr::stat_cor()

half_lives %>% 
  gather(key = condition, 
         value = half_life, 
         siControl_half_life:log2FC_halflife) %>% 
  filter(condition == "siControl_half_life") %>% 
  # mutate(stabilisation = case_when(half_life < 0 ~"destabilised",
  #                                  half_life > 0 ~ "stabilised",
  #                                  .default = "unchanged")) %>% 
  inner_join(longest %>% select(Gene, GC3_cont),
             by = c("Gene_ID"  = "Gene")) %>% 
  ggplot(aes(x = GC3_cont, y = log10(half_life))) +
  ggtitle("Ctrl conditions") +
  geom_point() +
  ylim(-0.5,2.5) +
  ggpubr::stat_cor()

# csc scores recalculated

codon_data_prepared <- 
longest %>% 
  filter(Gene %in% unique(half_lives$Gene_ID)) %>% 
  select(Gene, ENST, tct:tga) %>% 
  gather(key = "codon_id", value = "percentage", tct:tga) %>% 
  filter(!codon_id %in% c("atg", "taa", "tga", "tag")) %>% 
  inner_join(half_lives, by = c("Gene" = "Gene_ID"))

codon_data_prepared_split <- 
split(codon_data_prepared, codon_data_prepared$codon_id)

plot_codon_data <- function(data){
  
  final_gplot <- 
  data %>% 
    ggplot(aes(x = percentage, y = siControl_half_life))+
    geom_point() +
    ggtitle(paste0(unique(data$codon_id))) +
    ggpubr::stat_cor()
  
  return(final_gplot)
}

get_CSCs <- function(data){
  
  ctrl_cor <- cor(data$percentage, data$siControl_half_life)
  delta_cor <- cor(data$percentage, data$log2FC_halflife)
  
  df <- 
  data.frame(codon = unique(data$codon_id),
             ctrl_cor = ctrl_cor,
             delta_cor = delta_cor)
  
  return(df)
}

gplots <- lapply(codon_data_prepared_split, plot_codon_data)
new_CSCs <- do.call("rbind",lapply(codon_data_prepared_split, get_CSCs))
rownames(new_CSCs) <- NULL

new_CSCs <- 
new_CSCs %>% 
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A", "U", "G","C")))

new_CSCs %>% 
  ggplot(aes(x = reorder(codon, -delta_cor),
             y = delta_cor,
             fill = wobble)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

new_CSCs %>% 
  ggplot(aes(x = reorder(codon, -ctrl_cor),
             y = ctrl_cor,
             fill = wobble)) +
  geom_col()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# analysis 20/09/2025
gc3_cnot3selRP <- 
deseq %>% 
  filter(gene_sym %in% unique(half_lives$Gene_ID)) %>% 
  inner_join(features %>% select(ENST, GC3_cont),
             by = c("transcript" = "ENST")) %>% 
  mutate(class = case_when(padj < 0.05 & log2FoldChange > 0 ~"enriched",
                           padj < 0.05 & log2FoldChange < 0 ~ "depleted",
                           .default = "unchanged")) 

gc3_cnot3selRP$class <- factor(gc3_cnot3selRP$class,
                               levels = c("depleted",
                                          "unchanged",
                                          "enriched"))

gc3_cnot3selRP %>% 
  ggplot(aes(x = class, y = GC3_cont, fill = class)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(0, NA)) +
  publication_theme()
  

gc3_cnot3selRP %>% 
  group_by(class) %>% 
  tally() %>% 
  pull(n) %>% sum()

hl_data_c3 <- 
half_lives %>% 
  gather(key = condition, 
         value = half_life, 
         siControl_half_life:log2FC_halflife) %>% 
  filter(condition == "siControl_half_life") %>% 
  inner_join(deseq %>% filter(gene_sym %in% unique(half_lives$Gene_ID)),
             by = c("Gene_ID"  ="gene_sym")) %>% 
  mutate(enrichment = case_when(padj < 0.05 & log2FoldChange >0 ~"enriched",
                                padj < 0.05 & log2FoldChange <0 ~ "depleted",
                                .default = "unchanged")) 
  
hl_data_c3$enrichment <- factor(hl_data_c3$enrichment, 
                                levels = c("depleted",
                                           "unchanged",
                                           "enriched"),
                                labels = c("Dep.", "Unch.", "Enr.")) 
hl_data_c3 %>% 
ggplot(aes(x = enrichment, y = half_life, fill = enrichment)) +
  geom_violin() +
  scale_fill_brewer(palette = "Dark2") +
  geom_boxplot(width = 0.2)

ctrl_gplot <- 
hl_data_c3 %>% 
  ggplot(aes(x = enrichment, y = log10(half_life), fill = enrichment)) +
  geom_violin() +
  scale_fill_brewer(palette = "Dark2") +
  geom_boxplot(width = 0.2) +
  ylab(bquote(t[1/2])) +
  ggtitle("Ctrl") +
  publication_theme() +
  theme(legend.position = "none")

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/ctrl_half_lives_boxplots.png", height = 1250, width = 1250, res = 300)
print(ctrl_gplot)
dev.off()

hl_data_c3 <- 
  half_lives %>% 
  gather(key = condition, 
         value = half_life, 
         siControl_half_life:log2FC_halflife) %>% 
  filter(condition == "log2FC_halflife") %>% 
  inner_join(deseq %>% filter(gene_sym %in% unique(half_lives$Gene_ID)),
             by = c("Gene_ID"  ="gene_sym")) %>% 
  mutate(enrichment = case_when(padj < 0.05 & log2FoldChange >0 ~"enriched",
                                padj < 0.05 & log2FoldChange <0 ~ "depleted",
                                .default = "unchanged")) 


hl_data_c3$enrichment <- factor(hl_data_c3$enrichment, 
                                levels = c("depleted",
                                           "unchanged",
                                           "enriched"),
                                labels = c("Dep.", "Unch.", "Enr.")) 
siCNOT1 <- 
hl_data_c3 %>% 
  ggplot(aes(x = enrichment, y = half_life, fill = enrichment)) +
  geom_violin() +
  scale_fill_brewer(palette = "Dark2") +
  geom_boxplot(width = 0.2) +
  ylab(bquote(t[1/2])) +
  publication_theme() +
  ggtitle("siCNOT1") +
  theme(legend.position = "none",
        axis.title.x = element_blank())

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/siCNOT1_half_lives_boxplots.png", height = 1250, width = 1250, res = 300)
print(siCNOT1)
dev.off()
