library(tidyverse)

my_colours <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", "#D3D3D3")



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

most_abundant_transcript <- read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Riboseq/Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv")

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
         labels = c("Both down", "Both up",
                    "Totals down", "Totals up",
                    "RPFs down", "RPFs up",
                    "No change"))

summary(data_joined$RPFs_group)

RColorBrewer::brewer.pal.info

# get_stats("GC3_cont")
# get_stats("GC_cont_fp")
# get_stats("GC_cont_cds")
# get_stats("GC_cont_tp")
# get_stats("AG_content")
# get_stats("fputr_length")
# get_stats("cds_length")
# get_stats("tputr_length")
# get_stats("miRNA_binding_sites")


  
publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = 30, hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = 22),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(size = 16, color = "black", 
                                       angle = 90,
                                       hjust = 1,
                                       vjust = 0.5),
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


boxplot_width = 0.2

GC_5P <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = GC_cont_fp * 100, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = boxplot_width) +
  ylab("GC%") +
  ggtitle("5'UTR") +
  scale_y_continuous(limits = c(0,100)) +
  publication_theme() +
  scale_fill_manual(values = my_colours) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "black") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(margin = margin(b = 30)))

GC_CDS <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = GC_cont_cds * 100, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = boxplot_width) +
  ylab("GC%") +
  ggtitle("CDS") +
  scale_y_continuous(limits = c(0,100)) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "black") +
  publication_theme() +
  scale_fill_manual(values = my_colours) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(margin = margin(b = 30)))

GC_3P <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = GC_cont_tp * 100, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = boxplot_width) +
  ylab("GC%") +
  ggtitle("3'UTR") +
  scale_fill_manual(values = my_colours) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "black") +
  publication_theme() +
  scale_y_continuous(limits = c(0,100)) +
  scale_colour_manual(values = colours) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(margin = margin(b = 30)))


length_FPUTR <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = fputr_length, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  ylab("Length (nts)") +
  publication_theme() +
  ggtitle("5'UTR") +
  scale_y_continuous(trans = 'log10') +
  scale_fill_manual(values = my_colours) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "black") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(margin = margin(b = 30)))


length_CDS <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = cds_length, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = boxplot_width) +
  ylab("Length (nts)") +
  publication_theme() +
  ggtitle("CDS") +
  scale_y_continuous(trans = 'log10') +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "black") +
  scale_fill_manual(values = my_colours) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(margin = margin(b = 30)))


length_TPUTR <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = tputr_length, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  ylab("Length (nts)") +
  publication_theme() +
  ggtitle("3'UTR") +
  scale_y_continuous(trans = 'log10') +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "black") +
  scale_fill_manual(values = my_colours) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(margin = margin(b = 30)))

GC3_CDS <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = GC3_cont * 100, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = boxplot_width) +
  ylab("GC3%") +
  #gtitle("CDS") +
  publication_theme() +
  scale_y_continuous(limits = c(0,100)) +
  scale_fill_manual(values = my_colours) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "black") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(margin = margin(b = 30)))

AG <- 
  data_joined %>% 
  ggplot(aes(x = RPFs_group, y = AG_content * 100, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = boxplot_width) +
  ylab("AG%") +
  publication_theme() +
  #ggtitle("CDS") +
  scale_y_continuous(limits = c(0,100)) +
  scale_fill_manual(values = my_colours) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "black") +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        plot.title = element_text(margin = margin(b = 30)))

miRNA_sites <- 
  data_joined %>%
  mutate(miRNA_binding_sites = (miRNA_binding_sites + 1)) %>% 
  ggplot(aes(x = RPFs_group, y = miRNA_binding_sites, fill = RPFs_group)) +
  geom_violin() +
  geom_boxplot(width = boxplot_width) +
  #ylab("miRNA binding sites") +
  publication_theme() +
  #ggtitle("miRNAs") +
  scale_y_continuous(trans = 'log10') +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "black") +
  scale_fill_manual(values = my_colours) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(margin = margin(b = 30)))

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/siCNOT3_riboseq_features_all.png",
    res = 300, height = 4000, width = 3000)
cowplot::plot_grid(length_FPUTR, length_CDS, length_TPUTR,
                   GC_5P, GC_CDS, GC_3P, 
                   AG, GC3_CDS, miRNA_sites,
                   nrow = 3, align = "v")
dev.off()


master <- read_csv("\\\\data.beatson.gla.ac.uk/data/R11/James/sequences/master.csv")

######## RSCUs

seq_info_with_sequences <-
  inner_join(data_joined %>% 
               dplyr::rename(ENSG = gene,
                             ENST = transcript), 
             master %>% 
               select(ENSG, ENST, nucleotide_sequence_CDS), 
             by = c("ENSG","ENST")) %>% 
  select(RPFs_group, ENSG, ENST, nucleotide_sequence_CDS)

seq_info_with_sequences$nucleotide_sequence_CDS <- 
  strsplit(seq_info_with_sequences$nucleotide_sequence_CDS, split = "")



downregulated <- seq_info_with_sequences %>%
  filter(RPFs_group == "Both down")

upregulated <- seq_info_with_sequences %>%
  filter(RPFs_group == "Both up")

downregulated_names <- pull(downregulated, var = ENST)
downregulated_seqs <- pull(downregulated, var = nucleotide_sequence_CDS)
names(downregulated_seqs) <- downregulated_names

upregulated_names <- pull(upregulated, var = ENST)
upregulated_seqs <- pull(upregulated, var = nucleotide_sequence_CDS)
names(upregulated_seqs) <- upregulated_names
RSCU_lyst <- 
  list(lapply(downregulated_seqs, uco, frame = 0, index = c("eff", "freq", "rscu"), as.data.frame = T, NA.rscu = NA),
       lapply(upregulated_seqs,uco, frame = 0, index = c("eff", "freq", "rscu"), as.data.frame = T, NA.rscu = NA))

names(RSCU_lyst) <- c("downregulated", "upregulated")

data_upreg <- bind_rows(purrr::imap(RSCU_lyst[["upregulated"]], ~mutate(.x, transcript = .y))) %>%
  mutate(group = "upregulated")
data_downreg <- bind_rows(purrr::imap(RSCU_lyst[["downregulated"]], ~mutate(.x, transcript = .y))) %>%
  mutate(group = "downregulated")
row.names(data_upreg) <- NULL
row.names(data_downreg) <- NULL

data_combined <- rbind(data_downreg, data_upreg)


RCSUs <- data_combined %>%
  group_by(group, codon) %>%
  summarise(mean_RSCU = mean(RSCU, na.rm = T))

RCSUs <- RCSUs %>% 
  filter(codon != "taa" & codon != "tga" & codon != "tag") %>%
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A", "U", "G","C")))

RCSUs$group <- factor(RCSUs$group, levels = c("downregulated", "upregulated"),
                      labels = c("Both down", "Both up"))

RCSUs <- spread(RCSUs, group, mean_RSCU)

RCSUs <- 
  RCSUs %>% 
  mutate(TTG_label = case_when(codon == "ttg" ~ "ttg",
                               .default = "not_TTG"))

gplot_RSCU_groupboth <- 
  RCSUs %>% 
  ggplot() +
  geom_point(aes(x = `Both down`, y = `Both up`, colour = wobble), 
             shape = 19, alpha = 0.75, size = 3) +
  ggrepel::geom_text_repel(data = (RCSUs %>% filter(TTG_label == "ttg")), 
                           aes(x = `Both down`, y = `Both up`, label = TTG_label),
                           size = 6, colour = "#7570B3",
                           min.segment.length = 0,
                           segment.color = "black",
                           segment.size = 1,
                           nudge_y = +0.25,
                           nudge_x = -0.5
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), name = "Wobble Base Identity" )+
  geom_abline(slope = 1, intercept = 0, linewidth = 1, linetype = "dashed") +
  #ggtitle("Monosomes") +
  #xlim(0,3)+
  #ylim(0,2.5)+
  xlab("Both down") + 
  ylab("Both up") +
  ggtitle("Mean RSCU") +
  xlim(0,2.5) +
  ylim(0, 2.5) +
  publication_theme() +
  theme(legend.position = "none") +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(siz = 22),
    axis.title.x = element_text(siz = 22))

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/siCNOT3_riboseq_features_RSCU_bothup_bothdown.png",
    res = 300, height = 1250, width = 1250)
print(gplot_RSCU_groupboth)
dev.off()

#########
seq_info_with_sequences <-
  inner_join(data_joined %>% 
               dplyr::rename(ENSG = gene,
                             ENST = transcript), 
             master %>% 
               select(ENSG, ENST, nucleotide_sequence_CDS), 
             by = c("ENSG","ENST")) %>% 
  select(RPFs_group, ENSG, ENST, nucleotide_sequence_CDS)

seq_info_with_sequences$nucleotide_sequence_CDS <- 
  strsplit(seq_info_with_sequences$nucleotide_sequence_CDS, split = "")



downregulated <- seq_info_with_sequences %>%
  filter(RPFs_group == "RPFs down")

upregulated <- seq_info_with_sequences %>%
  filter(RPFs_group == "RPFs up")

downregulated_names <- pull(downregulated, var = ENST)
downregulated_seqs <- pull(downregulated, var = nucleotide_sequence_CDS)
names(downregulated_seqs) <- downregulated_names

upregulated_names <- pull(upregulated, var = ENST)
upregulated_seqs <- pull(upregulated, var = nucleotide_sequence_CDS)
names(upregulated_seqs) <- upregulated_names
RSCU_lyst <- 
  list(lapply(downregulated_seqs, uco, frame = 0, index = c("eff", "freq", "rscu"), as.data.frame = T, NA.rscu = NA),
       lapply(upregulated_seqs,uco, frame = 0, index = c("eff", "freq", "rscu"), as.data.frame = T, NA.rscu = NA))

names(RSCU_lyst) <- c("downregulated", "upregulated")

data_upreg <- bind_rows(purrr::imap(RSCU_lyst[["upregulated"]], ~mutate(.x, transcript = .y))) %>%
  mutate(group = "upregulated")
data_downreg <- bind_rows(purrr::imap(RSCU_lyst[["downregulated"]], ~mutate(.x, transcript = .y))) %>%
  mutate(group = "downregulated")
row.names(data_upreg) <- NULL
row.names(data_downreg) <- NULL

data_combined <- rbind(data_downreg, data_upreg)


RCSUs <- data_combined %>%
  group_by(group, codon) %>%
  summarise(mean_RSCU = mean(RSCU, na.rm = T))

RCSUs <- RCSUs %>% 
  filter(codon != "taa" & codon != "tga" & codon != "tag") %>%
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A", "U", "G","C")))

RCSUs$group <- factor(RCSUs$group, levels = c("downregulated", "upregulated"),
                      labels = c("RPFs down", "RPFs up"))

RCSUs <- spread(RCSUs, group, mean_RSCU)

RCSUs <- 
  RCSUs %>% 
  mutate(TTG_label = case_when(codon == "ttg" ~ "ttg",
                               .default = "not_TTG"))

gplot_RSCU_groupboth <- 
  RCSUs %>% 
  ggplot() +
  geom_point(aes(x = `RPFs down`, y = `RPFs up`, colour = wobble), 
             shape = 19, alpha = 0.75, size = 3) +
  ggrepel::geom_text_repel(data = (RCSUs %>% filter(TTG_label == "ttg")), 
                           aes(x = `RPFs down`, y = `RPFs up`, label = TTG_label),
                           size = 6, colour = "#7570B3",
                           min.segment.length = 0,
                           segment.color = "black",
                           segment.size = 1,
                           nudge_y = +0.25,
                           nudge_x = -0.5
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), name = "Wobble Base Identity" )+
  geom_abline(slope = 1, intercept = 0, linewidth = 1, linetype = "dashed") +
  #ggtitle("Monosomes") +
  #xlim(0,3)+
  #ylim(0,2.5)+
  xlab("RPFs down") + 
  ylab("RPFs up") +
  ggtitle("Mean RSCU") +
  xlim(0,2.5) +
  ylim(0, 2.5) +
  publication_theme() +
  theme(legend.position = "none") +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(siz = 22),
    axis.title.x = element_text(siz = 22))

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/siCNOT3_riboseq_features_RSCU_RPFsupdown.png",
    res = 300, height = 1250, width = 1250)
print(gplot_RSCU_groupboth)
dev.off()

#########
seq_info_with_sequences <-
  inner_join(data_joined %>% 
               dplyr::rename(ENSG = gene,
                             ENST = transcript), 
             master %>% 
               select(ENSG, ENST, nucleotide_sequence_CDS), 
             by = c("ENSG","ENST")) %>% 
  select(RPFs_group, ENSG, ENST, nucleotide_sequence_CDS)

seq_info_with_sequences$nucleotide_sequence_CDS <- 
  strsplit(seq_info_with_sequences$nucleotide_sequence_CDS, split = "")



downregulated <- seq_info_with_sequences %>%
  filter(RPFs_group == "Totals down")

upregulated <- seq_info_with_sequences %>%
  filter(RPFs_group == "Totals up")

downregulated_names <- pull(downregulated, var = ENST)
downregulated_seqs <- pull(downregulated, var = nucleotide_sequence_CDS)
names(downregulated_seqs) <- downregulated_names

upregulated_names <- pull(upregulated, var = ENST)
upregulated_seqs <- pull(upregulated, var = nucleotide_sequence_CDS)
names(upregulated_seqs) <- upregulated_names
RSCU_lyst <- 
  list(lapply(downregulated_seqs, uco, frame = 0, index = c("eff", "freq", "rscu"), as.data.frame = T, NA.rscu = NA),
       lapply(upregulated_seqs,uco, frame = 0, index = c("eff", "freq", "rscu"), as.data.frame = T, NA.rscu = NA))

names(RSCU_lyst) <- c("downregulated", "upregulated")

data_upreg <- bind_rows(purrr::imap(RSCU_lyst[["upregulated"]], ~mutate(.x, transcript = .y))) %>%
  mutate(group = "upregulated")
data_downreg <- bind_rows(purrr::imap(RSCU_lyst[["downregulated"]], ~mutate(.x, transcript = .y))) %>%
  mutate(group = "downregulated")
row.names(data_upreg) <- NULL
row.names(data_downreg) <- NULL

data_combined <- rbind(data_downreg, data_upreg)


RCSUs <- data_combined %>%
  group_by(group, codon) %>%
  summarise(mean_RSCU = mean(RSCU, na.rm = T))

RCSUs <- RCSUs %>% 
  filter(codon != "taa" & codon != "tga" & codon != "tag") %>%
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A", "U", "G","C")))

RCSUs$group <- factor(RCSUs$group, levels = c("downregulated", "upregulated"),
                      labels = c("Totals down", "Totals up"))

RCSUs <- spread(RCSUs, group, mean_RSCU)

RCSUs <- 
  RCSUs %>% 
  mutate(TTG_label = case_when(codon == "ttg" ~ "ttg",
                               .default = "not_TTG"))

gplot_RSCU_groupboth <- 
  RCSUs %>% 
  ggplot() +
  geom_point(aes(x = `Totals down`, y = `Totals up`, colour = wobble), 
             shape = 19, alpha = 0.75, size = 3) +
  ggrepel::geom_text_repel(data = (RCSUs %>% filter(TTG_label == "ttg")), 
                           aes(x = `Totals down`, y = `Totals up`, label = TTG_label),
                           size = 6, colour = "#7570B3",
                           min.segment.length = 0,
                           segment.color = "black",
                           segment.size = 1,
                           nudge_y = +0.25,
                           nudge_x = -0.5
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), name = "Wobble Base Identity" )+
  geom_abline(slope = 1, intercept = 0, linewidth = 1, linetype = "dashed") +
  #ggtitle("Monosomes") +
  #xlim(0,3)+
  #ylim(0,2.5)+
  xlab("Totals down") + 
  ylab("Totals up") +
  ggtitle("Mean RSCU") +
  xlim(0,2.5) +
  ylim(0, 2.5) +
  publication_theme() +
  theme(legend.position = "none") +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(siz = 22),
    axis.title.x = element_text(siz = 22))

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/siCNOT3_riboseq_features_RSCU_Totalsupdown.png",
    res = 300, height = 1250, width = 1250)
print(gplot_RSCU_groupboth)
dev.off()


#########




tmhmms <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/TM_domains_FINAL.csv") %>% 
  select(ENST, class) %>% unique()

# read in merged data
data_full <- read_csv(file.path(parent_dir, "Analysis/DESeq2_output/FINAL_siCNOT3_merged_DESeq2.csv"))

secretory_info <- 
  data_full %>% 
  inner_join(most_abundant_transcript, by = c("gene", "gene_sym")) %>% 
  left_join(tmhmms, by = c("transcript" = "ENST")) %>% 
  select(RPFs_group, transcript, class) %>% 
  mutate(class = case_when(is.na(class) ~ "cytosolic",
                           .default = "ER")) 

secretory_info2 <- secretory_info %>% filter(RPFs_group != "NS")

secretory_info2$RPFs_group <- factor(secretory_info2$RPFs_group, 
                                    levels = c("both down", "both up",
                                               "Totals down", "Totals up",
                                               "RPFs down", "RPFs up",
                                               "no change"),
                                    labels = c("Both\ndown", "Both\nup",
                                               "Totals\ndown", "Totals\nup",
                                               "RPFs\ndown", "RPFs\nup",
                                               "No\nchange"))

publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = 24, hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = 22),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(size = 16, color = "black"),
            axis.text.y = element_text(size = 16, color = "black"),
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

tmhmm_data <- 
  inner_join(
    secretory_info2 %>% 
      group_by(RPFs_group, class) %>% 
      tally(),
    secretory_info2 %>% 
      group_by(RPFs_group) %>% 
      tally() %>% 
      dplyr::rename(sum_n = n),
    by = "RPFs_group") %>% 
  mutate(per = n/sum_n * 100) %>% 
  filter(class == "ER") %>% 
  filter(RPFs_group != "NS") %>% 
  ggplot(aes(x = RPFs_group, y = per, fill = RPFs_group, label = n)) +
  geom_bar(stat = "identity") +
  geom_text(vjust = -0.1, fontface = "bold") +
  scale_fill_manual(values = my_colours) +
  ylab("Proportion of transcripts (%)") +
  ggtitle("SP/TM presence") +
  publication_theme() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_Riboseq_SPTMs_barchart.png",
    res = 300, height = 1500, width = 1750)
print(tmhmm_data)
dev.off()


secretory_info$class <- factor(secretory_info$class, levels = c("cytosolic", "ER"))

#test$class <- factor(test$class, levels = c("cytosolic", "ER")) 

publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = 24, hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = 22),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(size = 18, color = "black"),
            axis.text.y = element_text(size = 18, color = "black"),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = 18),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

final_secretory_info <- 
data_full %>% 
  inner_join(most_abundant_transcript, by = c("gene", "gene_sym")) %>% 
  inner_join(secretory_info, by = c("transcript", "RPFs_group"))

gplt <- 
  ggplot() +
  geom_point(data = final_secretory_info %>% filter(class == "cytosolic"),
             aes(x = totals_log2FC, RPFs_log2FC), colour = "black", alpha = 0.5) +
  geom_point(data = final_secretory_info %>% filter(class == "ER"),
             aes(x = totals_log2FC, RPFs_log2FC), colour = "blue", alpha = 0.5)+
  xlab("Log2FoldChange mRNA") +
  ylab("Log2FoldChange RPFs") +
  ggtitle("SP/TM presence") +
  publication_theme()

gplt_lgd <- 
ggplot() +
  geom_point(data = final_secretory_info %>% 
               mutate(class = case_when(class == "cytosolic" ~ "Cytosolic",
                                        class == "ER" ~ "Contains SP/TM")),
             aes(x = totals_log2FC, RPFs_log2FC, color = class), size = 2) +
  scale_color_manual(values = c(alpha("blue",0.5), alpha("black",0.5))) +
  xlab("Log2FoldChange mRNA") +
  ylab("Log2FoldChange RPFs") +
  ggtitle("SP/TM presence") +
  publication_theme() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank())

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_Riboseq_SPTMs_scatterplot_legend.png",
    res = 300, height = 200, width = 750)
cowplot::ggdraw(cowplot::get_legend(gplt_lgd))
dev.off()

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CNOT3_Riboseq_SPTMs_scatterplot.png",
    res = 300, height = 1750, width = 1750)
print(gplt)
dev.off()

final_secretory_info %>% filter(RPFs_group == "both down")

df <- data.frame(
  group = c("both_down", "both_up", "totals_down", "totals_up",
            "RPFs_down", "RPFs_up", "no_change", "NS"),
  contain = c(706, 149, 28, 70, 59, 7, 203, 859),
  total = c(1668, 1587, 368, 331, 118, 239, 1114, 4034)
)

results <- map_dfr(seq_len(nrow(df)), function(i) {
  
  # focal group
  a <- df$contain[i]
  b <- df$total[i] - a
  
  # all other groups combined
  c <- sum(df$contain[-i])
  d <- sum(df$total[-i]) - c
  
  test <- fisher.test(
    matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  )
  
  data.frame(
    group = df$group[i],
    proportion = a / (a + b),
    background_proportion = c / (c + d),
    odds_ratio = unname(test$estimate),
    p_value = test$p.value
  )
}) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH")
  )
