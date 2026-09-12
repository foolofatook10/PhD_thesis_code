library(tidyverse)
library(data.table)

axis_label_y = "CSC score"

dts <- read_csv("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/dwell_times_CNOT3SelRP_entire_transcriptome.csv") %>% 
  select(-`...1`)

dts_summarised <- 
  dts %>% 
  group_by(sample_id, codon) %>% 
  summarise(`-4` = sum(minus4),
            `-3` = sum(minus3),
            `E` = sum(E_site),
            `P` = sum(P_site),
            `A` = sum(A_site),
            `+1` = sum(plus1),
            `+2` = sum(plus2),
            total = sum(total))

dts_summarised_split <- split(dts_summarised, dts_summarised$sample_id)

downstream_processing <- function(df){
  
  df <-
    df %>%
    filter(codon != "tga" & codon != "tag" & codon != "taa") %>%
    mutate(minus_4 = (log2(`-4`)) - log2(total / 7),
           minus_3 = (log2(`-3`)) - log2(total / 7),
           E_site = (log2(E)) - log2(total / 7),
           P_site = (log2(P)) - log2(total / 7),
           A_site = (log2(A)) - log2(total / 7),
           plus_1 = (log2(`+1`)) - log2(total / 7),
           plus_2 = (log2(`+2`)) - log2(total / 7)) %>%
    select(codon, minus_4, minus_3, E_site, P_site, A_site, plus_1, plus_2, sample_id) %>%
    dplyr::rename(sample = sample_id)
  
  return(df)
}



data <- 
  do.call("rbind", lapply(dts_summarised_split, downstream_processing)) %>% 
  ungroup()



CSCs <- 
read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/CSCs.csv") %>% 
  mutate(codon = tolower(gsub("U", "T", codon))) %>% 
  select(codon, 
         CSC_ctrl = ctrl_cor_estimate,
         CSC_siCNOT1 = delta_cor_estimate)

# publication theme
publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(
                                      size = 22, hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = 20),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(size = 18, color = "black"),
            axis.text.y = element_text(size = 18, color = "black"),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = 16),
            legend.key = element_rect(colour = NA),
            legend.position = "none",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(5,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


A_site_dts <- 
  data %>% 
  filter(sample %in% c("REP2_M_Tot", "REP3_M_Tot", "REP4_M_Tot")) %>% 
  gather(key = ribosome_position, value = dwell_time, minus_4:plus_2) %>% 
  group_by(codon, ribosome_position) %>% 
  summarise(mean_dt = mean(dwell_time)) %>% 
  filter(ribosome_position == "A_site") %>% 
  arrange(-mean_dt)

Totals_A_site <- 
  inner_join(A_site_dts, CSCs, by = "codon") %>% 
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A", "U", "G","C")))
  #filter(wobble %in% c("A", "U")) %>%  

Totals_cor_gplot <- 
Totals_A_site %>%   
ggplot(aes(x = mean_dt, y = CSC_ctrl, colour = wobble, label = codon)) +
  geom_point(size = 2) +
  #geom_text() +
  # ggpubr::stat_cor(method="pearson", 
  #                  inherit.aes = F, 
  #                  data = Totals_A_site, 
  #                  aes(x = mean_dt, y = CSC_ctrl)) +
  scale_color_brewer(palette = "Set1") +
  ggtitle("All codons") +
  xlab("A site dwell time") +
  ylab(axis_label_y) +
  publication_theme()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CSC_A_site_Totals_all.png",
    res = 300, height = 1500, width = 1500)
print(Totals_cor_gplot)
dev.off()

RColorBrewer::brewer.pal(4, "Set1")

Totals_GC <-  
Totals_A_site %>% 
  filter(wobble %in% c("G", "C")) %>%
  filter(codon != "ttg") 
  # filter(!codon %in% c("ttg", "cgg", "cgc", "agc", "tgc", "ccc",
  #                      "atc", "gtc", "aac", "aag", "gac", "tac"))  

Totals_all_GC_codons_gplot <- 
Totals_GC %>% 
  #filter(wobble == "C") %>% 
  ggplot(aes(x = mean_dt, y = CSC_ctrl,colour = wobble, label = codon)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(show.legend = F, size = 5) +
  scale_color_manual(values = c("#4DAF4A", "#984EA3")) +
  # ggpubr::stat_cor(method="pearson", 
  #                  inherit.aes = F, 
  #                  data = Totals_GC, 
  #                  aes(x = mean_dt, y = CSC_ctrl),
  #                  label.x.npc = "left", label.y.npc = "top",
  #                  size = 3) +
  #scale_color_brewer(palette = "Set1") +
  ggtitle("GC3 codons") +
  xlab("A site dwell time") +
  ylab(axis_label_y) +
  publication_theme()

# png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CSC_A_site_Totals_GC_codons.png",
#     res = 300, height = 1500, width = 1500)
# print(Totals_all_GC_codons_gplot)
# dev.off()

# select_GC_codons <- 
# Totals_A_site %>% 
#   filter(wobble %in% c("G", "C")) %>%
#   filter(codon != "ttg") %>% 
#  filter(!codon %in% c("ttg", "cgg", "cgc", "agc", "tgc", "ccc",
#                       "atc", "gtc", "aac", "aag", "gac", "tac")) 
# 
# Totals_select_GC_codons_gplot <- 
# select_GC_codons %>% 
#   ggplot(aes(x = mean_dt, y = CSC_ctrl,colour = wobble, label = codon)) +
#   geom_point() +
#   geom_text(show.legend = F) +
#   scale_color_manual(values = c("#4DAF4A", "#984EA3")) +
#   ggpubr::stat_cor(method="pearson", 
#                    inherit.aes = F, 
#                    data = select_GC_codons, 
#                    aes(x = mean_dt, y = CSC_ctrl)) +
#   #scale_color_brewer(palette = "Set1") +
#   ggtitle("Totals") +
#   xlab("A site dwell time") +
#   ylab("Ctrl CSC") +
#   publication_theme()
# 
# png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CSC_A_site_Totals_selectGC_codons.png",
#     res = 300, height = 1500, width = 1500)
# print(Totals_select_GC_codons_gplot)
# dev.off()

Totals_AU <- 
  Totals_A_site %>% 
  filter(wobble %in% c("A", "U") | codon == "ttg")

Totals_all_AU_codons_gplot <- 
Totals_AU %>% 
  # filter(!codon %in% c("ttg", "cgg", "cgc", "agc", "tgc", "ccc",
  #                      "atc", "gtc", "aac", "aag", "gac", "tac")) %>%  
  ggplot(aes(x = mean_dt, y = CSC_ctrl,colour = wobble, label = codon)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(show.legend = F, size = 5) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  # ggpubr::stat_cor(method="pearson",
  #                  inherit.aes = F,
  #                  data = Totals_AU,
  #                  aes(x = mean_dt, y = CSC_ctrl)) +
  scale_color_brewer(palette = "Set1") +
  ggtitle("AU3 Codons") +
  xlab("A site dwell time") +
  ylab(axis_label_y) +
  publication_theme()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CSC_A_site_Totals_AU_codons.png",
    res = 300, height = 1500, width = 1500)
print(Totals_all_AU_codons_gplot)
dev.off()


A_site_dts <- 
  data %>% 
  filter(sample %in% c("REP2_M_CNOT3", "REP3_M_CNOT3", "REP4_M_CNOT3")) %>% 
  gather(key = ribosome_position, value = dwell_time, minus_4:plus_2) %>% 
  group_by(codon, ribosome_position) %>% 
  summarise(mean_dt = mean(dwell_time)) %>% 
  filter(ribosome_position == "A_site") %>% 
  arrange(-mean_dt)





CNOT3_A_site <- 
  inner_join(A_site_dts, CSCs, by = "codon") %>% 
  mutate(wobble = 
           factor(str_sub(codon, 3,3), 
                  levels = c("a", "t", "g", "c"), 
                  labels = c("A", "U", "G","C")))
#filter(wobble %in% c("A", "U")) %>%  
#### here ####
CNOT3_cor_gplot <- 
  CNOT3_A_site %>%   
  ggplot(aes(x = mean_dt, y = CSC_ctrl, colour = wobble, label = codon)) +
  geom_point(size = 2) +
  #geom_text() +
  # ggpubr::stat_cor(method="pearson", 
  #                  inherit.aes = F, 
  #                  data = CNOT3_A_site, 
  #                  aes(x = mean_dt, y = CSC_ctrl)) +
  scale_color_brewer(palette = "Set1") +
  ggtitle("All Codons") +
  xlab("A site dwell time") +
  ylab(axis_label_y) +
  publication_theme()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CSC_A_site_CNOT3_all.png",
    res = 300, height = 1500, width = 1500)
print(CNOT3_cor_gplot)
dev.off()

RColorBrewer::brewer.pal(4, "Set1")

CNOT3_GC <-  
  CNOT3_A_site %>% 
  filter(wobble %in% c("G", "C")) %>%
  filter(codon != "ttg") 
# filter(!codon %in% c("ttg", "cgg", "cgc", "agc", "tgc", "ccc",
#                      "atc", "gtc", "aac", "aag", "gac", "tac"))  

CNOT3_all_GC_codons_gplot <- 
  CNOT3_GC %>% 
  ggplot(aes(x = mean_dt, y = CSC_ctrl,colour = wobble, label = codon)) +
  geom_point() +
  ggrepel::geom_text_repel(show.legend = F, size = 5) +
  scale_color_manual(values = c("#4DAF4A", "#984EA3")) +
  # ggpubr::stat_cor(method="pearson", 
  #                  inherit.aes = F, 
  #                  data = CNOT3_GC, 
  #                  aes(x = mean_dt, y = CSC_ctrl)) +
  #scale_color_brewer(palette = "Set1") +
  ggtitle("GC3 Codons") +
  xlab("A site dwell time") +
  ylab(axis_label_y) +
  publication_theme()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CSC_A_site_CNOT3_GC_codons.png",
    res = 300, height = 1500, width = 1500)
print(CNOT3_all_GC_codons_gplot)
dev.off()

# select_GC_codons <- 
#   CNOT3_A_site %>% 
#   filter(wobble %in% c("G", "C")) %>%
#   filter(codon != "ttg") %>% 
#   filter(!codon %in% c("ttg", "cgg", "cgc", "agc", "tgc", "ccc",
#                        "atc", "gtc", "aac", "aag", "gac", "tac")) 
# 
# CNOT3_select_GC_codons_gplot <- 
#   select_GC_codons %>% 
#   ggplot(aes(x = mean_dt, y = CSC_ctrl,colour = wobble, label = codon)) +
#   geom_point() +
#   geom_text(show.legend = F) +
#   scale_color_manual(values = c("#4DAF4A", "#984EA3")) +
#   ggpubr::stat_cor(method="pearson", 
#                    inherit.aes = F, 
#                    data = select_GC_codons, 
#                    aes(x = mean_dt, y = CSC_ctrl)) +
#   #scale_color_brewer(palette = "Set1") +
#   ggtitle("CNOT3") +
#   xlab("A site dwell time") +
#   ylab("Ctrl CSC") +
#   publication_theme()
# 
# png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CSC_A_site_CNOT3_selectGC_codons.png",
#     res = 300, height = 1500, width = 1500)
# print(CNOT3_select_GC_codons_gplot)
# dev.off()

CNOT3_AU <- 
  CNOT3_A_site %>% 
  filter(wobble %in% c("A", "U") | codon == "ttg")

CNOT3_all_AU_codons_gplot <- 
  CNOT3_AU %>% 
  # filter(!codon %in% c("ttg", "cgg", "cgc", "agc", "tgc", "ccc",
  #                      "atc", "gtc", "aac", "aag", "gac", "tac")) %>%  
  ggplot(aes(x = mean_dt, y = CSC_ctrl,colour = wobble, label = codon)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(show.legend = F, size = 5) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  # ggpubr::stat_cor(method="pearson", 
  #                  inherit.aes = F, 
  #                  data = CNOT3_AU, 
  #                  aes(x = mean_dt, y = CSC_ctrl)) +
  #scale_color_brewer(palette = "Set1") +
  ggtitle("AU3 Codons") +
  xlab("A site dwell time") +
  ylab(axis_label_y) +
  publication_theme()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CSC_A_site_CNOT3_AU_codons.png",
    res = 300, height = 1500, width = 1500)
print(CNOT3_all_AU_codons_gplot)
dev.off()

png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/CSC_versus_A_site_dwell_time.png",
    res = 300, height = 3500, width = 3000)
cowplot::plot_grid(Totals_cor_gplot, CNOT3_cor_gplot,
                   Totals_all_AU_codons_gplot, CNOT3_all_AU_codons_gplot,
                   Totals_all_GC_codons_gplot, CNOT3_all_GC_codons_gplot, 
                   ncol = 2, nrow = 3, 
                   align = "vh")
dev.off()
