### signal peptide investigations ####

tmhmms <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/TM_domains_FINAL.csv") %>% 
  select(ENST, class) %>% unique()

# read in merged data
data_full <- read_csv(file.path(parent_dir, "Analysis/DESeq2_output/FINAL_sieS25_merged_DESeq2.csv"))

data_full <- data_full %>% filter(!RPFs_group %in% c("Totals up", "Totals down"))

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
                                                #"Totals down", "Totals up",
                                                "RPFs down", "RPFs up",
                                                "no change"),
                                     labels = c("Both\ndown", "Both\nup",
                                                #"Totals\ndown", "Totals\nup",
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

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/sieS25_Riboseq_SPTMs_barchart.png",
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

png(filename = "\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/eS25_Riboseq_SPTMs_scatterplot.png",
    res = 300, height = 1750, width = 1750)
print(gplt)
dev.off()