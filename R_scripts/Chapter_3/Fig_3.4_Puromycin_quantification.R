library(tidyverse)

data <- 
read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/get_shit_done/chapter3/puromycin_quantification.csv")

publication_theme <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = 28, hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background = element_rect(colour = NA, fill = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.title.x = element_text(size = 16, vjust = -0.2),
            axis.text.y = element_text(size = 16), 
            axis.text.x = element_text(size = 16),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.text = element_text(size = 12),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.75, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(2,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

data_final <- 
data %>% 
  spread(key = band, value = adjusted_volume) %>% 
  mutate(normalised_PURO = PURO/VINC) %>% 
  mutate(replicate = as.character(replicate)) 

gplot_final <- 
data_final %>% 
  ggplot(aes(x = condition, y = normalised_PURO, colour = replicate)) +
  stat_summary(
    aes(group = condition),
    fun = mean,
    geom = "bar",
    fill = "grey70",
    colour = "black",
    alpha = 0.5,
    width = 0.6
  ) +
  geom_point(position = position_jitter(height = 0, width = 0.1)) +
  scale_y_continuous(limits = c(0,NA)) +
  scale_color_brewer(palette = "Set1") +
  ylab("PURO/VINC Ratio") +
  publication_theme() +
  theme(axis.title.x = element_blank())
  
png("\\\\data.beatson.gla.ac.uk/data/JETTLES/thesis_figures/Puro_quantification.png",
    res = 300, height = 1000, width = 1250)
print(gplot_final)
dev.off()

# stats
data_wide <- data_final %>%
  select(replicate, condition, normalised_PURO) %>% 
  pivot_wider(
    names_from = condition,
    values_from = normalised_PURO
  )

test_siCNOT3 <- t.test(
  data_wide$siCNOT3,
  data_wide$NTC,
  paired = T
)

test_sieS25 <- t.test(
  data_wide$sieS25,
  data_wide$NTC,
  paired = T
)

test_siCNOT3
test_sieS25

results <- tibble(
  comparison = c("siCNOT3 vs NTC", "sieS25 vs NTC"),
  p_value = c(
    test_siCNOT3$p.value,
    test_sieS25$p.value
  ),
  mean_difference = c(
    mean(data_wide$siCNOT3 - data_wide$NTC),
    mean(data_wide$sieS25 - data_wide$NTC)
  )
) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  )

results
