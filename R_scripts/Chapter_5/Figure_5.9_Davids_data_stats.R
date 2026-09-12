data <- read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/davis_species_data.csv")

data_split <- split(data, data$species)


TukeyHSD(aov(value~ codon, data_split$`MDCK.II - SHH`))
TukeyHSD(aov(value~ codon, data_split$`MDCK.II - COL`))

TukeyHSD(aov(value~ codon, data_split$`3T3 - SHH`))
TukeyHSD(aov(value~ codon, data_split$`3T3 - COL`))
