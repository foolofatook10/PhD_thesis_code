#load libraries
library(DESeq2)
library(tidyverse)
library(vsn)

#read in common variables
source("\\\\data.beatson.gla.ac.uk/data/R11/James/CNOT3_Sel_RiboSeq/FINAL/R_scripts/common_variables.R")

#create a variable for what the treatment is----
control <- "Tot"
treatment <- "CNOT3"

#read in the most abundant transcripts per gene csv file----
most_abundant_transcripts <- 
  read_csv(file = file.path("\\\\data.beatson.gla.ac.uk/data/R11/external_sequencing_data/Gillen_2021_RiboSeq/Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))

#read in data----
#the following for loop reads in each final CDS counts file and renames the counts column by the sample name and saves each data frame to a list
data_list <- list()
for (sample in RPF_sample_names) {
  df <- read_csv(file = file.path(parent_dir, "Analysis/CDS_counts", paste0(sample, "_pc_final_counts_all_frames.csv")), col_names = T)
  colnames(df) <- c("transcript", sample)
  data_list[[sample]] <- df
}

#merge all data within the above list using reduce
#Using full join retains all transcripts, but the NAs need to be replaced with 0 as these are transcipts that had 0 counts in that sample
#DESeq2 needs the transcripts to be as rownames, not as a column
data_list %>%
  purrr::reduce(full_join, by = "transcript") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames("transcript") -> RPF_counts

#create a data frame with the condition/replicate information----
#you need to make sure this data frame is correct for your samples, the below creates one for a n=3 with EFT226 treatment.
sample_info <- data.frame(row.names = RPF_sample_names_M,
                          condition = factor(rep(c(control,treatment),3)),
                          replicate = factor(c(rep(2:4, each = 2))))

#print the data frame to visually check it has been made as expected
sample_info

#make a DESeq data set from imported data----
DESeq2data <- DESeqDataSetFromMatrix(countData = RPF_counts,
                                     colData = sample_info,
                                     design = ~ replicate + condition)

#pre-filter to remove genes with less than an average of 10 counts across all samples----
keep <- rowMeans(counts(DESeq2data)) >= 10
table(keep)
DESeq2data <- DESeq2data[keep,]

#make sure levels are set appropriately so that Ctrl is "untreated"
DESeq2data$condition <- relevel(DESeq2data$condition, ref = control)

#run DESeq on DESeq data set----
dds <- DESeq(DESeq2data)

#extract results for each comparison----
res <- results(dds, contrast=c("condition", treatment, control))

#summarise results----
summary(res)

#write summary to a text file
sink(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("RPFs_", treatment, "_DEseq2_summary_M_28nts.txt")))
summary(res)
sink()

#apply LFC shrinkage for each comparison----
lfc_shrink <- lfcShrink(dds, coef=paste("condition", treatment, "vs", control, sep = "_"), type="apeglm")

#write reslts to csv----
as.data.frame(lfc_shrink[order(lfc_shrink$padj),]) %>%
  rownames_to_column("transcript") %>%
  inner_join(most_abundant_transcripts, by = "transcript") -> DEseq2_output
write_csv(DEseq2_output, file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("RPFs_", treatment, "_DEseq2_apeglm_LFC_shrinkage_M_28nts.csv")))

#extract normalised counts and plot SD vs mean----
ntd <- normTransform(dds) #this gives log2(n + 1)
vsd <- vst(dds, blind=FALSE) #Variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) #Regularized log transformation

meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

#write out normalised counts data----
#Regularized log transformation looks preferable for this data. Check for your own data and select the appropriate one
#The aim is for the range of standard deviations to be similar across the range of abundances, i.e. for the red line to be flat
as.data.frame(assay(vsd)) %>%
  rownames_to_column("transcript") %>%
  inner_join(most_abundant_transcripts, by = "transcript") -> normalised_counts
write_csv(normalised_counts, file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("RPFs_", treatment, "_normalised_counts_M_28nts.csv")))

#plot PCA----
pcaData <- plotPCA(vsd, intgroup=c("condition", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

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
            legend.text = element_text(size = 16),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", size = 16),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

png(filename = file.path(parent_dir, "plots/PCAs", paste0(treatment, "_RPFs_PCA_M_28nts.png")), width = 400, height = 350)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=6) +
  scale_shape_manual(values = c("2", "3", "4")) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
    guides(shape = F) +
  publication_theme() +
  ggtitle("Monosomes")
dev.off()

#apply batch correct and re-plot heatmap and PCA----
mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, vsd$replicate)
assay(vsd) <- mat

#PCA
pcaData <- plotPCA(vsd, intgroup=c("condition", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png(filename = file.path(parent_dir, "plots/PCAs", paste0(treatment, "_RPFs_batch_corrected_PCA_M_28nts.png")), width = 1500, height = 1500, res = 300)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=6) +
  scale_shape_manual(values = c("1", "2", "3")) +
  scale_color_manual(values = c("#FFC20A", "#0C7BDC")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  guides(shape = F) +
  publication_theme() +
  ggtitle("Monosomes")
dev.off()


