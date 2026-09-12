library(tidyverse)
library(data.table)
library(DESeq2)
library(vsn)

print("loaded libraries")

home <- "\\\\data.beatson.gla.ac.uk"
#home <- "/home/local/BICR/jettles"

print("reading in CPMs")

# read in CPMs
CPMs <- fread(file = file.path(home, "data/R11/James/CNOT3_Sel_RiboSeq/FINAL/CPMs/CDS_CPMs.csv"), header = T)

# filter data for just the monosomes and tidy
counts <- 
  CPMs %>% 
  filter(MD30 == "M") %>% 
  dplyr::select(replicate, MD30, IP, transcript, codon, summed_counts_codon) %>% 
  dplyr::rename(counts = summed_counts_codon)

#filter out CNOT3 transcript
C3_transcripts <- read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/master.csv") %>% 
  filter(Gene == "CNOT3") %>% pull(ENST)

C3 <- counts %>% filter(transcript %in% C3_transcripts) %>% pull(transcript) %>% unique()

counts <- counts %>% filter(transcript != C3)

# make new column ID that fuses codon identifier to each transcript ID. i.e transcript_codon
counts$ID <- paste0(counts$transcript, "_", counts$codon)

# get rid of unecessary information
counts <- counts %>% dplyr::select(-transcript, -codon)

# create sample identifiers
counts$sample <- paste0(counts$replicate, "_", counts$IP)

# drop the redundant columns
counts <- counts %>% dplyr::select(- MD30, -replicate, -IP)

# spread the data for deseq 2
counts <- counts %>% spread(key = sample, value = counts)

# make the ID column the rownames
counts <- counts %>% column_to_rownames("ID")

print("finished pre processing")

#read in common variables
source(file.path(home, "data/R11/James/CNOT3_Sel_RiboSeq/FINAL/R_scripts/common_variables.R"))

# adjust sample names so that they're good and proper
RPF_sample_info_M$sample <- str_replace(RPF_sample_info_M$sample, "_M_", "_")

# convert the column "sample" to rownames
RPF_sample_info_M <- RPF_sample_info_M %>% column_to_rownames("sample")

# perform check to make sure the rownames of the metadata are in the same order as the column names of the counts data
all(rownames(RPF_sample_info_M) == colnames(counts))

# create indices to reorder
idx <- base::match(colnames(counts), rownames(RPF_sample_info_M))

# reorder the meta data
RPF_sample_info_M <- RPF_sample_info_M[idx,]

all(rownames(RPF_sample_info_M) == colnames(counts))

#create a variable for what the treatment is----
control <- "Tot"
treatment <- "CNOT3"

#RPF_sample_info_M <- RPF_sample_info_M %>% column_to_rownames("sample")

#read in the most abundant transcripts per gene csv file----
most_abundant_transcripts <- 
  read_csv(file = file.path(home,"data/R11/external_sequencing_data/Gillen_2021_RiboSeq/Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))

print("creating deseq2 object")

DESeq2data <- DESeqDataSetFromMatrix(counts,
                       colData = RPF_sample_info_M,
                       design = ~ replicate + condition)

#pre-filter to remove genes with less than an average of 10 counts across all samples----
keep <- rowMeans(counts(DESeq2data)) >= 10
table(keep)
DESeq2data <- DESeq2data[keep,]

#make sure levels are set appropriately so that Ctrl is "untreated"
DESeq2data$condition <- relevel(DESeq2data$condition, ref = control)

print("running deseq2")

#run DESeq on DESeq data set----
dds <- DESeq(DESeq2data)

#extract results for each comparison----
res <- results(dds, contrast=c("condition", treatment, control))

#summarise results----
summary(res)

#apply LFC shrinkage for each comparison----
lfc_shrink <- lfcShrink(dds, coef=paste("condition", treatment, "vs", control, sep = "_"), type="apeglm")

#read in the most abundant transcripts per gene csv file----
most_abundant_transcripts <- 
  read_csv(file = file.path("\\\\data.beatson.gla.ac.uk/data/R11/external_sequencing_data/Gillen_2021_RiboSeq/Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))

#write reslts to csv----
as.data.frame(lfc_shrink[order(lfc_shrink$padj),]) %>%
  rownames_to_column("transcript") %>%
  mutate(codon = str_remove(str_extract(string = transcript, pattern = "_\\d+"), "_"),
         transcript = str_remove(string = transcript, pattern = "_\\d+")) %>% 
  inner_join(most_abundant_transcripts, by = "transcript") -> DEseq2_output
write_csv(DEseq2_output, file = file.path(parent_dir, "Analysis/DESeq2_output", "M_RPFs_28nts.csv"))

#extract normalised counts and plot SD vs mean----
ntd <- normTransform(dds) #this gives log2(n + 1)
vsd <- vst(dds, blind=FALSE) #Variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) #Regularized log transformation

meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

as.data.frame(assay(rld)) %>%
  rownames_to_column("transcript") %>%
  mutate(codon = str_remove(str_extract(string = transcript, pattern = "_\\d+"), "_"),
         transcript = str_remove(string = transcript, pattern = "_\\d+")) %>% 
  inner_join(most_abundant_transcripts, by = "transcript") -> normalised_counts
write_csv(normalised_counts, file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("RPFs_", treatment, "_normalised_counts_M_28nts.csv")))

#plot PCA----
pcaData <- plotPCA(rld, intgroup=c("condition", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

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
            axis.text.x = element_text(size = 16, color = "black"),
            axis.text.y = element_text(size = 16, color = "black"),
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
mat <- assay(rld)
mat <- limma::removeBatchEffect(mat, rld$replicate)
assay(rld) <- mat

#PCA
pcaData <- plotPCA(rld, intgroup=c("condition", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData <- pcaData %>% mutate(replicate2 = case_when(replicate == "2" ~ "1",
                                                     replicate == "3" ~ "2",
                                                     replicate == "4" ~ "3"))

png(filename = file.path(parent_dir, "plots/PCAs", paste0(treatment, "_codons_deseq.png")), width = 1200, height = 1200, res = 300)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate2)) +
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
  ggtitle("DE Codons")
dev.off()


