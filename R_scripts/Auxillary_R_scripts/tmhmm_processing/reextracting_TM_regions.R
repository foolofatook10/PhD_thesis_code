library(tidyverse)
library(data.table)

# set working directory
#setwd("\\\\data.beatson.gla.ac.uk/data/JETTLES/CTG_codons/bioinformatic_analysis/plots")

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

# read in THMHMM data
data <- 
  read_csv("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/merged_deepTMHMM_filtered.csv") %>% 
  dplyr::rename(codon_start = start,
                codon_end= end) %>% 
  dplyr::rename(ENSP = protein)


# read in the Ensembl identifiers
conversion_table <- 
  fread("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/gencode.v38.pc_transcripts_protein_IDs.csv", header = F)
# name the columns
colnames(conversion_table) <- c("ENSP", "ENST", "ENSG", "Gene")

# read in features table
features <- as_tibble(fread("C:/Users/Jettles/OneDrive - University of Glasgow/Documents/features.csv", drop = "V1")) 

# read in master table
master <- as_tibble(fread(file = "C:/Users/Jettles/OneDrive - University of Glasgow/Documents/master.csv", header = T, drop = "V1")) %>% 
  select(Gene, ENSG, ENST, nucleotide_sequence_CDS)

# select the longest transcript per gene
# longest <- master %>% 
#   mutate(polypeptide_length = nchar(master$amino_acid_sequence)) %>%
#   group_by(ENSG) %>%
#   mutate(the_rank  = rank(-polypeptide_length, ties.method = "first")) %>%
#   filter(the_rank == 1) %>% 
#   select(-the_rank) %>% 
#   select(Gene, ENSG, ENST, nucleotide_sequence_CDS)

# join together all the data
data <- 
  inner_join(data, conversion_table, by = "ENSP") %>% 
  select(Gene, ENSG, ENST, ENSP, location, codon_start, codon_end) %>% 
  inner_join(master, by = c("Gene", "ENSG", "ENST")) %>% 
  select(-ENSG, -ENSP)

# first step is to idenify classes
data$location <- factor(data$location)
summary(data)

# some transcripts contain the ID tag "Beta sheet" and "periplasm"
weird_transcripts <- 
  data %>% 
  filter(location %in% c("Beta sheet", "periplasm")) %>% 
  pull(var = "ENST") %>% 
  unique()

#
length(data$ENST %>% unique())

# filter out the weird transcripts
data <- data %>% filter(!ENST %in% weird_transcripts)

# check numbers of unique transcripts
length(weird_transcripts)
length(unique(data$ENST))

# drop the unused levels
data$location <- droplevels(data$location)

# recall summary
summary(data) # See it. Say it. Sorted.


# define the classes of transcript
classes <- 
  data %>%
  filter(location %in% c("signal", "TMhelix")) %>% 
  group_by(ENST, location ) %>% 
  tally() %>% 
  spread(key = location, value = n)%>% 
  replace(is.na(.), 0) %>% 
  mutate(class = case_when(signal == 1 & TMhelix == 0 ~ "sponly",
                           signal == 1 & TMhelix == 1 ~ "type1",
                           signal == 0 & TMhelix  == 1 ~ "type2/3",
                           signal == 0 & TMhelix > 1  ~ "multipassT2/T3",
                           signal == 1 & TMhelix > 1 ~ "multipassT1",.default = "wtf")) %>% 
  select(ENST, class)

classes$class <- factor(classes$class)
summary(classes)

# join the data table with classes
data_with_classes <- inner_join(data, classes, by = "ENST")

# with this data, we should be able to ascertain which of the single pass membrane transcripts are Type II orientation (N terminus: cytosol, C terminus: Lumen) or Type III (N Terminus: Lumen, C terminus: Cytosol).
# Furthermore type III should have preceding N terminal region of < 50 codons#
# Can also delineate tail anchored versus type II 

# filter for single pass transcripts
single_pass <- 
  data_with_classes%>% 
  filter(class == "type2/3")

# split this data into a list
single_pass_split <- split(single_pass, f = single_pass$ENST)

# example observation in list
single_pass_split[[1]]

# create a function to tag a transcript if it is Type II or not
t2_orientation <- function(x){
  
  
  N_term <- 
    x %>% 
    filter(codon_start == 1) %>% 
    pull(var = location)
  
  result <- N_term == "inside"
  
  return(result)
  
}

# vector of transcripts with a type II orientation
type2_orientation <- 
  data.frame(type2 = sapply(single_pass_split, t2_orientation)) %>% 
  rownames_to_column("ENST") %>% 
  filter(type2 == T) %>% 
  pull(ENST)

length(type2_orientation) * 3

TA_transcripts <- 
data_with_classes %>% 
  filter(ENST %in% type2_orientation) %>% 
  select(Gene, ENST, location, codon_start, codon_end, class, nucleotide_sequence_CDS) %>% 
  mutate(CDS_length_nt = nchar(nucleotide_sequence_CDS)) %>% 
  mutate(CDS_length_codons = (CDS_length_nt/3)) %>% 
  select(-nucleotide_sequence_CDS, -CDS_length_nt) %>% 
  mutate(stpcod_dist = CDS_length_codons - codon_end) %>% 
  filter(location == "TMhelix" & stpcod_dist < 60) %>% 
  pull(ENST)

true_T2s <- type2_orientation[!type2_orientation %in% TA_transcripts]

length(TA_transcripts) 
length(true_T2s)
length(type2_orientation)

type3_orientation <- 
  data.frame(type2 = sapply(single_pass_split, t2_orientation)) %>% 
  rownames_to_column("ENST") %>% 
  filter(type2 == F) %>% 
  pull(ENST)

type3_orientation

length(unique(single_pass$ENST)) == length(type3_orientation) + length(true_T2s) + length(TA_transcripts)

# tag with type 2 or type 3 orientation
data_with_classes <- 
  data_with_classes %>% 
  mutate(class = factor(case_when(class == "type2/3" & ENST %in% true_T2s ~ "type2",
                           class == "type2/3" & ENST %in% type3_orientation ~"type3",
                           class == "type2/3" & ENST %in% TA_transcripts ~ "TA",
                           .default = class)))

summary(data_with_classes)

multipass_T2T3 <- data_with_classes %>% filter(class == "multipassT2/T3")
multipass_T2T3_split <- split(multipass_T2T3, multipass_T2T3$ENST)

identify_multipass_subtype <- 
  function(x){
  
  orient <- x %>% filter(codon_start == 1) %>% pull(location)
  torf <- orient == "inside"
  df <- data.frame(ENST = unique(x$ENST),
                   multi_t2_orient = torf)
  
  return(df)
  
}

multipass_T2orT3 <- 
do.call("rbind",lapply(multipass_T2T3_split, identify_multipass_subtype))

rownames(multipass_T2orT3) <- NULL

multipass_T2s <- multipass_T2orT3 %>% filter(multi_t2_orient == T) %>% pull(ENST)
multipass_T3s <- multipass_T2orT3 %>% filter(multi_t2_orient == F) %>% pull(ENST)

summary(data_with_classes)

data_with_classes <- 
data_with_classes %>% 
  mutate(class = factor(case_when(ENST %in% multipass_T2s ~ "multipassT2",
                           ENST %in% multipass_T3s ~ "multipassT3",
                           .default = class)))

summary(data_with_classes$class)

# make column region length
data_with_classes <- data_with_classes %>% mutate(region_length = codon_end - codon_start)

# get rid of nucleotide sequence for now
data_with_classes <- data_with_classes %>% select(-nucleotide_sequence_CDS)

# get the transcript length (in nts) of the transcripts from features
data_with_classes <- data_with_classes %>% inner_join(features %>% select(ENST, cds_length))

# the transcript length is in nts so  convert to codon length. 
#Note: the value of the cds_length is +1 to the final coordinate of codon_end as I think THMHMM has not included the stop codon
data_with_classes$cds_length <- data_with_classes$cds_length/3

data_with_classes %>% filter(class == "TA" & location == "TMhelix")
data_with_classes %>% filter(location != "TMhelix" & codon_end == (cds_length-1))

# none of the codon end coordinates are stop codons
data_with_classes %>% 
  mutate(is_equiv = case_when(codon_end == cds_length~T, 
                              .default = F)) %>% 
  pull(is_equiv) %>% 
  table()

#  split the data
data_split <- split(data_with_classes, data_with_classes$ENST)

sponly_transcript <- 
  (data_with_classes %>% filter(class == "sponly") %>% pull(ENST) %>% unique())[1]
t1_transcript <- 
  (data_with_classes %>% filter(class == "type1") %>% pull(ENST) %>% unique())[1]
t2_transcript <- 
  (data_with_classes %>% filter(class == "type2") %>% pull(ENST) %>% unique())[1]
t3_transcript <- 
  (data_with_classes %>% filter(class == "type3") %>% pull(ENST) %>% unique())[1]
ta_transcript <- 
  (data_with_classes %>% filter(class == "TA") %>% pull(ENST) %>% unique())[1]
multiT1_transcript <- 
  (data_with_classes %>% filter(class == "multipassT1") %>% pull(ENST) %>% unique())[1]
multiT2_transcript <- 
  (data_with_classes %>% filter(class == "multipassT2") %>% pull(ENST) %>% unique())[1]
multiT3_transcript <- 
  (data_with_classes %>% filter(class == "multipassT3") %>% pull(ENST) %>% unique())[1]


sp_ex <- data_split[[{{sponly_transcript}}]]
t1_ex <- data_split[[{{t1_transcript}}]]
t2_ex <- data_split[[{{t2_transcript}}]]
t3_ex <- data_split[[{{t3_transcript}}]]
ta_ex <- data_split[[{{ta_transcript}}]]

multipassT1_ex <- data_split[[{{multiT1_transcript}}]]
multipassT2_ex <- data_split[[{{multiT2_transcript}}]]
multipassT3_ex <- data_split[[{{multiT3_transcript}}]]

summary(data_with_classes$class)




#create function for the TMs
TM_processing <- function(a, dat){
  
  if(unique(dat$class) == "type2" | unique(dat$class) == "type3" | unique(dat$class) == "TA"){
    
    prefix <- dat[a-1,] %>% mutate(region = "prefix")
    TM <- dat[a,] %>% mutate(region = "TS")
    suffix <- dat[a+1,] %>% mutate(region = "suffix")
  } else{
    
    
    # For multi pass membrane proteins
    
    prefix <- dat[a-1,] %>% mutate(codon_start = case_when(
      region_length == 0 ~ codon_start, # This accounts for when there is no prefix 
      region_length %% 2 == 0 ~codon_end - floor((region_length-1)/2), # This accounts for when the linker region is incredibly short
      .default = codon_end - floor(region_length/2)), # This is the default behaviour which splits linkers into 2 equal parts
      region = "prefix")
    TM <- dat[a,] %>% mutate(region = "TS")
    suffix <- dat[a+1,] %>%  
      mutate(codon_end = case_when(region_length == 0 ~ codon_end,
                                   .default = codon_start + floor(region_length/2)),
             region = "suffix")
    
  }
  
  
  
  res <- rbind(rbind(prefix, TM),suffix)
  
  return(res)
  
}

TM_extraction <- function(x){
  
  idx <- which(x$location == "TMhelix")  # obtain indices of the TM domains
  final <- do.call("rbind", lapply(idx,TM_processing, x)) # for each TM domain 
  TMs <- final %>% filter(location == "TMhelix")
  final <-   
    final %>% 
    mutate(ID = rep(1:nrow(TMs), each = 3)) %>% 
    mutate(location = case_when(location == "TMhelix"~"membrane",
                                .default = location)) 
  
  last_observation <- final %>% filter(ID == max(final$ID) & region == "suffix") %>% 
    mutate(codon_end = (cds_length -1))
  
  new_final <- rbind(head(final, n = -1), last_observation)
  
  if(x$class[1] %in% c("multipassT1", "multipassT2")){
    
    new_final <- 
      new_final %>% mutate(codon_start = case_when(ID == 1 & region == "prefix" ~ 1,
                                                   .default = codon_start))
    
  }
  
  return(new_final)
  
}

SP_extraction <- function(dat){
  
  signal <- 
    dat %>% 
    filter(location == "signal") %>% 
    mutate(region = "TS",
           ID = 0)
  
  signal_suffix_idx <- which(dat$location == "signal") +1
  signal_suffix <- 
    dat[signal_suffix_idx,] %>% 
    mutate(codon_end = codon_start + floor(region_length/2),
           region = "suffix",
           ID = 0)
  
  final <- rbind(signal, signal_suffix) %>% 
    mutate(location = case_when(location == "signal" ~"membrane",
                                .default = location))
  
  return(final)
  
}

SP_only <- function(dat){
  
  signal <- 
    dat %>% 
    filter(location == "signal") %>% 
    mutate(region = "TS",
           ID = 0)
  
  signal_suffix_idx <- which(dat$location == "signal") +1
  signal_suffix <- 
    dat[signal_suffix_idx,] %>% 
    mutate(
      region = "suffix",
      ID = 0)
  
  final <- rbind(signal, signal_suffix) %>% 
    mutate(location = case_when(location == "signal" ~"membrane",
                                .default = location))
  
  return(final)
  
}

process_SPonlys <- function(q){
  
  SP <- SP_only(q)
  
  
  fin <- rbind(SP)
  return(fin)
  
}

process_t1 <- function(q){
  
  SP <- SP_extraction(q)
  TM <- TM_extraction(q)
  
  fin <- rbind(SP, TM)
  return(fin)
  
}

process_t2 <- function(q){
  
  fin <- TM_extraction(q)
  
  return(fin)
  
}

process_t3 <- function(q){
  
  fin <- TM_extraction(q)
  
  return(fin)
  
}

process_ta <- function(q){
  
  fin <- TM_extraction(q)
  
  return(fin)
  
}

process_multipass_T2 <- function(q){
  
  fin <- TM_extraction(q)
  
  return(fin)
  
}

process_multipass_T3 <- function(q){
  
  fin <- TM_extraction(q)
  
  return(fin)
  
}

process_multipass_T1 <- function(q){
  
  SP <- SP_extraction(q)
  TM <- TM_extraction(q)
  
  fin <- rbind(SP, TM)
  return(fin)
  
}


process_t1(t1_ex)
process_t2(t2_ex)
process_t3(t3_ex)
process_ta(ta_ex)
#View(process_t4(t4_ex))
process_multipass_T1(multipassT1_ex)
process_multipass_T2(multipassT2_ex)
process_multipass_T2(multipassT3_ex)

process_TMs <- function(x){
  
  if(unique(x$class) == "type1"){
    
    final <- process_t1(x)
    
  } else if(unique(x$class) %in% c("type2", "type3", "TA")){
    
    final <- process_t2(x)
    
  } else if(unique(x$class) %in% c("multipassT2", "multipassT3")) {
    
    final <- process_multipass_T2(x)
    
  } else if(unique(x$class) == "multipassT1"){
    
    final <- process_multipass_T1(x)
    
  } else{
    
    final <- process_SPonlys(x)
  }
  
  return(final)
  
}

data_with_classes$class <- factor(data_with_classes$class) 
summary(data_with_classes$class)

#data_with_classes %>% filter(class == "sp_only")
#transmembranes <- data_with_classes %>% filter(class != "sp_only")

TSs_split <- split(data_with_classes, f = data_with_classes$ENST)
extracted_regions <- do.call("rbind", lapply(TSs_split, process_TMs))

#extracted_regions$codon_start <- ifelse(extracted_regions$codon_start == 1, 2, extracted_regions$codon_start)
extracted_regions[!c(extracted_regions$codon_start <= extracted_regions$codon_end),] # A sanity check

# sanity checks
# 1. Are all codon_start positions less than codon-end positions?
table(extracted_regions$codon_start <= extracted_regions$codon_end)

#2. Are all codon end positions less than the length of the transcript (i.e are stop codons removed?)
table(extracted_regions$codon_end <= extracted_regions$cds_length)

extracted_regions$class <- factor(extracted_regions$class)
summary(extracted_regions$class)

# All TA final suffixes stop 1 codon before stop codon.
extracted_regions %>% 
  filter(class == "TA" & region == "suffix") %>% 
  mutate(is1less = case_when(codon_end == (cds_length-1) ~T,
                             .default = F)) %>% 
  pull(is1less) %>% table()

# Note: initiator ATG still needs to be removed

# rejoin with nucleotide sequences
extracted_regions <- 
  extracted_regions %>% 
  inner_join(master %>% select(ENST, nucleotide_sequence_CDS), by = "ENST")

# create the nt start and stop codons and extract out desired sequence
extracted_regions <- 
  extracted_regions %>% 
  mutate(nt_start = codon_start * 3 -2,
         nt_end = codon_end * 3,
         extracted_seq = substr(nucleotide_sequence_CDS, nt_start, nt_end))

# prefixes
regions_split <- split(extracted_regions, extracted_regions$region)

is.na(regions_split[[1]]$nucleotide_sequence_CDS) %>% table()
regions_split[[1]]


# make function split_into_codons
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


codon_usage <- function(x){
  
  if(x$region[1] == "prefix"){
    start_codon <-  -1
    end_codon <-  -100
  } else if(x$region[1] == "suffix"){
    
    start_codon <- 1
    end_codon <- 100
    
  } else{
    
    start_codon <-  1
    end_codon <-  max(unlist(
      lapply(lapply(x$extracted_seq, 
                    split_into_codons), length)))
  }
  
  
  # split the extracted sequence into codons
  y <- lapply(x$extracted_seq, split_into_codons)
  
  # create unique IDs for each sequence
  ids <- paste0(x$ENST, "_", x$class, "_", x$region, "_", x$location, "_", x$ID, "_", x$codon_start, "_", x$codon_end)
  
  # add these IDs
  names(y) <- ids
  
  if(x$region[1] == "prefix"){
    y <- lapply(y, rev)
  }
  
  
  # make each sequence the same length by appending NAs to the end of the sequence
  seqs_with_NA <- lapply(y, `length<-`, max(lengths(y)))
  
  if(x$region[1] %in% c("prefix", "suffix")){
    
    
    seqs_with_NA <- lapply(seqs_with_NA, head, n = 100)
    
  }
  
  # do call, rbind
  seqs_with_NA <- do.call("rbind", seqs_with_NA)
  
  # make into a data frame 
  codonsN_table <- data.frame(seqs_with_NA)
  
  #max column value
  #max_codon <- dim(codonsN_table)[2]
  
  # make column names nicer
  colnames(codonsN_table) <- c(paste0("codon", start_codon:end_codon))
  
  # make ENST column
  codonsN_table <- codonsN_table %>% rownames_to_column("ENST")
  
  
  # convert table to long format
  codonsN_table_long <- 
    codonsN_table %>% 
    gather(key = "position", 
           value =  "codon", 
           paste0("codon",start_codon):paste0("codon",end_codon))
  
  id_tags <- strsplit(codonsN_table_long$ENST, "_")
  
  codonsN_table_long$ENST <- unlist(lapply(id_tags, function(x){x[[1]]}))
  codonsN_table_long$class <- unlist(lapply(id_tags, function(x){x[[2]]}))
  codonsN_table_long$region <- unlist(lapply(id_tags, function(x){x[[3]]}))
  codonsN_table_long$location <- unlist(lapply(id_tags, function(x){x[[4]]}))
  codonsN_table_long$ID <- unlist(lapply(id_tags, function(x){x[[5]]}))
  codonsN_table_long$codon_start <- unlist(lapply(id_tags, function(x){x[[6]]}))
  codonsN_table_long$codon_end <- unlist(lapply(id_tags, function(x){x[[7]]}))
  
  # convert positions to numerics
  codonsN_table_long$position <- 
    as.numeric(str_remove(codonsN_table_long$position, pattern = "codon"))
  
  codonsN_table_long <- codonsN_table_long %>% arrange(ENST,position)
  
  return(codonsN_table_long)
}


all_data <- do.call("rbind", lapply(regions_split, codon_usage))
rownames(all_data) <- NULL

all_data <- na.omit(all_data)

write.csv(all_data, "C:/Users/Jettles/OneDrive - University of Glasgow/Documents/get_shit_done/TM_intermediate_table.csv")

