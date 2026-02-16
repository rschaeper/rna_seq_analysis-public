# script to retrieve metadata and count data for each subject in SYS01 and SYS03
# outputs dds.rds, dds_subsets.rds, l2fc_data.rds, l2fc_data_imputed.rds, l2fc_imputed_subsets.rds

# Initialization
library(dplyr)
library(stringr)
library(openxlsx)
library(ggplot2)
library(scales)  # for controlling axis labels
library(impute)
library(DESeq2)

# load count matrix and metadata ---------------------------------------------------

# read in count data for SYS01
count_matrix_sys01 <- read.csv("data/counts/SYS01_salmon.merged.gene.counts.tsv", sep = "\t")
# read in count data for SYS03
count_matrix_sys03 <- read.csv("data/counts/SYS03_salmon.merged.gene.counts.tsv", sep = "\t")

# join them
count_matrix <- inner_join(count_matrix_sys01, count_matrix_sys03, by = c("gene_id", "gene_name"))

# note: 61 gene IDs differ from their respective gene name
# for now, we will use the gene IDs as rows, as duplicate gene names are not allowed as rows anyways

# convert gene_id column to rownames
rownames(count_matrix) <- count_matrix$gene_id

# remove gene_id and gene_name columns
count_matrix <- count_matrix[, -c(1, 2)]

# counts are non-integer values because Salmon does not align on base-by-base directly, but performs so-called pseudoalignment
# if a transcript detected could be ambiguously mapped to multiple isoforms or genes because of similarity, a read is fractionally distributed using an
# expectation maximalization algorithm

# round counts to nearest integer value
# uses "banker's rounding", meaning we round to the nearest even integer in the case of "0.5"
count_matrix <- round(count_matrix)

head(count_matrix)

saveRDS(count_matrix, "output/raw_count_matrix_rounded.rds")

metadata <- readRDS("output/metadata.rds")

# construct DESeq2 dataset -------------------------------------------------------------

# reorder metadata rows (CeGatIDs) to match the order of CeGatID columns in the count_matrix
metadata <- metadata[match(colnames(count_matrix), metadata$CEGAT_ID), ]
# note: the line of code above will remove all the samples for which we have no count data

# for SYS01 we have 295 recorded samples, but only obtained data for 259 samples (technical problems at CeGaT)
# for SYS03 we have 513 recorded samples, but only obtained data for 508 samples (inadequate RNA amounts)
# B_016 is receiving isatuximab, which targets CD38 and will likely affect vaccination response in a drastic way. We will exclude this subject from the analysis for now.
# B_104 is missing V1. We will use the other available samples for all analyses except for the subject- and gene-wise log2FoldChange related plots
# this leaves us with a total of 767 samples

# verify order
if (!all(colnames(count_matrix) == metadata$CEGAT_ID)) {
  stop("Error: Count matrix columns and metadata rows are not in the same order for dds!")
}

message("Count matrix columns and metadata rows are in order for dds.") # returns true if successfully ordered

# DESeqDataSet inherits from RangedSummarizedExperiment, so we can easily subset metadata and the count_matrix in sync

# initialize DESeq Data Set Object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ STUDY_ID + VISIT) # specify a generic design here, but be aware that it is used for vst or rlog when blind=FALSE

# set reference level for VISIT
dds$VISIT <- relevel(dds$VISIT, ref = "V1")

# set reference level for DIAGNOSIS
dds$DIAGNOSIS <- relevel(dds$DIAGNOSIS, ref = "HC")

message("Filtering genes not reaching a count of 10 or more in at least 5 samples")

# require genes to have a count of 10 or more in at least 5 samples to be included in the analysis
keep <- rowSums(counts(dds) >= 10) >= 5

# log the discarded genes
kept_genes <- rownames(dds)[keep]
discarded_genes <- rownames(dds)[!keep]

message("Number of genes kept for subsequent analysis: ", length(kept_genes))
message("Number of genes discarded before analysis: ", length(discarded_genes))

# subset dds genes
dds <- dds[keep,]

# apply normalization for RNA library size
# we do this here, because we will (and should) pretty much always work with data adjusted by library size
dds <- estimateSizeFactors(dds)

# save dds
saveRDS(dds, file = "output/dds.rds")

# calculate subject- and gene-wise log2FoldChange matrix with V1 as reference ------------------

# function to calculate log2FoldChanges on a sample level
calculate_sample_log2foldchanges <- function(count_matrix, metadata) {
  
  common_samples <- intersect(colnames(count_matrix), rownames(metadata))
  count_matrix <- count_matrix[, common_samples]
  metadata <- metadata[common_samples, ]
  
  # get all samples (including V1)
  all_samples <- rownames(metadata)
  
  # create an empty matrix to store the log2 fold changes
  log2foldchange_matrix <- matrix(NA_real_, nrow = nrow(count_matrix), ncol = length(all_samples)) # use NA_real_ for numeric NA
  rownames(log2foldchange_matrix) <- rownames(count_matrix)
  colnames(log2foldchange_matrix) <- all_samples
  
  # loop through each gene
  for (gene in rownames(count_matrix)) {
    
    # loop through each sample
    for (sample in all_samples) {
      
      # get the subject ID for the current sample
      subject <- metadata[sample, "STUDY_ID"]
      
      # find the corresponding V1 sample for the subject
      v1_sample <- rownames(metadata[metadata$STUDY_ID == subject & metadata$VISIT == "V1", ])
      
      # check if the sample is V1
      if(metadata[sample, "VISIT"] == "V1"){
        log2foldchange_matrix[gene, sample] <- 0 # V1 is reference, so log2FC is 0
      } else {
        # check if both samples exist
        if (length(v1_sample) == 1) {
          
          # get counts
          v1_count <- count_matrix[gene, v1_sample]
          sample_count <- count_matrix[gene, sample]
          
          # calculate log2 fold change, avoid -Inf
          if (v1_count == 0 && sample_count == 0) {
            log2fc <- 0 # if both V1 and VX are 0, log2FC is 0
          } else if (v1_count == 0 || sample_count == 0) {
            log2fc <- NA_real_ # if either V1 or VX are 0, log2FC is NA
            # we can then later set these NA values to an arbitrarily low value or the gene-wise median or apply clever imputation
          } else {
            log2fc <- log2(sample_count / v1_count)
          }
          
          # store result
          log2foldchange_matrix[gene, sample] <- log2fc
        } else {
          # if the V1 sample doesn't exist, set the log2fc to NA
          log2foldchange_matrix[gene, sample] <- NA_real_
        }
      }
    }
  }
  return(log2foldchange_matrix)
}

# check if l2fc_data exists
l2fc_data_file <- "output/l2fc_data.rds"

if (file.exists(l2fc_data_file)) {
  
  # load existing results
  l2fc_data <- readRDS(l2fc_data_file)
  
  l2fc_metadata <- l2fc_data$metadata
  
  l2fc_count_matrix <- l2fc_data$count_matrix
  
  message(paste0("Log2FoldChange data already exists at: ", l2fc_data_file))
  
} else {
  
  message("Calculating Log2FoldChange count matrix.")
  
  metadata <- as.data.frame(colData(dds))
  
  # make sure to calculate subject- and gene-wise log2FoldChanges with counts that were normalized (for RNA library size)
  normalized_count_matrix <- as.data.frame(counts(dds, normalized = TRUE))
  
  # run l2fc calc
  l2fc_count_matrix <- calculate_sample_log2foldchanges(normalized_count_matrix, metadata)
  
  # find columns with all NA values, those were missing V1 (donor B_107)
  all_na_cols <- colSums(is.na(l2fc_count_matrix)) == nrow(l2fc_count_matrix)
  
  # get the names of the columns with all NA values
  na_col_names <- colnames(l2fc_count_matrix)[all_na_cols]
  
  # remove the columns with all NA values
  l2fc_count_matrix <- l2fc_count_matrix[, !all_na_cols]
  
  # remove the corresponding rows from the metadata
  l2fc_metadata <- metadata[!metadata$CEGAT_ID %in% na_col_names, ]
  
  # relevel l2fc_metadata
  #l2fc_metadata <- as.data.frame(l2fc_metadata) # make sure it is a dataframe
  for (col in colnames(l2fc_metadata)) {
    # relevel all factor columns
    if (is.factor(l2fc_metadata[[col]])) {
      l2fc_metadata[[col]] <- droplevels(l2fc_metadata[[col]])
    }
  }
  
  # print the names of the removed columns
  message("Samples removed from l2fc_count_matrix and metadata due to all NA values: ", paste(na_col_names, collapse = ", "))
  
  # pair count data and associated metadata into one R object if they match
  if (all.equal(colnames(l2fc_count_matrix), as.character(l2fc_metadata$CEGAT_ID))) {
    
    l2fc_data <- list(count_matrix = l2fc_count_matrix, metadata = l2fc_metadata)
    
    # calculating this l2fc_count_matrix takes really long, save it for later!
    saveRDS(l2fc_data, l2fc_data_file)
    
  } else {
    
    # error
    stop("Error: Column names of l2fc_count_matrix do not match l2fc_metadata$CEGAT_ID!")
  }
}

# check if l2fc_data_imputed exists
l2fc_data_imputed_file <- "output/l2fc_data_imputed.rds"

if (file.exists(l2fc_data_imputed_file)) {
  
  # load existing results
  l2fc_data_imputed <- readRDS(l2fc_data_imputed_file)
  
  l2fc_count_matrix_imputed <- l2fc_data_imputed$count_matrix
  
  message(paste0("Log2FoldChange data with imputation applied already exists at: ", l2fc_data_imputed_file))
  
} else {
  
  message("Imputing NA values in log2FoldChange data.")
  
  # impute NA data
  l2fc_count_matrix_imputed <- impute.knn(l2fc_data$count_matrix)$data
  
  l2fc_data_imputed <- list(count_matrix = l2fc_count_matrix_imputed, metadata = l2fc_metadata)
  
  # save the imputed matrix
  saveRDS(l2fc_data_imputed, l2fc_data_imputed_file)
}

message("Data extraction and preparation completed.")
# score pathways subject-wise using GSVA ---------------------------------------------

pathway_data_file <- "output/pathway_data.rds"

if (file.exists(pathway_data_file)) {
  
  message(paste0("Subject-level pathway data already exists at: ", pathway_data_file))
  
  #pathway_data <- readRDS(pathway_data_file)
  
} else {
  
  param <- gsvaParam(exprData = as.matrix(count_matrix),
                     geneSets = PATHWAYS)
  
  gsva_matrix <- gsva(param)
  
  pathway_data <- list()
  
  pathway_data$gsva_matrix <- gsva_matrix
  pathway_data$metadata <- metadata
  
  saveRDS(pathway_data, pathway_data_file)
  
}

# score l2fc pathways subject-wise using GSVA ---------------------------------------------

l2fc_pathway_data_file <- "output/l2fc_pathway_data.rds"

if (file.exists(l2fc_pathway_data_file)) {
  
  message(paste0("Subject-level l2fc pathway data already exists at: ", l2fc_pathway_data_file))
  
  #l2fc_pathway_data <- readRDS(l2fc_pathway_data_file)
  
} else {
  
  param <- gsvaParam(exprData = as.matrix(l2fc_data_imputed$count_matrix),
                     geneSets = PATHWAYS)
  
  gsva_matrix <- gsva(param)
  
  l2fc_pathway_data <- list()
  
  l2fc_pathway_data$gsva_matrix <- gsva_matrix
  l2fc_pathway_data$metadata <- l2fc_data_imputed$metadata
  
  saveRDS(l2fc_pathway_data, l2fc_pathway_data_file)
}

