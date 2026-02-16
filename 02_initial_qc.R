# script to check for outliers and confirm our most basic assumptions about the data
# note: look at nf-core qc files first

library(variancePartition)
library(DESeq2)
library(PCAtools)
library(ggplot2)
library(stringr)
library(scales)
library(tidyr)
library(dplyr)

library(ComplexHeatmap)
library(circlize)
library(grid)
library(Hmisc)

# prepare data ------------------------------------------------------------
dds <- readRDS("output/dds.rds")

# apply variance stabilizing transformation
vsd <- assay(vst(dds, blind = TRUE))

# get normalized counts
count_matrix <- counts(dds, normalized = TRUE)

# get metadata
metadata <- as.data.frame(colData(dds))

# convert count matrix to long-format data frame for boxplot
long_format_data <- as.data.frame(count_matrix) %>%
  # convert the row names of the data frame to a new column called GENE
  tibble::rownames_to_column(var = "GENE") %>%
  # pivot all columns except GENE to long-format
  tidyr::pivot_longer(cols = -GENE, names_to = "CEGAT_ID", values_to = "COUNT")
# the result is a data frame in long-format, where each row represents a single count for a single gene in a single sample

# filter counts that are zero for plotting
long_format_data <- long_format_data %>% filter(COUNT > 0)

# merge metadata with long_format_data
long_format_data <- long_format_data %>% 
  left_join(metadata, by = c("CEGAT_ID" = "CEGAT_ID"))  # ensure "CEGAT_ID" matches the metadata column

l2fc_data <- readRDS("output/l2fc_data.rds")
l2fc_data_imputed <- readRDS("output/l2fc_data_imputed.rds")

pathway_data <- readRDS("output/pathway_data.rds")

l2fc_pathway_data <- readRDS("output/l2fc_pathway_data.rds")

# check whether STUDY can explain variance in our data --------------------

# let us model STUDY_ID + VISIT
form <- ~ STUDY

# the varPart model enforces us to model either all or none of the categorical variables as random effects
# it is recommended by the developers to model categorical variables with many levels, like inidividual (STUDY_ID) as random effects

varPart <- fitExtractVarPartModel(vsd, form, metadata)

# create violin plot of contribution of each variable to total variance
plotVarPart(varPart)

# save the plot
ggsave("figs/qc/normalized_counts_STUDY_variance_violins.png")

# see also the pca pairs plot in a later sectioned colored by STUDY

# check how much actual days deviate from visit ---------------------------

# convert DAY to numeric
metadata$DAY <- as.numeric(as.character(metadata$DAY))
metadata$NOMINAL_DAY <- as.numeric(as.character(metadata$NOMINAL_DAY))
metadata <- metadata %>% mutate(CEGAT_ID = rownames(.))

for (visit in VISITS) {
  dfv <- metadata %>%
    filter(
      VISIT     == visit,
      DIAGNOSIS %in% names(DIAGNOSIS_COLKEY)
    )
  if (nrow(dfv) == 0) next
  
  nd      <- unique(dfv$NOMINAL_DAY)
  dev_max <- max(abs(dfv$DAY - nd), na.rm = TRUE)
  pad     <- 5
  x_limits <- c(nd - dev_max - pad, nd + dev_max + pad)
  
  p <- ggplot(dfv, aes(y = CEGAT_ID)) +
    # grey line: nominal → actual
    geom_segment(aes(x = NOMINAL_DAY, xend = DAY),
                 colour = "grey80", alpha = 0.5) +
    # # grey anchor at nominal day
    # geom_point(aes(x = NOMINAL_DAY),
    #            colour = "grey40", size = 1.5) +
    # colored point at actual day
    geom_point(aes(x = DAY, colour = DIAGNOSIS),
               size = 2) +
    scale_x_continuous(
      name   = "Day",
      limits = x_limits,
      breaks = pretty(x_limits)
    ) +
    scale_colour_manual(values = DIAGNOSIS_COLKEY) +
    labs(
      title = paste0("Visit ", visit, " (Nominal = ", nd, ")"),
      y     = NULL
    ) +
    theme_bw() +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title   = element_text(hjust = 0.5)
    )
  
  fname <- paste0("figs/qc/DAY_by_DIAGNOSIS_dumbbell_", visit, ".png")
  ggsave(fname, p, width = 6, height = 4, dpi = 300)
  message("Saved: ", fname)
}


# quality control of l2fc data ---------------------------------------------------------------

# find the positions of the NA values in the original matrix
na_positions <- is.na(l2fc_data$count_matrix)

# calculate the total number of values in the matrix
total_values <- length(l2fc_data$count_matrix)

# calculate the number of NA values
num_na_values <- sum(na_positions)

# calculate the percentage of NA values
percentage_na <- round((num_na_values / total_values) * 100, digits = 2)

# print the percentage of NA values
message(paste("Percentage of NA values in log2FoldChange data:", percentage_na, "%"))

# extract the imputed values
imputed_values <- l2fc_data_imputed$count_matrix[na_positions]

# visualize which kind of values were imputed with histogram
p <- ggplot(data.frame(imputed_values), aes(x = imputed_values)) +
  geom_histogram(fill = "lightcoral", color = "darkred", bins = 50) +
  scale_x_continuous(labels = scales::label_number()) +
  labs(title = "Distribution of Imputed Values In Log2FoldChange Data", 
       x = "Imputed Log2FoldChange", 
       y = "Frequency") +
  theme_minimal()

# save histogram of imputed values
ggsave("figs/qc/imputed_values_hist.png", plot = p, width = 8, height = 6)

# PCA and boxplots of normalized count data ---------------------------------------------------------------

# run pca
p <- pca(vsd, metadata = metadata, removeVar = 0.1)

# define variables to color by
colby_vars <- c("STUDY", "AGE_BELOW_60", "VISIT", "REGULAR_MEDICATION", "DIAGNOSIS", "RNA_SAMPLE_INPUT", "RNA_SAMPLE_DNA_DIGESTION", "RNA_SAMPLE_MONTHS_SINCE_VISIT")

# create pairs plot to see whether metadata variables could explain any trends in the data
for (colby in colby_vars) {
  
  # set custom colkey if applicable
  custom_colkey <- NULL
  if (colby == "DIAGNOSIS") {
    custom_colkey <- DIAGNOSIS_COLKEY
  } else if (colby == "VISIT") {
    custom_colkey <- VISIT_COLKEY
  }
  
  png(filename = paste0("figs/qc/COHORT_pca_pairsplot_colby_", colby, ".png"), width = 3600, height = 3600, res = 200)
  p_pairsplot <- pairsplot(p,
                          components = getComponents(p, c(1:10)),
                          triangle = TRUE, trianglelabSize = 12,
                          hline = 0, vline = 0,
                          pointSize = 0.4,
                          gridlines.major = FALSE, gridlines.minor = FALSE,
                          colby = colby,
                          colkey = custom_colkey,
                          title = 'Cohort pca pairs plot of gene counts',
                          plotaxes = FALSE,
                          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))
  
  print(p_pairsplot)  # explicitly render the plot
  dev.off()
}

# save a biplot showing the outlier A_007 and peripheral samples
png(filename = "figs/qc/COHORT_pc1_pc2_pca_biplot_colby_DIAGNOSIS.png", width = 2600, height = 2400, res = 200)
biplot(p,
       lab = p$metadata$SAMPLE_ID,
       x = "PC1",
       y = "PC2",
       colby = "DIAGNOSIS",
       colkey = DIAGNOSIS_COLKEY,
       title = 'Cohort pca biplot of gene counts',
       hline = 0, vline = 0,
       legendPosition = 'right')
dev.off()

# save a biplot for pc3 and pc4 as well
png(filename = "figs/qc/COHORT_pc3_pc4_pca_biplot_colby_DIAGNOSIS.png", width = 2600, height = 2400, res = 200)
biplot(p,
       lab = p$metadata$SAMPLE_ID,
       x = "PC3",
       y = "PC4",
       colby = "DIAGNOSIS",
       colkey = DIAGNOSIS_COLKEY,
       title = 'Cohort pca biplot of gene counts',
       hline = 0, vline = 0,
       legendPosition = 'right')
dev.off()

# loop through colby variables and save boxplots
for (colby in colby_vars) {
  
  # base plot
  boxplot_gene <- ggplot(long_format_data, aes(x = CEGAT_ID, y = COUNT, fill = .data[[colby]])) +
    geom_boxplot(outlier.shape = 16) +
    scale_y_log10() +
    labs(
      title = paste("Normalized Gene Count Distribution per Sample - Colored by", colby),
      x = "Sample",
      y = "Gene Count (log10 scale)",
      fill = str_to_title(colby)  # Legend title
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  # apply custom colors if colby is DIAGNOSIS or VISIT
  if (colby == "DIAGNOSIS") {
    boxplot_gene <- boxplot_gene + scale_fill_manual(values = DIAGNOSIS_COLKEY)
  } else if (colby == "VISIT") {
    boxplot_gene <- boxplot_gene + scale_fill_manual(values = VISIT_COLKEY)
  }
  
  # save the boxplot
  ggsave(paste0("figs/qc/COHORT_normalized_genecount_boxplot_colby_", colby, ".pdf"),
         plot = boxplot_gene, limitsize = FALSE, width = 80, height = 8)
}


# PCA of log2FoldChange count data ----------------------------------------

# run PCA
p <- pca(l2fc_data_imputed$count_matrix, center = TRUE, scale = TRUE, removeVar = 0.1, metadata = l2fc_data_imputed$metadata)

# define variables to color by
colby_vars <- c("STUDY", "AGE_BELOW_60", "VISIT", "DIAGNOSIS", "RNA_SAMPLE_INPUT", "RNA_SAMPLE_DNA_DIGESTION", "RNA_SAMPLE_MONTHS_SINCE_VISIT")

# create pairs plot to see whether metadata variables could explain any trends in the data
for (colby in colby_vars) {
  png(filename = paste0("figs/qc/COHORT_l2fc_pca_pairsplot_colby_", colby, ".png"), width = 3600, height = 3600, res = 200)
  
  # set custom colkey if applicable
  custom_colkey <- NULL
  if (colby == "DIAGNOSIS") {
    custom_colkey <- DIAGNOSIS_COLKEY
  } else if (colby == "VISIT") {
    custom_colkey <- VISIT_COLKEY
  }
  
  # plot
  p_pairsplot <- pairsplot(p,
                           components = getComponents(p, c(1:10)),
                           triangle = TRUE, trianglelabSize = 12,
                           hline = 0, vline = 0,
                           pointSize = 0.4,
                           gridlines.major = FALSE, gridlines.minor = FALSE,
                           colby = colby,
                           colkey = custom_colkey,
                           title = 'Cohort pairs plot of V1vsVX log2FoldChange gene counts',
                           plotaxes = FALSE,
                           margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))
  
  print(p_pairsplot)
  dev.off()
}


# sample distances heatmap --------------------------------------------------------

##### euclidean distance #####

# check if euclidean distances exist
euclid_dists_file <- "output/euclid_dists.rds"
if (file.exists(euclid_dists_file)) {
  # load existing results
  euclid_dists <- readRDS(euclid_dists_file)
  message("Loaded existing sample euclidean distances from file.")
} else {
  message("Calculating euclidean distances between samples")
  euclid_dists <- list()
  # to ensure roughly equal contribution from all genes, we will use variance stabilized counts
  # distances can be calculated with R's dist() function and we transpose, because the function expects rows to be samples
  euclid_dists$matrix <- as.matrix(dist(t(vsd)))
  euclid_dists$metadata <- metadata
  
  # save euclid_dists
  saveRDS(euclid_dists, file = euclid_dists_file)
}

# convert euclid dists to matrix
euclid_dists_matrix <- as.matrix(euclid_dists)

# function to visualize euclidean distances
plot_sample_distance_heatmap <- function(dists_matrix,
                                         metadata,
                                         diagnoses,                   # Either a vector of diagnoses or a single diagnosis (e.g., "NSCLC")
                                         visits,                      # e.g., "V2" or c("V1", "V2")
                                         sample_set,                  # A fixed vector of samples to plot (e.g., c("", "", ""))
                                         cluster_by = "rows",         # Accepts: "rows", "columns", or "both"
                                         n_row_clusters = 4,          # Number of row clusters
                                         n_column_clusters = 4,       # Number of column clusters
                                         additional_annotation_col = NULL, # Name of an extra metadata column (e.g., "AGE")
                                         value_name = "euclidean distance" # Name of the value that is plotted and described in the legend
) {
  
  # required packages
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  if (!is.null(additional_annotation_col)) {
    library(RColorBrewer)
  }
  
  ## subset data ##
  
  # use %in% to allow a vector of visits.
  metadata <- metadata %>% 
    filter(VISIT %in% visits) %>%
    filter(DIAGNOSIS %in% diagnoses) %>%
    droplevels() %>%                      # drop unused factor levels
    mutate(CEGAT_ID = as.character(CEGAT_ID), # convert to character to avoid problems with factor levels
           STUDY_ID = as.character(STUDY_ID)) %>%
    arrange(DIAGNOSIS, match(CEGAT_ID, colnames(dists_matrix)))
  
  # subset dists_matrix to include only columns corresponding to filtered metadata
  dists_matrix <- dists_matrix[metadata$CEGAT_ID, metadata$CEGAT_ID, drop = FALSE]
  
  # create a named vector mapping from CEGAT_IDs to SAMPLE_IDs
  id_mapping <- setNames(as.character(metadata$SAMPLE_ID), metadata$CEGAT_ID)
  new_colnames <- id_mapping[colnames(dists_matrix)]
  
  # update matrix column/row names
  colnames(dists_matrix) <- new_colnames
  rownames(dists_matrix) <- new_colnames
  
  ## build the top annotation data frame ##
  
  # use the actual visits from metadata for annotation
  ann_df <- data.frame(Visit = metadata$VISIT,
                       stringsAsFactors = FALSE)
  
  # drop any unused factor levels
  ann_df$Visit <- factor(as.character(ann_df$Visit),
                         levels = unique(as.character(ann_df$Visit)))
  
  if (length(diagnoses) == 1) {
    ann_df$Diagnosis <- factor(rep(diagnoses[1], nrow(metadata)), levels = diagnoses)
  } else {
    ann_df$Diagnosis <- factor(metadata$DIAGNOSIS, levels = unique(metadata$DIAGNOSIS))
  }
  
  if (!is.null(additional_annotation_col)) {
    if (!(additional_annotation_col %in% colnames(metadata))) {
      stop(paste("Column", additional_annotation_col, "not found in metadata."))
    }
    ann_df[[additional_annotation_col]] <- metadata[[additional_annotation_col]]
  }
  
  ## build color mapping for annotations ##
  ann_colors <- list()
  
  ann_colors[["Visit"]] <- VISIT_COLKEY
  
  ann_colors[["Diagnosis"]] <- DIAGNOSIS_COLKEY
  
  if (!is.null(additional_annotation_col)) {
    extra_vals <- ann_df[[additional_annotation_col]]
    if (is.numeric(extra_vals)) {
      ann_colors[[additional_annotation_col]] <- colorRamp2(
        c(min(extra_vals, na.rm = TRUE), max(extra_vals, na.rm = TRUE)),
        c("blue", "red")
      )
    } else {
      extra_levels <- unique(as.character(extra_vals))
      n_levels <- length(extra_levels)
      if (n_levels < 3) {
        extra_pal <- brewer.pal(3, "Set1")[seq_len(n_levels)]
      } else {
        extra_pal <- brewer.pal(n_levels, "Set1")
      }
      ann_colors[[additional_annotation_col]] <- setNames(extra_pal, extra_levels)
    }
  }
  
  top_annotation <- HeatmapAnnotation(
    df = ann_df,
    col = ann_colors,
    annotation_name_gp = gpar(fontsize = 10),
    annotation_legend_param = list(
      Diagnosis = list(title = "Diagnosis", 
                       at = levels(ann_df$Diagnosis), 
                       labels = levels(ann_df$Diagnosis))
    )
  )
  
  
  ## set up the heatmap color function ##

  col_fun <- colorRamp2(
    c(min(dists_matrix, na.rm = TRUE), 0, max(dists_matrix, na.rm = TRUE)),
    c("blue", "white", "red")
  )
  
  
  ## clustering and heatmap creation with split clusters ##

  ht <- NULL
  clusters <- list()
  
  if (cluster_by == "rows") {
    row_dist <- dist(dists_matrix, method = "euclidean")
    row_clust <- hclust(row_dist, method = "ward.D2")
    row_order <- seriation::get_order(seriation::seriate(row_dist, method = "OLO"))
    clusters_row <- cutree(row_clust, k = n_row_clusters)
    
    ht <- Heatmap(
      dists_matrix,
      name = value_name,
      col = col_fun,
      
      cluster_rows = TRUE, # perform clustering internally
      clustering_distance_rows = function(m) dist(m, method = "euclidean"),
      clustering_method_rows = "ward.D2",
      row_order = row_order,
      row_split = clusters_row,        # split rows into n_row_clusters
      row_dend_reorder = FALSE,
      
      cluster_columns = FALSE,
      show_column_names = TRUE,
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 6),
      top_annotation = top_annotation,
      row_dend_side = "right",
      row_dend_width = unit(3, "cm"),
      show_parent_dend_line = TRUE,
      use_raster = TRUE
    )
    clusters$row_clusters <- clusters_row
    
  } else if (cluster_by == "columns") {
    col_dist <- dist(t(dists_matrix), method = "euclidean")
    col_clust <- hclust(col_dist, method = "ward.D2")
    col_order <- seriation::get_order(seriation::seriate(col_dist, method = "OLO"))
    clusters_col <- cutree(col_clust, k = n_column_clusters)
    
    ht <- Heatmap(
      dists_matrix,
      name = value_name,
      col = col_fun,
      cluster_rows = FALSE,
      
      cluster_columns = TRUE,         # perform clustering internally 
      clustering_distance_columns = function(m) dist(m, method = "euclidean"),
      clustering_method_columns = "ward.D2",
      column_order = col_order,
      column_split = clusters_col,     # split columns into n_column_clusters
      column_dend_reorder = FALSE,
      
      show_column_names = TRUE,
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 6),
      top_annotation = top_annotation,
      row_dend_side = "right",
      row_dend_width = unit(3, "cm"),
      show_parent_dend_line = TRUE,
      use_raster = TRUE
    )
    clusters$column_clusters <- clusters_col
    
  } else if (cluster_by == "both") {
    # Row clustering
    row_dist <- dist(dists_matrix, method = "euclidean")
    row_clust <- hclust(row_dist, method = "ward.D2")
    row_order <- seriation::get_order(seriation::seriate(row_dist, method = "OLO"))
    clusters_row <- cutree(row_clust, k = n_row_clusters)
    
    # Column clustering
    col_dist <- dist(t(dists_matrix), method = "euclidean")
    col_clust <- hclust(col_dist, method = "ward.D2")
    col_order <- seriation::get_order(seriation::seriate(col_dist, method = "OLO"))
    clusters_col <- cutree(col_clust, k = n_column_clusters)
    
    ht <- Heatmap(
      dists_matrix,
      name = value_name,
      col = col_fun,
      
      cluster_rows = TRUE, # perform clustering internally
      clustering_distance_rows = function(m) dist(m, method = "euclidean"),
      clustering_method_rows = "ward.D2",
      row_order = row_order,
      row_split = clusters_row,        # split rows into n_row_clusters
      row_dend_reorder = FALSE,
      
      cluster_columns = TRUE, # perform clustering internally
      clustering_distance_columns = function(m) dist(m, method = "euclidean"),
      clustering_method_columns = "ward.D2",
      column_order = col_order,
      column_split = clusters_col,     # split columns into n_column_clusters
      column_dend_reorder = FALSE,
      
      show_column_names = TRUE,
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 6),
      top_annotation = top_annotation,
      row_dend_side = "right",
      row_dend_width = unit(3, "cm"),
      show_parent_dend_line = TRUE,
      use_raster = TRUE
    )
    clusters$row_clusters <- clusters_row
    clusters$column_clusters <- clusters_col
    
  } else {
    stop("Invalid cluster_by parameter. Please choose 'rows', 'columns', or 'both'.")
  }
  
  return(list(ht = ht, clusters = clusters))
}


sample_distance_heatmap <- plot_sample_distance_heatmap(euclid_dists$matrix,
                                                        euclid_dists$metadata,
                                                        diagnoses = c("HC", "NSCLC", "SMM"),
                                                        visits = c("V3", "V4"),
                                                        cluster_by = "both",
                                                        n_column_clusters = 12,
                                                        n_row_clusters = 12,
                                                        additional_annotation_col = "AGE")
sample_distance_heatmap$ht

# save the heatmap as a PNG file
# width = 3400, height = 2600 for single diagnosis
# width = 6800, height = 5200 for cohort single visit
# width = 18600, height = 16400 for cohort all visits
png("figs/qc/COHORT_sample_euclidean_distance_heatmap.png", width = 18600, height = 16400, res = 300) # adjust width, height, and resolution
draw(sample_distance_heatmap$ht)
dev.off()

# save heatmap as pdf file
pdf("figs/qc/COHORT_sample_euclidean_distance_heatmap.pdf", width = 46, height = 46)  # width and height in inches
draw(sample_distance_heatmap$ht)
dev.off()

# checking potential sample contamination with sex genes ---------------------

# S11947Nr54 is female
# S11947Nr55 is male

# define your sex genes of interest
sex_genes <- c("XIST", "USP9Y", "UTY")

sexes <- c("f", "m")

for (sex in sexes) {
  
  # subset the counts for the current sex
  counts_sex <- as.data.frame(count_matrix[sex_genes, dds$SEX == sex])
  
  # convert the normalized counts matrix into long format for boxplot
  counts_long <- counts_sex %>%
    tibble::rownames_to_column("GENE") %>%
    pivot_longer(cols = -GENE, names_to = "CEGAT_ID", values_to = "COUNT")
  
  # flag the samples we want to highlight (based on CEGAT_ID)
  if (sex == "f") {
    
    # highlight samples
    counts_long <- counts_long %>%
      mutate(HIGHLIGHT = ifelse(CEGAT_ID == "S11947Nr54", "S11947Nr54", "other samples"))
    
    # create boxplot
    p <- ggplot(counts_long, aes(x = GENE, y = COUNT)) +
      geom_boxplot(outlier.shape = NA) +  # draw boxplot without showing the default outliers
      geom_jitter(aes(color = HIGHLIGHT), width = 0.2, size = 2) +  # Overlay individual points
      scale_color_manual(values = c("S11947Nr54" = "red", "S11947Nr55" = "blue", "other samples" = "black")) +
      labs(title = "Normalized Counts for Sex Genes in FEMALE COHORT",
           x = "Gene",
           y = "Normalized Count") +
      theme_bw()
    
    # save plot
    ggsave(paste0("figs/qc/FEMALE_sex_genes_genecount_boxplot_colby_.png"), plot = p, width = 8, height = 8)
    
  } else {
    
    # highlight samples
    counts_long <- counts_long %>%
      mutate(HIGHLIGHT = ifelse(CEGAT_ID == "S11947Nr55", "S11947Nr55", "other samples"))
    
    # create boxplot
    p <- ggplot(counts_long, aes(x = GENE, y = COUNT)) +
      geom_boxplot(outlier.shape = NA) +  # draw boxplot without showing the default outliers
      geom_jitter(aes(color = HIGHLIGHT), width = 0.2, size = 2) +  # Overlay individual points
      scale_color_manual(values = c("S11947Nr54" = "red", "S11947Nr55" = "blue", "other samples" = "black")) +
      labs(title = "Normalized Counts for Sex Genes in MALE COHORT",
           x = "Gene",
           y = "Normalized Count") +
      theme_bw()
    
    # save plot
    ggsave(paste0("figs/qc/MALE_sex_genes_genecount_boxplot_colby_.png"), plot = p, width = 8, height = 8)
  }
}

# these plots indicate that S11947Nr54 is not contaminated with S11947Nr55 or vice versa, we keep both samples

# excluding specific samples ---------------------

# our QC has revealed that A_007 is a potential outlier (PCA biplot and euclidean distance heatmap). We could remove it from the subsequent analysis to denoise our effects of interest

# while some C subjects, like C_008, have quite a different gene expression profile from all other samples, we will keep them, because we suspect a biological effect of Chemotherapy

# subject B_107 has no V1 RNAseq sample. We will include the other samples from this subject whenever we work with absolute counts, but remove it whenever looking at comparisons with V1 as the reference

# subject B_016 receives isatuximab, which could influence the immune system a lot, we will remove this sample because it is more likely to be just a noise factor (could we check this somehow?)


# excluded samples (after qc results)
EXCLUDED_SUBJECTS = c("SYS01_A_007", "SYS01_B_016")

# update dds, log2FoldChange data, distance data and pathway data
for (subject in EXCLUDED_SUBJECTS) {
  
  dds <- subset_dds(dds, condition = paste0("STUDY_ID !=", "'", subject, "'"))

  l2fc_data <- subset_l2fc_data(l2fc_data, condition = paste0("STUDY_ID !=", "'", subject, "'"))

  l2fc_data_imputed <- subset_l2fc_data(l2fc_data_imputed, condition = paste0("STUDY_ID !=", "'", subject, "'"))

  euclid_dists <- subset_euclid_dists(euclid_dists, condition = paste0("STUDY_ID !=", "'", subject, "'"))
  
  pathway_data <- subset_pathway_data(pathway_data, condition = paste0("STUDY_ID !=", "'", subject, "'"))
  
  l2fc_pathway_data <-subset_pathway_data(l2fc_pathway_data, condition = paste0("STUDY_ID !=", "'", subject, "'"))
}

EXCLUDED_SAMPLES = c("SYS03_B_107_V1") # should not be present in the first place, maybe in metadata

for (sample in EXCLUDED_SAMPLES) {
  
  # exclude sample from dds
  dds <- subset_dds(dds, condition = paste0("SAMPLE_ID !=", "'", sample, "'"))

  # exclude sample from l2fc data
  l2fc_data <- subset_l2fc_data(l2fc_data, condition = paste0("SAMPLE_ID !=", "'", sample, "'"))

  l2fc_data_imputed <- subset_l2fc_data(l2fc_data_imputed, condition = paste0("SAMPLE_ID !=", "'", sample, "'"))

  euclid_dists <- subset_euclid_dists(euclid_dists, condition = paste0("SAMPLE_ID !=", "'", sample, "'"))
  
  pathway_data <- subset_pathway_data(pathway_data, condition = paste0("SAMPLE_ID !=", "'", sample, "'"))
  
  l2fc_pathway_data <-subset_pathway_data(l2fc_pathway_data, condition = paste0("SAMPLE_ID !=", "'", sample, "'"))
}

saveRDS(dds, "output/dds_qc.rds")
saveRDS(l2fc_data, "output/l2fc_data_qc.rds")

# ideally, one would perform the imputation after the exclusion, but it won't change much
saveRDS(l2fc_data_imputed, "output/l2fc_data_imputed_qc.rds")

saveRDS(euclid_dists, "output/euclid_dists_qc.rds")

saveRDS(pathway_data, "output/pathway_data.rds")
saveRDS(l2fc_pathway_data, "output/l2fc_pathway_data.rds")


message("Quality control and exclusion of samples completed.")