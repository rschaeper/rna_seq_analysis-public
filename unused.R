# script to store parts of the analysis that were not used in the final thesis

#### Quality Control ####
# unused quality control libraries ----------------------------------------
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

# unused quality control prepare data ------------------------------------------------------------
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


# checking whether IgG levels are correlated with STUDY within diagnosis groups -------------------
# 1) Recode STUDY to −1/+1
metadata$STUDY_num <- ifelse(metadata$STUDY == "SYS01", -1, 1)

# 2) Compute Spearman ρ & p-values across DIAGNOSIS x VISIT
igg_variables <- c("IgG_V1", "IgG_V3", "IgG_V4")
diagnoses     <- unique(metadata$DIAGNOSIS)

corr_mat <- matrix(
  NA_real_,
  nrow = length(diagnoses),
  ncol = length(igg_variables),
  dimnames = list(diagnoses, igg_variables)
)
p_mat    <- corr_mat

for (igg_variable in igg_variables) {
  visit <- sub("^IgG_", "", igg_variable)
  for (d in diagnoses) {
    idx <- which(
      metadata$DIAGNOSIS == d &
        metadata$VISIT     == visit &
        !is.na(metadata[[igg_variable]])
    )
    if (length(idx) < 3) next
    
    ct <- cor.test(
      metadata$STUDY_num[idx],
      metadata[[igg_variable]][idx],
      method = "spearman",
      exact  = FALSE
    )
    corr_mat[d, igg_variable] <- unname(ct$estimate)
    p_mat[d,    igg_variable] <- ct$p.value
  }
}

# 3) Build a diverging color function from -1 to 0 to +1
col_fun <- colorRamp2(
  breaks = c(-1, 0, 1),
  colors = c("blue", "white", "red")
)

# 4) Draw heatmap with fixed scale
png("figs/qc/IgG_STUDY_correlation_within_DIAGNOSIS_heatmap.png", width = 2400, height = 2400, res = 300)
Heatmap(
  corr_mat,
  name       = "Spearman ρ\n(+1=SYS03,-1=SYS01)",
  col        = col_fun,
  na_col     = "grey90",
  row_title     = "Diagnosis",
  column_title  = "Correlation of IgG with STUDY within diagnosis subsets",
  heatmap_legend_param = list(
    at     = seq(-1, 1, by = 0.5),
    labels = seq(-1, 1, by = 0.5),
    title_gp  = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 8)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(p_mat[i, j]) && p_mat[i, j] < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 14, col = "black"))
    }
  }
)
dev.off()


##### poisson distance #####
library("PoiClaClu")

# check if poisson distances exist
pois_dists_file <- "output/pois_dists.rds"
if (file.exists(pois_dists_file)) {
  # load existing results
  pois_dists <- readRDS(pois_dists_file)
  message("Loaded existing sample poisson distances from file.")
} else {
  message("Calculating poisson distances between samples")
  # run pois
  pois_dists <- list()
  # use original counts because poisson takes the inherent variance structure into account
  # we transpose, because the function expects rows to be samples
  pois_dists$matrix <- as.matrix(PoissonDistance(t(counts(dds)))$dd)
  rownames(pois_dists$matrix) <- metadata$CEGAT_ID
  colnames(pois_dists$matrix) <- metadata$CEGAT_ID
  
  pois_dists$metadata <- metadata
  # save pois
  saveRDS(pois_dists, file = pois_dists_file)
}



sample_distance_heatmap <- plot_sample_distance_heatmap(pois_dists$matrix,
                                                        pois_dists$metadata,
                                                        diagnoses = DIAGNOSES,
                                                        visits = VISITS,
                                                        cluster_by = "both",
                                                        n_column_clusters = 12,
                                                        n_row_clusters = 12,
                                                        additional_annotation_col = "AGE",
                                                        value_name = "poisson distance")
sample_distance_heatmap$ht

# save the heatmap as a PNG file
# width = 3400, height = 2600 for single diagnosis
# width = 6800, height = 5200 for cohort single visit
# width = 18600, height = 16400 for cohort all visits
png("figs/qc/COHORT_sample_poisson_distance_heatmap.png", width = 18600, height = 16400, res = 300) # adjust width, height, and resolution
draw(sample_distance_heatmap$ht)
dev.off()

# save heatmap as pdf file
pdf("figs/qc/COHORT_sample_poisson_distance_heatmap.pdf", width = 30, height = 28)  # width and height in inches
draw(sample_distance_heatmap$ht)
dev.off()


#### Exploratory Analysis ####
# unused exploratory analysis libraries --------------------------------------------

# load required libraries
library(ggplot2)
library(ggpubr)
#library(scales)  # needed for pretty_breaks()
library(PCAtools)
library(treemap)
library(variancePartition)
library(dplyr)
library(impute)
library(reshape2)
library(ComplexHeatmap)
library(circlize)        # for colorRamp2()
library(grid)            # for unit()
library(GSVA)
library(DESeq2)

# MineICA and dependencies
library(Biobase)
library(plyr)
library(foreach)
library(xtable)
library(biomaRt)
library(GOstats)
library(cluster)
library(marray)
library(mclust)
library(RColorBrewer)
library(igraph)
library(Rgraphviz)
library(graph)
library(colorspace)
library(annotate)
library(scales)
library(gtools)

library(MineICA)

# unused exploratory analysis prepare data ---------------------------------------------------------------

# load dds
dds <- readRDS("output/dds_qc.rds")

# set to TRUE for debugging (much faster)
vsd <- assay(vst(dds, blind = TRUE)) # yields a DESeqTransform object and extracts counts matrix

# get metadata
metadata <- as.data.frame(colData(dds))

# get normalized counts
count_matrix <- as.data.frame(counts(dds, normalized = TRUE))

# define V1 subset
V1_dds <- subset_dds(dds, condition = "VISIT == 'V1'")

# use the V1 metadata of all subjects as subject metadata (no double counting of individuals)
V1_metadata <- as.data.frame(colData(V1_dds))

V1_count_matrix <- as.data.frame(counts(V1_dds, normalized = TRUE))

# apply variance-stabilizing transformation
# we will not perform this blindly, because we expect our experimental design will lead to large differences in counts
# this makes best use of our signal for downstream analysis and visualization
# however, for quality control (QC) we wanted a completely unbiased look and will used blind transformation (the default)

# set to TRUE for debugging (much faster)
V1_vsd <- as.data.frame(assay(vst(V1_dds, blind = TRUE)))

# load l2fc data
l2fc_data_imputed <- readRDS("output/l2fc_data_imputed_qc.rds")

# remove V1 samples from l2fc count data for this script
l2fc_data_imputed <- subset_l2fc_data(l2fc_data_imputed, condition = "VISIT != 'V1'")

# load euclidean distances for dispersion study
euclid_dists <- readRDS("output/euclid_dists_qc.rds")

V1_euclid_dists <- subset_euclid_dists(euclid_dists, condition = "VISIT == 'V1'")



# visualize metadata variable distributions --------------------------------------------

# visualize study as a whole

# 1. prepare the data for the treemap
# We need to count the occurrences of each DIAGNOSIS-VISIT combination
treemap_data <- V1_metadata %>%
  group_by(DIAGNOSIS, VISIT) %>%
  dplyr::summarise(Count = n(), .groups = 'drop') %>%
  ungroup() # Ensure data is ungrouped for treemap function

# calculate total samples
total_samples <- sum(treemap_data$Count)

# calculate samples per main category (DIAGNOSIS)
diagnosis_counts <- V1_metadata %>%
  group_by(DIAGNOSIS) %>%
  dplyr::summarise(DiagnosisCount = n(), .groups = 'drop')

# merge diagnosis counts back into treemap_data for labeling
treemap_data <- treemap_data %>%
  left_join(diagnosis_counts, by = "DIAGNOSIS")

# 2. create a custom label for DIAGNOSIS
# We will create a new column that combines DIAGNOSIS with its total count
treemap_data <- treemap_data %>%
  mutate(
    DIAGNOSIS_LABEL = paste0(DIAGNOSIS, "\n(N=", DiagnosisCount, ")"),
    # turn it into a factor in the exact order of your colour key:
    DIAGNOSIS = factor(DIAGNOSIS, levels = names(DIAGNOSIS_COLKEY))
  )

# 3. create the treemap
png("figs/exploration/COHORT_treemap.png", width = 3600, height = 2800, res = 300) # adjust width/height as needed
treemap(treemap_data,
        index = c("DIAGNOSIS_LABEL"), # defines the hierarchy with custom label
        vSize = "Count",                 # variable determining the size of the rectangles
        vColor = "DIAGNOSIS",            # variable determining the color of the rectangles (still based on original DIAGNOSIS for consistent coloring)
        palette = unname(DIAGNOSIS_COLKEY),       # specify how colors are mapped
        type = "categorical",                  # type of coloring (colors by the first index variable, which is now DIAGNOSIS_LABEL)
        title = paste0("Sample Distribution by Diagnosis and Visit (Total N=", total_samples, ")"), # Title with total samples
        fontsize.title = 18,
        fontsize.labels = c(8, 2),     # font sizes for main (DIAGNOSIS_LABEL) and sub-labels (VISIT)
        fontcolor.labels = c("black", "black"), # font colors for main and sub-labels
        align.labels = list(c("left", "top"), c("right", "bottom")), # alignment of labels
        overlap.labels = 0.5,            # how much labels can overlap
        inflate.labels = TRUE,           # inflate labels to fill space
        border.col = c("white", "white"), # border colors for main and sub-rectangles
        border.lwds = c(2, 1),           # border line widths
        force.print.labels = TRUE,
        fontface.labels = "bold",
)
dev.off()

# do not plot all factors with as many levels as samples (unique for every row)
# do not plot factors with uniformative constants
vars_to_plot <- names(V1_metadata)[
  !sapply(V1_metadata, function(col) {
    n_unique <- length(unique(col))
    n_unique == 1 || n_unique == nrow(V1_metadata)
  })
]

# loop through each column in the metadata dataframe and create bar or hist
for (var in vars_to_plot) {
  
  # skip certain variables
  if (var %in% c("VISIT", "SUBJECT_GROUP_NUMBER", "VISIT_DATE", "YEAR_OF_BIRTH", "IgG", "IgG_CENTERED_SCALED", "CANCER_THERAPY_AT_V1", "HSCT_DATE")) next
  
  # grab the column vector
  vec <- V1_metadata[[var]]
  
  # set up common aesthetics
  aes_mapping <- aes(x = .data[[var]])
  base_plot <- ggplot(V1_metadata, mapping = aes_mapping) +
    labs(title = paste0(if (is.numeric(vec)) "Histogram of " else "Bar chart of ", var),
         x = var,
         y = "Count") +
    theme_minimal()
  
  # choose geom based on type
  if (is.numeric(vec)) {
    # numeric -> histogram
    # choose a sensible binwidth (e.g. Sturges' rule)
    binwidth <- diff(range(vec, na.rm = TRUE)) / (6.0 * nclass.Sturges(vec))
    p <- base_plot +
      geom_histogram(binwidth = binwidth, fill = "steelblue", color = "black")
    
  } else if (is.factor(vec) || is.character(vec)) {
    # categorical -> bar chart
    p <- base_plot +
      geom_bar(fill = "steelblue", color = "black")
    
  } else {
    # skip other types (e.g. dates)
    next
  }
  
  # save the plot
  fname <- file.path("figs/exploration/metadata_distributions", paste0("COHORT_", var, if (is.numeric(vec)) "_histogram.png" else "_bar.png"))
  ggsave(filename = fname, plot = p, width = 6, height = 4)
}




# COHORT exploration ------------------------------------------------------
##### normalized gene count PCA #####

p <- pca(vsd, metadata = colData(dds), removeVar = 0.1)

# determine elbow point and plot in scree
elbow <- findElbowPoint(p$variance)

png(filename = "figs/exploration/pca/COHORT/normalized_counts/screeplot.png", width = 2400, height = 2400, res = 200)
screeplot(p,
          getComponents(p, 1:(elbow + 4)),
          title = "COHORT PCA Scree Plot in normalized counts space",
          axisLabSize = 12,
          titleLabSize = 22,
          xlabAngle = 90,
          vline = elbow) +
  annotate("text", x = elbow, y = 50,
           label = 'Elbow method', vjust = -5, hjust = 1.2, size = 6)
dev.off()

# perform eigencorrelation of metadata with PCs to identify most interesting PCA biplots and how to color them
png("figs/exploration/pca/COHORT/normalized_counts/eigencor.png", width = 1600, height = 1000, res = 150)
eigencorplot(p,
             components = getComponents(p, 1:13),
             metavars = c("DIAGNOSIS",
                          "SEX",
                          "AGE_BELOW_60",
                          "VISIT",
                          "INFLUENZA_VACCINATIONS_LAST_5_YEARS",
                          "COVID19_VACCINATIONS",
                          "HISTORY_OF_CANCER",
                          "LENALIDOMIDE_AT_V1",
                          "CHEMOTHERAPY_AT_V1"),
             col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
             cexCorval = 0.7,
             colCorval = 'white',
             fontCorval = 2,
             posLab = 'bottomleft',
             rotLabX = 45,
             posColKey = 'top',
             cexLabColKey = 1.5,
             scale = TRUE,
             main = 'COHORT PC1-13 metadata correlations in normalized counts space',
             corFUN = "pearson",
             corMultipleTestCorrection = "BH",
             colFrame = 'white',
             plotRsquared = FALSE,
             cexLabY = 0.8)
dev.off()

# create PCA biplots
biplot_configs <- list(
  list(x = "PC1", y = "PC2", colby = "DIAGNOSIS"),
  list(x = "PC1", y = "PC2", colby = "CHEMOTHERAPY_AT_V1"),
  list(x = "PC1", y = "PC2", colby = "LENALIDOMIDE_AT_V1"),
  list(x = "PC1", y = "PC2", colby = "HISTORY_OF_CANCER"),
  list(x = "PC1", y = "PC2", colby = "AGE_BELOW_60"),
  
  list(x = "PC3", y = "PC4", colby = "LENALIDOMIDE_AT_V1"),
  list(x = "PC3", y = "PC4", colby = "CHEMOTHERAPY_AT_V1"),
  list(x = "PC3", y = "PC4", colby = "COVID19_VACCINATIONS"),
  list(x = "PC3", y = "PC4", colby = "INFLUENZA_VACCINATIONS_LAST_5_YEARS"),
  list(x = "PC3", y = "PC4", colby = "SEX"),
  list(x = "PC3", y = "PC4", colby = "DIAGNOSIS"),
  list(x = "PC3", y = "PC4", colby = "HISTORY_OF_CANCER"),
  
  list(x = "PC5", y = "PC6", colby = "CHEMOTHERAPY_AT_V1"),
  list(x = "PC5", y = "PC6", colby = "LENALIDOMIDE_AT_V1"),
  list(x = "PC5", y = "PC6", colby = "SEX"),
  list(x = "PC5", y = "PC6", colby = "DIAGNOSIS"),
  list(x = "PC5", y = "PC6", colby = "HISTORY_OF_CANCER"),
  
  list(x = "PC7", y = "PC8", colby = "CHEMOTHERAPY_AT_V1"),
  list(x = "PC7", y = "PC8", colby = "LENALIDOMIDE_AT_V1"),
  
  list(x = "PC9", y = "PC10", colby = "CHEMOTHERAPY_AT_V1"),
  list(x = "PC9", y = "PC10", colby = "LENALIDOMIDE_AT_V1"),
  list(x = "PC9", y = "PC10", colby = "DIAGNOSIS"),
  
  list(x = "PC11", y = "PC13", colby = "LENALIDOMIDE_AT_V1"),
  list(x = "PC11", y = "PC13", colby = "AGE_BELOW_60"),
  list(x = "PC11", y = "PC13", colby = "DIAGNOSIS")
)

for (conf in biplot_configs) {
  
  # set custom colkey if applicable
  custom_colkey <- NULL
  if (conf$colby == "DIAGNOSIS") {
    custom_colkey <- DIAGNOSIS_COLKEY
  } else if (conf$colby == "VISIT") {
    custom_colkey <- VISIT_COLKEY
  }
  
  png(filename = paste0("figs/exploration/pca/COHORT/normalized_counts/biplot_", conf$x, "_", conf$y, "_colby_", conf$colby, ".png"), width = 2600, height = 2400, res = 200)
  p_biplot <- biplot(p,
                     title = "COHORT V1 PCA Biplot in normalized counts space",
                     x = conf$x,
                     y = conf$y,
                     showLoadings = TRUE,
                     lab = NULL,
                     selectLab = NULL,
                     legendIconSize = 5,
                     labSize = 5,
                     pointSize = 5,
                     sizeLoadingsNames = 5,
                     colby = conf$colby,
                     colLegendTitle = conf$colby,
                     colkey = custom_colkey,
                     legendPosition = "top")
  print(p_biplot)  # explicitly render the plot
  dev.off()
}


##### normalized gene count ICA #####

# set the number of components (adjust based on data complexity)
n_comp <- 18  #  experiment with different values

# perform ICA
ica_result <- runICA(X = as.data.frame(vsd), 
                     nbComp = n_comp, 
                     method = "JADE", # alternative: "FastICA" and clusterFastICARuns to get a consensus
                     maxit = 10000)

# with our ICA result ready, we create the mineICAParams object to include it in our final IcaSet instance
params <- buildMineICAParams(resPath = "figs/exploration/ica/COHORT/normalized_counts/", selCutoff = 3, pvalCutoff = 0.05, annotfile = "NULL")

ica_set_build <- buildIcaSet(params=params, A=data.frame(ica_result$A), S = data.frame(ica_result$S), dat = as.data.frame(vsd),
                             pData = metadata, alreadyAnnot = TRUE)


ica_set <- ica_set_build$icaSet
params <- ica_set_build$params


# restrict the metadata to the variables of interest
keepVar <- c("DIAGNOSIS",
             "SEX",
             "AGE_BELOW_60",
             "VISIT",
             "INFLUENZA_VACCINATIONS_LAST_5_YEARS",
             "COVID19_VACCINATIONS")

# run the analysis of the ICA decomposition
try(runAn(params=params, icaSet=ica_set, writeGenesByComp = TRUE,
          keepVar=keepVar))

# this way we can find out which samples contribute to a specific IC (which samples are responsible for that signal)

# analysis function can also be called individually
# ## heatmap with dendrogram
# resH <- plot_heatmapsOnSel(icaSet = ica_set, selCutoff = 5, level = "genes",
#                            keepVar = keepVar,
#                            doSamplesDendro = TRUE, doGenesDendro = TRUE, keepComp = 2,
#                            file = "figs/exploration/ica/TEST/heatmap.pdf", annot2col=annot2col(params))


# let us try to understand what each independent component could represent by enriching with the genes that have the strongest projections onto ICs

# extract the gene projections (loadings) onto ICs, only take those with an absolute projection of 5 or higher
contrib <- selectContrib(ica_set, cutoff = 5, level="genes")

# loop over ICs and perform enrichment for BTMs
for (comp in 1:length(contrib)) {
  comp_genes <- names(sort(abs(contrib[[comp]]), decreasing=TRUE))
  ego <- dotplot_enrichment(comp_genes)
  
  if (!is.null(ego$plot)) {
    enriched_plot <- ego$plot + ggtitle(paste0("COHORT - BTM enrichment for genes projecting onto IC", comp, " in normalized count space"))
    ggsave(
      filename = paste0("figs/exploration/ica/COHORT/normalized_counts/enrichment/IC", comp, "_enrichment_dotplot.png"),
      plot = enriched_plot,
      height = 8,
      width = 10
    )
  } else {
    message(paste("No enrichment found for component", comp))
  }
}


##### normalized genecount boxplots of genes and module scores #####

### Eigengene calculation ###

# example pathway
pathway <- "plasma cells & B cells, immunoglobulins (M156.0)"

## negative control pathway with random genes
# # get genes associated with that pathway
# pathway_genes <- PATHWAYS[[pathway]]
# 
# # generate negative pathway (random gene set)
# # how many genes in the pathway?
# n_genes <- length(pathway_genes)
# 
# # all genes available in your count matrix
# all_genes <- rownames(vsd_subset)
# 
# # exclude the pathway genes if you want no overlap
# background_genes <- setdiff(all_genes, pathway_genes)
# 
# # sample random genes as a negative control set
# set.seed(24)
# random_genes <- sample(background_genes, n_genes, replace = FALSE)
# 
# # assign to a new pathway name for plotting, etc.
# negative_control_pathway <- "random control (same size as M156.0)"
# PATHWAYS[[negative_control_pathway]] <- random_genes

# get genes associated with that pathway
pathway_genes <- PATHWAYS[[pathway]]

# subset count matrix to genes of pathway
# remove those genes defined in the pathway but not present in our count_matrix
vsd_pathway_subset <- drop_na(V1_vsd[pathway_genes, ])

# run PCA on the module
# removeVar=0 means keep all genes (no variance filtering)
p <- pca(
  mat       = as.matrix(vsd_pathway_subset),
  metadata  = V1_metadata,
  removeVar = 0
)

# extract the first PC (the "eigengene")
# p$rotated is a samples×PCs matrix
eigengene_df <- data.frame(
  sample    = rownames(p$rotated),
  PC1 = p$rotated[, 1],
  stringsAsFactors = FALSE
)

# merge with metadata
# assume rownames(metadata) are sample names
eigengene_df <- merge(
  eigengene_df,
  data.frame(sample = rownames(V1_metadata), V1_metadata, stringsAsFactors = FALSE),
  by = "sample"
)

# boxplot of module eigengene by diagnosis
ggplot(eigengene_df, aes(x = DIAGNOSIS, y = PC1, colour = DIAGNOSIS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(
    title = paste0("Module eigengene for ‘", pathway, "’"),
    x     = "Diagnosis",
    y     = "PC1 (eigengene) score"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.text.x  = element_text(angle = 45, hjust = 1)
  ) + scale_colour_manual(values = DIAGNOSIS_COLKEY)

### GSVA boxplot ###

# extract the pathway's gsva score
pathway_df <- data.frame(
  CEGAT_ID    = rownames(V1_pathway_data$metadata),
  pathway_gsva_score = V1_pathway_data$gsva_matrix[pathway, ],
  stringsAsFactors = FALSE
)

# merge with metadata
pathway_df <- merge(
  pathway_df,
  data.frame(V1_pathway_data$metadata, stringsAsFactors = FALSE),
  by = "CEGAT_ID"
)

# boxplot of module gsva score by diagnosis
ggplot(pathway_df, aes(x = DIAGNOSIS, y = pathway_gsva_score, colour = DIAGNOSIS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(
    title = paste0("Module GSVA score for ‘", pathway, "’"),
    x     = "Diagnosis",
    y     = "GSVA score"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.text.x  = element_text(angle = 45, hjust = 1)
  ) +
  scale_colour_manual(values = DIAGNOSIS_COLKEY) +
  stat_compare_means(method = "kruskal.test", label.y = max(pathway_df$pathway_gsva_score) * 1.05)



##### l2fc gene count variance partitioning analysis #####

# Let us investigate how well the metadata variables we assume to explain a lot of differences in our data actually perform

# Now in log2Fold change space we believe VISIT to explain much more variance, while subject baseline differences should be reduced (or removed alltogether)

# let us model STUDY_ID + VISIT
form <- ~ (1 | STUDY_ID) + (1 | VISIT)

# the varPart model enforces us to model either all or none of the categorical variables as random effects
# we will model both STUDY_ID and VISIT as fixed effects and not assume that STUDY_ID comes from a larger population (i.e. has probability distribution)

# takes around 5 minutes
varPart <- fitExtractVarPartModel(l2fc_data_imputed$count_matrix, form, l2fc_data_imputed$metadata)

# create violin plot of contribution of each variable to total variance
plotVarPart(varPart)

# save the plot
ggsave("figs/exploration/variance_partitioning_analysis/COHORT/l2fc_counts_STUDY_ID_VISIT_variance_violins.png")

for (var in c("STUDY_ID", "VISIT", "Residuals")) {
  
  # order varPart with descending variance explained by variable
  ordered_varPart <- varPart[order(varPart[[var]], decreasing = TRUE), ]
  
  # define number of genes to enrich
  num_of_top_genes = 200
  
  # select the top num_of_top_genes rows (genes)
  top_genes <- head(ordered_varPart, num_of_top_genes)
  
  # get num_of_top_genes genes in which variable explains the most variance
  top_gene_names <- rownames(top_genes)
  
  ego <- dotplot_enrichment(top_gene_names)
  
  if (!is.null(ego$plot)) {
    enriched_plot <- ego$plot + ggtitle(paste0("COHORT - BTM enrichment for top genes with most variance explained by ", var, " in l2fc counts space"))
    ggsave(
      filename = paste0("figs/exploration/variance_partitioning_analysis/COHORT/l2fc_counts_variance_violins_", var, "_top_genes_enrichment.png"),
      plot = enriched_plot,
      height = 8,
      width = 10
    )
  } else {
    message(paste("No enrichment found for var", var))
  }
}

##### l2fc gene count PCA #####

p <- pca(l2fc_count_matrix_scaled, metadata = l2fc_data_imputed$metadata, removeVar = 0.1)

# determine elbow point and plot in scree
elbow <- findElbowPoint(p$variance)

png(filename = "figs/exploration/pca/COHORT/l2fc_counts/screeplot.png", width = 2400, height = 2400, res = 200)
screeplot(p,
          getComponents(p, 1:(elbow + 4)),
          title = "COHORT PCA Scree Plot in l2fc counts space",
          axisLabSize = 12,
          titleLabSize = 22,
          xlabAngle = 90,
          vline = elbow) +
  annotate("text", x = elbow, y = 50,
           label = 'Elbow method', vjust = -5, hjust = 1.2, size = 6)
dev.off()

# perform eigencorrelation of metadata with PCs to identify most interesting PCA biplots and how to color them
png("figs/exploration/pca/COHORT/l2fc_counts/eigencor.png", width = 1600, height = 1000, res = 150)
eigencorplot(p,
             components = getComponents(p, 1:4),
             metavars = c("DIAGNOSIS",
                          "SEX",
                          "AGE_BELOW_60",
                          "VISIT",
                          "INFLUENZA_VACCINATIONS_LAST_5_YEARS",
                          "COVID19_VACCINATIONS",
                          "HISTORY_OF_CANCER",
                          "LENALIDOMIDE_AT_V1",
                          "CHEMOTHERAPY_AT_V1"),
             col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
             cexCorval = 0.7,
             colCorval = 'white',
             fontCorval = 2,
             posLab = 'bottomleft',
             rotLabX = 45,
             posColKey = 'top',
             cexLabColKey = 1.5,
             scale = TRUE,
             main = 'COHORT PC1-4 metadata correlations in l2fc counts space',
             corFUN = "pearson",
             corMultipleTestCorrection = "BH",
             colFrame = 'white',
             plotRsquared = FALSE,
             cexLabY = 0.8)
dev.off()

# create PCA biplots
biplot_configs <- list(
  list(x = "PC1", y = "PC2", colby = "VISIT"),
  list(x = "PC1", y = "PC2", colby = "CHEMOTHERAPY_AT_V1"),
  
  list(x = "PC3", y = "PC4", colby = "LENALIDOMIDE_AT_V1"),
  list(x = "PC3", y = "PC4", colby = "SEX"),
  list(x = "PC3", y = "PC4", colby = "VISIT"),
  list(x = "PC3", y = "PC4", colby = "AGE_BELOW_60"),
  list(x = "PC3", y = "PC4", colby = "VISIT"),
  list(x = "PC3", y = "PC4", colby = "HISTORY_OF_CANCER")
)

for (conf in biplot_configs) {
  png(filename = paste0("figs/exploration/pca/COHORT/l2fc_counts/biplot_", conf$x, "_", conf$y, "_colby_", conf$colby, ".png"), width = 2600, height = 2400, res = 200)
  p_biplot <- biplot(p,
                     title = "COHORT V1 PCA Biplot in l2fc counts space",
                     x = conf$x,
                     y = conf$y,
                     showLoadings = TRUE,
                     lab = NULL,
                     selectLab = NULL,
                     legendIconSize = 5,
                     labSize = 5,
                     pointSize = 5,
                     sizeLoadingsNames = 5,
                     colby = conf$colby,
                     colLegendTitle = conf$colby,
                     legendPosition = "top")
  print(p_biplot)  # explicitly render the plot
  dev.off()
}

##### l2fc gene count ICA #####

# let's keep just visits V2, V3 and V4 for now as these are probably the most interesting
#l2fc_data_imputed_subset <- subset_l2fc_data(l2fc_data_imputed, condition = "VISIT %in% c('V2', 'V3', 'V4')")

n_comp <- 10  #  experiment with different values


# perform ICA
ica_result <- runICA(X = l2fc_data_imputed$count_matrix, 
                     nbComp = n_comp, 
                     method = "JADE", # alternative: "FastICA" and clusterFastICARuns to get a consensus
                     maxit = 10000)

# with our ICA result ready, we create the mineICAParams object to include it in our final IcaSet instance
params <- buildMineICAParams(resPath = "figs/exploration/ica/COHORT/l2fc_counts/", selCutoff = 3, pvalCutoff = 0.05, annotfile = "NULL")

ica_set_build <- buildIcaSet(params=params, A=data.frame(ica_result$A), S = data.frame(ica_result$S), dat = as.data.frame(l2fc_data_imputed$count_matrix),
                             pData = l2fc_data_imputed$metadata, alreadyAnnot = TRUE)


ica_set <- ica_set_build$icaSet
params <- ica_set_build$params


# restrict the metadata to the variables of interest
keepVar <- c("DIAGNOSIS",
             "SEX",
             "AGE_BELOW_60",
             "VISIT",
             "INFLUENZA_VACCINATIONS_LAST_5_YEARS",
             "COVID19_VACCINATIONS")

# run the analysis of the ICA decomposition
try(runAn(params=params, icaSet=ica_set, writeGenesByComp = TRUE,
          keepVar=keepVar))
# this way we can find out which samples contribute to a specific IC (which samples are responsible for that signal)

# let us try to understand what each independent component could represent by enriching with the genes that have the strongest projections onto ICs
# extract the gene projections (loadings) onto ICs, only take those with an absolute projection of 3 or higher
contrib <- selectContrib(ica_set, cutoff = 3, level="genes")

# loop over ICs and perform enrichment for BTMs
for (comp in 1:length(contrib)) {
  comp_genes <- names(sort(abs(contrib[[comp]]), decreasing=TRUE))
  ego <- dotplot_enrichment(comp_genes)
  
  if (!is.null(ego$plot)) {
    enriched_plot <- ego$plot + ggtitle(paste0("COHORT - BTM enrichment for genes projecting onto IC", comp, " in l2fc count space"))
    ggsave(
      filename = paste0("figs/exploration/ica/COHORT/l2fc_counts/enrichment/IC", comp, "_enrichment_dotplot.png"),
      plot = enriched_plot,
      height = 8,
      width = 10
    )
  } else {
    message(paste("No enrichment found for component", comp))
  }
}

##### l2fc gene count boxplots of genes and module scores #####

# subset data to particular visit to show kinetics
l2fc_data_imputed_subset <- subset_l2fc_data(l2fc_data_imputed, condition = "VISIT == 'V4'")

count_matrix_subset <- as.data.frame(l2fc_data_imputed_subset$count_matrix)

metadata_subset <- l2fc_data_imputed_subset$metadata

# example pathway
pathway <- "enriched in NK cells (I) (M7.2)"

# get genes associated with that pathway
pathway_genes <- PATHWAYS[[pathway]]

# subset count matrix to genes of pathway
# remove those genes defined in the pathway but not present in our count_matrix
count_matrix_pathway_subset <- drop_na(count_matrix_subset[pathway_genes, ])

# run PCA on the module
# removeVar=0 means keep all genes (no variance filtering)
p <- pca(
  mat       = as.matrix(count_matrix_pathway_subset),
  metadata  = metadata_subset,
  removeVar = 0
)

# extract the first PC (the "eigengene")
# p$rotated is a samples×PCs matrix
eigengene_df <- data.frame(
  sample    = rownames(p$rotated),
  PC1 = p$rotated[, 1],
  stringsAsFactors = FALSE
)

# merge with metadata (which has a column DIAGNOSIS)
# assume rownames(metadata) are sample names
eigengene_df <- merge(
  eigengene_df,
  data.frame(sample = rownames(metadata_subset), metadata_subset, stringsAsFactors = FALSE),
  by = "sample"
)

# boxplot of module eigengene by diagnosis
ggplot(eigengene_df, aes(x = DIAGNOSIS, y = PC1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(
    title = paste0("Module eigengene for ‘", pathway, "’"),
    x     = "Diagnosis",
    y     = "PC1 (eigengene) score"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.text.x  = element_text(angle = 45, hjust = 1)
  ) + scale_fill_manual(values = DIAGNOSIS_COLKEY)








#### Differential Gene Expression Analysis ####

# unused differential gene expression analysis libraries ------------------
# core tidyverse packages
library(dplyr)       # Data manipulation (before tidyr, purrr, and tibble)
library(tidyr)       # Data tidying
library(purrr)       # Functional programming (works with dplyr)
library(tibble)      # Tidy data frames (works with dplyr)

# differential expression analysis
library(BiocParallel)
library(DESeq2)      # RNA-seq analysis

# ANOVA model testing
library(AICcmodavg)

# gene set enrichment analysis
library(fgsea)       # GSEA

# visualization
library(RColorBrewer) # Color palettes
library(reshape2)
library(ggplot2)     # Core plotting package
library(ggh4x)       # for fixed pnel sizes
library(ggpubr)      # Builds on ggplot2 (should come after it)
library(UpSetR)      # Upset plot
library(cowplot)     # for the fancy 2d grid in the degs bar plot
library(PCAtools)
library(colorspace)  # for changing color saturation in code

# heatmaps and clustering
library(ComplexHeatmap) # Detailed heatmaps
library(circlize)       # Color scaling for heatmaps
library(grid)           # for gpar(), grid.text()
library(seriation)      # improved row ordering in heatmaps
library(dendextend)     # dendrogram customization
library(magick)         # image processing (for ComplexHeatmap)

# enhanced volcano plot
library(EnhancedVolcano)
library(patchwork) # create plots side by side
# unused differential gene expression analysis functions ------------------

# function to create volcano plot for specific DIAGNOSIS at specific VISIT
plot_degs_volcano <- function(dgea_result,
                              x = "log2FoldChange",
                              y = "padj",
                              p_cut = 0.05,
                              l2fc_cut = 0.2,
                              xlim = c(-2, 2)) {
  
  # set y-axis label
  if (y == "padj") {
    ylab <- bquote(-Log[10](italic(p)[adj]))
  } else if (y == "pvalue") {
    ylab <- bquote(-Log[10](italic(p)))
  }
  
  # create a new column that clamps the x-values to the boundaries defined in xlim
  dgea_result$plotX <- dgea_result[[x]]
  dgea_result$plotX[dgea_result$plotX > xlim[2]] <- xlim[2]
  dgea_result$plotX[dgea_result$plotX < xlim[1]] <- xlim[1]
  
  # call EnhancedVolcano using the capped values.
  EnhancedVolcano(dgea_result,
                  lab = dgea_result$row,
                  x = "plotX",
                  xlim = xlim,
                  xlab = bquote(~Log[2] ~ "fold change"),
                  y = y,
                  ylab = ylab,
                  pCutoff = p_cut,
                  FCcutoff = l2fc_cut,
                  subtitle = NULL)
}

# visualize DEGs with upset plot
plot_significant_degs_upset <- function(dgea_results, visit) {
  
  # Convert HC genes to character to ensure proper indexing
  #hc_genes <- as.character(dgea_results[["HC"]][[visit]][["significant_results"]][["row"]])
  #cat("Number of significant genes for HC:", length(unique(hc_genes)), "\n")
  
  # Initialize an empty list to store significant genes for each diagnosis
  significant_genes_list <- list()
  
  # Loop through all diagnoses and collect their significant genes
  for (diagnosis in names(dgea_results)) {
    # Convert gene names to character to avoid factor indexing issues
    current_genes <- as.character(dgea_results[[diagnosis]][[visit]][["significant_results"]][["row"]])
    significant_genes_list[[diagnosis]] <- current_genes
  }
  
  # Create a binary matrix for UpSetR based on the union of all genes
  all_genes <- unique(unlist(significant_genes_list))
  binary_matrix <- matrix(0, nrow = length(all_genes), ncol = length(significant_genes_list))
  rownames(binary_matrix) <- all_genes
  colnames(binary_matrix) <- names(significant_genes_list)
  
  # Fill the binary matrix with 1s where the gene is present for the diagnosis
  for (i in seq_along(significant_genes_list)) {
    binary_matrix[significant_genes_list[[i]], i] <- 1
  }
  
  # Extract genes for unique & combination subsets
  subset_genes <- list()
  for (i in 1:nrow(binary_matrix)) {
    active_sets <- names(which(binary_matrix[i, ] == 1))
    if (length(active_sets) > 0) {
      subset_name <- paste(sort(active_sets), collapse = "_")  # create a name like "HC_LC"
      subset_genes[[subset_name]] <- c(subset_genes[[subset_name]], rownames(binary_matrix)[i])
    }
  }
  
  
  # Create the UpSet plot with colors
  p <- upset(as.data.frame(binary_matrix),
             sets = c("SMM", "MM", "NSCLC", "HC"),  # force this specific order
             keep.order = TRUE,                  # preserve the given order
             mainbar.y.label = "Intersection Size", 
             sets.x.label = "Set Size",
             text.scale = 1.2,
             sets.bar.color = DIAGNOSIS_COLKEY[c("SMM", "MM", "NSCLC", "HC")], 
             matrix.color = "black")
  
  return(list(plot = p, subset_genes = subset_genes))
}

# to plot a log2FoldChange histogram for group "V3" and diagnosis "MM":
#plot_dgea_variable_distribution(dgea_results, group = "V3", condition = "NSCLC", variable_to_plot = "l2fc")

plot_l2fc_distribution_grid <- function(dgea_results, groups, condition, plot_range = c(-5, 5),
                                        hist_breaks = 500, thresh_vals = c(-0.2, 0.2)) {
  # Determine grid dimensions automatically
  n_plots <- length(groups)
  n_row <- ceiling(sqrt(n_plots))
  n_col <- ceiling(n_plots / n_row)
  
  # Save current graphic parameters and set new layout
  op <- par(mfrow = c(n_row, n_col))
  
  for (grp in groups) {
    if (!grp %in% names(dgea_results)) {
      warning(paste("Group", grp, "not found in dgea_results"))
      next
    }
    
    # Extract the results data frame for this group and condition
    group_data <- dgea_results[[grp]][[condition]][["results"]]
    
    if (is.null(group_data) || !("log2FoldChange" %in% colnames(group_data))) {
      warning(paste("No 'log2FoldChange' column found for", grp, "at", condition))
      next
    }
    
    l2fc_vals <- group_data$log2FoldChange
    
    # Identify any extreme values beyond the plotting range
    extreme_vals <- l2fc_vals[l2fc_vals < plot_range[1] | l2fc_vals > plot_range[2]]
    
    # Plot the histogram using density rather than frequency
    hist(l2fc_vals,
         breaks = hist_breaks,
         col = "skyblue",
         border = "white",
         main = paste("DESeq2 log2FoldChange for", grp, "at", condition),
         xlab = "log2FoldChange",
         xlim = plot_range,
         ylab = "Density",
         freq = FALSE)
    
    # Add vertical lines for the threshold values
    abline(v = thresh_vals, col = "red", lty = 2, lwd = 2)
    
    # If there are extreme values, annotate with a legend
    if (length(extreme_vals) > 0) {
      legend_text <- paste(length(extreme_vals), "outlier(s) beyond", paste(plot_range, collapse = " to "))
      legend("topright", legend = legend_text, bty = "n", cex = 0.8)
    }
  }
  
  # Restore original graphic parameters
  par(op)
}

# example:
#plot_l2fc_distribution_grid(dgea_results, groups = DIAGNOSES, condition = "V2")


plot_padj_distribution_grid <- function(dgea_results, groups, condition, x_limits = c(0, 0.1),
                                        hist_breaks = 500, padj_threshold = 0.05) {
  # Determine grid dimensions automatically
  n_plots <- length(groups)
  n_row <- ceiling(sqrt(n_plots))
  n_col <- ceiling(n_plots / n_row)
  
  # Save current graphic settings and set new layout
  op <- par(mfrow = c(n_row, n_col))
  
  for (grp in groups) {
    if (!grp %in% names(dgea_results)) {
      warning(paste("Group", grp, "not found in dgea_results"))
      next
    }
    
    # Extract the results data frame for this group and condition
    group_data <- dgea_results[[grp]][[condition]][["results"]]
    
    if (is.null(group_data) || !("padj" %in% colnames(group_data))) {
      warning(paste("No 'padj' column found for", grp, "at", condition))
      next
    }
    
    padj_vals <- group_data$padj
    
    # Plot the histogram for adjusted p-values
    hist(padj_vals,
         breaks = hist_breaks,
         col = "skyblue",
         border = "white",
         main = paste("DESeq2 padj distribution for", grp, "at", condition),
         xlab = "padj",
         xlim = x_limits,
         ylab = "Frequency")
    
    # Add a vertical line for the p-value threshold
    abline(v = padj_threshold, col = "red", lty = 2, lwd = 2)
  }
  
  # Restore original graphic parameters
  par(op)
}

# example
#plot_padj_distribution_grid(dgea_results, groups = VISITS, condition = "NSCLC")

# function to plot gene-level violins for dgea_results variable
plot_variable_violins <- function(dgea_results,
                                  variable_to_plot = "pval",  # options: "pval", "padj", "log2fc"
                                  groups,                    # e.g., c("HC", "NSCLC", "MM", "SMM")
                                  conditions,                # e.g., c("V1", "V2", "V3")
                                  group_label = "Group",     # label for the primary grouping (x-axis)
                                  condition_label = "Condition") {  # label for the secondary grouping (facets)
  
  # Initialize an empty list to store data frames
  results_list <- list()
  
  # Loop over each combination of group and condition
  for (grp in groups) {
    for (cond in conditions) {
      # Check if the nested result exists for this combination
      if (!is.null(dgea_results[[grp]][[cond]][["results"]])) {
        # Choose the appropriate column from the results
        col_to_use <- switch(variable_to_plot,
                             "pval"   = "pvalue",
                             "padj"   = "padj",
                             "log2fc" = "log2FoldChange",
                             stop("Invalid 'variable_to_plot' argument. Choose 'pval', 'padj', or 'log2fc'."))
        
        # Subset the data: assume each result data frame has a column "row" for gene names
        df <- dgea_results[[grp]][[cond]][["results"]][, c("row", col_to_use)]
        
        # Rename the chosen column to "y_value" for ease of plotting
        colnames(df)[2] <- "y_value"
        
        # Add columns to record the group and condition for each gene
        df$Group <- grp
        df$Condition <- cond
        
        # Store the subset data frame in the results list (using a combined name)
        results_list[[paste(grp, cond, sep = "_")]] <- df
      }
    }
  }
  
  # Combine all subset data frames into one data frame
  plot_data <- do.call(rbind, results_list)
  
  # Convert the plotting variable to numeric (if not already)
  plot_data$y_value <- as.numeric(plot_data$y_value)
  
  # Set labels and scaling parameters based on the variable of interest
  if (variable_to_plot == "pval") {
    y_label <- "p-value"
    file_suffix <- "pvalue"
    y_scale <- "log10"  # apply log scale for p-values
  } else if (variable_to_plot == "padj") {
    y_label <- "adjusted p-value"
    file_suffix <- "padj"
    y_scale <- "log10"  # apply log scale for adjusted p-values
  } else if (variable_to_plot == "log2fc") {
    y_label <- "log2 fold change"
    file_suffix <- "log2foldchange"
    y_scale <- "identity"  # no transformation for log2 fold change
  }
  
  # Create a violin plot (with an overlaid boxplot) grouped by the primary grouping variable
  # and faceted by the secondary condition.
  p <- ggplot(plot_data, aes(x = Group, y = y_value)) +
    geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    theme_minimal() +
    labs(title = paste(y_label, "Distribution per", condition_label),
         x = group_label,
         y = y_label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Apply a log10 scale if appropriate
  if (y_scale == "log10") {
    p <- p + scale_y_continuous(trans = "log10")
  }
  
  # Facet the plot by the secondary condition
  p <- p + facet_wrap(~ Condition, nrow = 1)
  
  # Optionally, save the plot to file. (Here 'study_subset' is assumed to be defined elsewhere.)
  # ggsave(paste0("figs/analysis/SYS_01_", study_subset, "_", file_suffix, "_violins.png"), plot = p)
  
  return(p)
}
# example:
# plot_variable_violins(dgea_results,
#                       variable_to_plot = "log2fc",
#                       groups = DIAGNOSES,
#                       conditions = VISITS[-1],
#                       group_label = "Diagnosis",
#                       condition_label = "Visit")
# unused differential gene expression analysis prepare data ---------------
# prepare data ---------------------------------------------------------------

# read in dds
dds <- readRDS("output/dds_qc.rds")

metadata <- as.data.frame(colData(dds))

# # set reference level for VISIT
# dds$VISIT <- relevel(dds$VISIT, ref = "V1")
# 
# # set reference level for DIAGNOSIS
# dds$DIAGNOSIS <- relevel(dds$DIAGNOSIS, ref = "HC")

# initialize list to store results of DGEA, access per DIAGNOSIS
dgea_results <- list()

# initialize list to store results of GSEA, access per DIAGNOSIS
gsea_results <- list()

# load l2fc_count_matrix from 03_exploratory_analysis
l2fc_data <- readRDS("output/l2fc_data_qc.rds")

# load l2fc_count_matrix_imputed from 03_exploratory_analysis
l2fc_data_imputed <- readRDS("output/l2fc_data_imputed_qc.rds")

# remove V1 for now
l2fc_data <- subset_l2fc_data(l2fc_data, condition = "VISIT != 'V1'")

l2fc_data_imputed <- subset_l2fc_data(l2fc_data_imputed, condition = "VISIT != 'V1'")

# load gsva pathway data
pathway_data <- readRDS("output/pathway_data.rds")

# visualize subject-wise pathway scores at V1
V1_pathway_data <- subset_pathway_data(pathway_data, condition = "VISIT == 'V1'")

# define V1 subset
V1_dds <- subset_dds(dds, condition = "VISIT == 'V1'")

# set to TRUE for debugging (much faster)
V1_vsd <- as.data.frame(assay(vst(V1_dds, blind = TRUE)))

V1_metadata <- as.data.frame(colData(V1_dds))

V1_count_matrix <- as.data.frame(counts(V1_dds, normalized = TRUE))

# ---currently not working loop to subset dds to 10 subjects from each diag
# # subset data to 10 random subjects for each diagnosis
# # get unique STUDY_IDs
# unique_study_ids <- unique(colData(dds)$STUDY_ID)

# for (diagnosis in DIAGNOSES) {
# # randomly select 10 STUDY_IDs if there are more than 10
#   if (length(unique_study_ids) > 10) {
#     selected_study_ids <- sample(unique_study_ids, 10)
#   
#     # subset dds_subset to only include the selected STUDY_IDs
#     dds <- dds[, colData(dds_subset)$STUDY_ID %in% selected_study_ids]
#     message(paste("Randomly selected 10 subjects for", diagnosis))
#   } else {
#     message(paste("Less than or equal to 10 subjects for", diagnosis, ". Using all subjects."))
#   }
# }
# ---currently not working loop to subset dds to 10 subjects from each diag


# DESeq2 analysis for cohort with interaction term ~ SEX + AGE + DIAGNOSIS + VISIT + DIAGNOSIS:VISIT ---------------------------------------------------------------

# set design
design(dds) <- ~ SEX + AGE + DIAGNOSIS + VISIT + DIAGNOSIS:VISIT

# check if DESeq2 run exists
dds_file <- "output/COHORT_DIAGNOSIS_VISIT_interaction_contrast_dds.rds"
if (file.exists(dds_file)) {
  # load existing results
  dds <- readRDS(dds_file)
  message("Loaded existing DESeq2 results from file.")
} else {
  message(paste("Running DESeq2 for Diagnosis:Visit interaction"))
  # run DESeq2
  dds <- DESeq(dds, parallel = TRUE, BPPARAM=MulticoreParam(2)) # two cores works quite well, more could have diminishing returns
  # save dds
  saveRDS(dds, file = dds_file)
}

# effect of V3 vs 1 in Healthy + the effects of V3 vs 1 in SMM (the additional coefficients)
#res <- results(dds, contrast = list("VISIT_V3_vs_V1", "DIAGNOSISMM.VISITV3"), tidy = TRUE, independentFiltering = TRUE)

# and this way we would look at just what comes on top or is missing as a VISIT effect in MM when compared to the effect it had in HC
#interaction_res <- results(dds, name = "DIAGNOSISMM.VISITV3", tidy = TRUE, independentFiltering = TRUE)

for (diagnosis in DIAGNOSES) {
  
  if (diagnosis != "HC") {
    
    for (visit in VISITS[-1]) {
      
      #res <- results(dds, contrast = list(paste0("VISIT_", visit, "_vs_V1"), paste0("DIAGNOSIS", diagnosis, ".VISIT", visit)), tidy = TRUE, independentFiltering = TRUE)
      res <- results(dds, name = paste0("DIAGNOSIS", diagnosis, ".VISIT", visit), tidy = TRUE, independentFiltering = TRUE)
      
      # remove NAs due to being extreme outliers (but not due to independent filtering)
      num_genes_before <- length(res$row)
      
      res <- res[!is.na(res$pvalue), ] # how many get removed here?
      
      num_genes_after <- length(res$row)
      
      # print how many genes were removed due to being extreme outliers
      removed_genes <- num_genes_before - num_genes_after
      message(paste("Removed", removed_genes, "genes due to being extreme outliers in", diagnosis, " ", visit))
      
      # filter most significant results
      significant_genes <- subset(res, (pvalue < 0.001) & (abs(log2FoldChange) > 0.2))
      
      # name them differently here
      dgea_results[[paste0(diagnosis, "_INTERACTION")]][[visit]][["results"]] <- res
      dgea_results[[paste0(diagnosis, "_INTERACTION")]][[visit]][["significant_results"]] <- significant_genes
      
    }
  }
}

# generate legend colors with desaturated downregulated versions
legend_colors <- unlist(lapply(names(DIAGNOSIS_COLKEY), function(diagnosis) {
  base_color <- DIAGNOSIS_COLKEY[[diagnosis]]
  # desaturate base color
  down_color <- hex(colorspace::desaturate(hex2RGB(base_color), 0.5))
  
  # create named vector with up/down colors
  c(
    setNames(down_color, paste0(diagnosis, "_INTERACTION_Downregulated")),
    setNames(base_color, paste0(diagnosis, "_INTERACTION_Upregulated"))
  )
}))


# plot differentially expressed genes by VISIT and DIAGNOSIS
deg_barplot <- plot_bar_by_group_and_condition(dgea_results,
                                               groups = paste0(DIAGNOSES[-1], "_INTERACTION"),
                                               conditions = VISITS[-1],
                                               ref_condition = "V1",
                                               group_label = "Diagnosis",
                                               condition_label = "Visit",
                                               legend_colors = legend_colors)

ggsave("figs/analysis/dgea/COHORT_degs_bar_by_DIAGNOSIS_INTERACTION_and_VISIT.png", plot = deg_barplot, width = 12, height = 8, dpi = 300)



# DESeq2 TRI_ELISA upper quartile High vs Low across all diagnoses + the cohort & visits ~ SEX + AGE + TRI_ELISA_HIGH ---------------

tri_variables = c("TRI_ELISA_V3", "TRI_ELISA_V4", "TRI_ELISA_V5")

metadata <- as.data.frame(colData(dds))

for (tri_variable in tri_variables) {
  # dynamically create the label for the new column based on the current tri_variable
  tri_high_colname <- paste0(tri_variable, "_HIGH")
  
  dgea_results_tri <- list()
  gsea_results_tri <- list()
  
  dgea_results_tri_cohort  <- list()
  gsea_results_tri_cohort  <- list()
  
  # initialize the dynamic TRI_ELISA_VX_HIGH column
  colData(dds)[[tri_high_colname]] <- NA_character_
  
  for (diag in DIAGNOSES) {
    message("Processing diagnosis: ", diag, " for TRI variable: ", tri_variable)
    
    # 1) compute threshold for tri within this diagnosis at V1 (to not count any subjects twice)
    idx_threshold_calc <- colData(dds)$DIAGNOSIS == diag & colData(dds)$VISIT == "V1"
    
    # calculate upper quartile threshold for High responders
    upper_quantile_threshold <- quantile(colData(dds)[[tri_variable]][idx_threshold_calc], probs = 0.75, na.rm = TRUE)
    message("Using upper quartile (", round(upper_quantile_threshold, 3), ") for ", tri_variable, " in ", diag)
    
    # Assign High/Low for V1 samples in this diag:
    # All values in the upper quartile get "High", otherwise "Low".
    colData(dds)[[tri_high_colname]][idx_threshold_calc] <- ifelse(
      colData(dds)[[tri_variable]][idx_threshold_calc] >= upper_quantile_threshold, "High", "Low"
    )
    
    # 3) propagate each STUDY_ID’s V1 label to all visits of that diag
    v1_labels <- unique(data.frame(
      STUDY_ID = colData(dds)$STUDY_ID[idx_threshold_calc],
      dynamic_tri_label = colData(dds)[[tri_high_colname]][idx_threshold_calc],
      stringsAsFactors = FALSE
    ))
    
    idx_all <- colData(dds)$DIAGNOSIS == diag
    colData(dds)[[tri_high_colname]][idx_all] <-
      v1_labels$dynamic_tri_label[match(
        colData(dds)$STUDY_ID[idx_all],
        v1_labels$STUDY_ID
      )]
    
    # 4) set factor (Low reference)
    if (! tri_high_colname %in% names(colData(dds))) {
      stop(paste0("Column ", tri_high_colname, " was not created properly."))
    }
    colData(dds)[[tri_high_colname]] <- factor(
      colData(dds)[[tri_high_colname]],
      levels = c("Low", "High")
    )
    
    # store results per diagnosis for later combined plotting
    if (! (tri_variable %in% names(dgea_results_tri))) {
      dgea_results_tri[[tri_variable]] <- list()
      gsea_results_tri[[tri_variable]] <- list()
    }
    dgea_results_tri[[tri_variable]][[diag]] <- list()
    gsea_results_tri[[tri_variable]][[diag]] <- list()
    
    # 5) loop over visits for DGEA and GSEA
    for (visit in VISITS[1:5]) {
      message(" ", diag, " | Visit = ", visit)
      
      dds_sub <- subset_dds(dds, condition = paste0("DIAGNOSIS == '", diag, "'"))
      dds_sub <- subset_dds(dds_sub, condition = paste0("VISIT == '", visit, "'"))
      # colData(dds_sub) <- droplevels(colData(dds_sub)) # should already be dropped by subset_dds function!
      
      tab <- colData(dds)[colData(dds)$DIAGNOSIS == diag & colData(dds)$VISIT == visit,
                          c("SEX", "AGE", tri_high_colname)]
      if (any(!complete.cases(tab)) && nrow(tab) > 0) {
        message("The following sample has an NA value in at least one of the design variables:")
        print(tab[!complete.cases(tab), ])
      }
      
      vars_for_design <- c("SEX", "AGE", tri_high_colname)
      complete_samples <- complete.cases(as.data.frame(colData(dds_sub)[, vars_for_design]))
      
      if (sum(!complete_samples) > 0) {
        message("Removing ", sum(!complete_samples), " sample(s) with missing design variables")
        dds_sub <- dds_sub[, complete_samples]
        colData(dds_sub) <- droplevels(colData(dds_sub))
      }
      
      if (ncol(dds_sub) < 2 || length(unique(colData(dds_sub)[[tri_high_colname]])) < 2) {
        message("Skipping: not enough valid samples or no contrast possible.")
        next
      }
      
      design_formula <- as.formula(paste("~ SEX + AGE +", tri_high_colname))
      design(dds_sub) <- design_formula
      
      dds_file <- paste0("output/", visit, "_", diag, "_", tri_variable, "_dds.rds")
      if (file.exists(dds_file)) {
        dds_sub <- readRDS(dds_file)
      } else {
        dds_sub <- DESeq(dds_sub, parallel = TRUE, BPPARAM = MulticoreParam(2))
        saveRDS(dds_sub, file = dds_file)
      }
      
      dgea_results_tri[[tri_variable]][[diag]][[visit]] <- run_dgea(
        dds = dds_sub,
        conditions = "High",
        method = "contrast",
        factor_name = tri_high_colname,
        ref = "Low"
      )
      
      nes_file <- paste0("output/", visit, "_", diag, "_", tri_variable, "_nes.rds")
      if (file.exists(nes_file)) {
        gsea_results_tri[[tri_variable]][[diag]][[visit]] <- readRDS(nes_file)
      } else {
        gsea_results_tri[[tri_variable]][[diag]][[visit]] <- run_gsea(
          dgea_result = dgea_results_tri[[tri_variable]][[diag]][[visit]],
          pathways = PATHWAYS,
          conditions = "High"
        )
        saveRDS(gsea_results_tri[[tri_variable]][[diag]][[visit]], file = nes_file)
      }
    }
  }
  
  # --- Add COHORT Analysis to the combined dotmatrix ---
  message("\n--- Running COHORT analysis for TRI variable: ", tri_variable, " ---")
  tri_high_colname_cohort <- paste0(tri_variable, "_HIGH")
  
  # step 1: Calculate COHORT threshold for High/Low classification
  idx_cohort_threshold_calc <- colData(dds)$VISIT == "V1" # All V1 samples, regardless of diagnosis
  
  if (sum(idx_cohort_threshold_calc, na.rm = TRUE) == 0) {
    message("  No V1 samples found for COHORT threshold calculation for ", tri_variable, ". Skipping COHORT analysis.")
    next
  }
  
  cohort_tri_values_v1 <- colData(dds)[[tri_variable]][idx_cohort_threshold_calc]
  
  # Calculate COHORT upper quartile threshold for High responders
  cohort_upper_quantile_threshold <- quantile(cohort_tri_values_v1, probs = 0.75, na.rm = TRUE)
  message("  COHORT upper quartile (", round(cohort_upper_quantile_threshold, 3), ") used for ", tri_variable, " (across all V1 samples).")
  
  # step 2: Assign COHORT High/Low labels to ALL relevant samples in colData(dds)
  colData(dds)[[tri_high_colname_cohort]] <- NA_character_ # clear previous per-diagnosis labels
  
  # Assign based on V1 values to the V1 samples:
  colData(dds)[[tri_high_colname_cohort]][idx_cohort_threshold_calc] <- ifelse(
    colData(dds)[[tri_variable]][idx_cohort_threshold_calc] >= cohort_upper_quantile_threshold, "High", "Low"
  )
  message("  COHORT classification: Values in upper quartile are labeled 'High'.")
  
  # propagate these COHORT V1 labels to all visits for each STUDY_ID
  cohort_v1_labels <- unique(data.frame(
    STUDY_ID = colData(dds)$STUDY_ID[idx_cohort_threshold_calc],
    dynamic_tri_label = colData(dds)[[tri_high_colname_cohort]][idx_cohort_threshold_calc],
    stringsAsFactors = FALSE
  ))
  
  colData(dds)[[tri_high_colname_cohort]] <- cohort_v1_labels$dynamic_tri_label[match(
    colData(dds)$STUDY_ID,
    cohort_v1_labels$STUDY_ID
  )]
  
  # set factor (Low reference)
  if (! all(c("Low", "High") %in% unique(na.omit(colData(dds)[[tri_high_colname_cohort]])))) {
    message("  Warning: After COHORT classification, not both 'Low' and 'High' levels exist for ", tri_variable, ". Skipping COHORT analysis for this variable.")
    next
  }
  colData(dds)[[tri_high_colname_cohort]] <- factor(
    colData(dds)[[tri_high_colname_cohort]],
    levels = c("Low", "High")
  )
  
  # Loop over visits for COHORT DGEA and GSEA
  for (visit in VISITS[1:5]) {
    message("  Analyzing COHORT | Visit = ", visit)
    
    dds_sub_cohort <- dds[, colData(dds)$VISIT == visit]
    colData(dds_sub_cohort) <- droplevels(colData(dds_sub_cohort))
    
    # one could experiment with adjusting for diagnosis here as well, but lets keep designs consistent for now
    vars_for_design_cohort <- c("SEX", "AGE", tri_high_colname_cohort)
    tab_cohort <- colData(dds_sub_cohort)[, vars_for_design_cohort]
    
    if (any(!complete.cases(tab_cohort)) && nrow(tab_cohort) > 0) {
      message("    The following COHORT sample(s) have an NA value in at least one of the design variables:")
      print(tab_cohort[!complete.cases(tab_cohort), ])
    }
    
    complete_samples_cohort <- complete.cases(as.data.frame(colData(dds_sub_cohort)[, vars_for_design_cohort]))
    
    if (sum(!complete_samples_cohort) > 0) {
      message("    Removing ", sum(!complete_samples_cohort), " COHORT sample(s) with missing design variables")
      dds_sub_cohort <- dds_sub_cohort[, complete_samples_cohort]
      colData(dds_sub_cohort) <- droplevels(colData(dds_sub_cohort))
    }
    
    if (ncol(dds_sub_cohort) < 2 || length(unique(colData(dds_sub_cohort)[[tri_high_colname_cohort]])) < 2) {
      message("    Skipping COHORT DGEA/GSEA for ", visit, ": Not enough valid samples (", ncol(dds_sub_cohort), ") or not enough levels for contrast (", length(unique(colData(dds_sub_cohort)[[tri_high_colname_cohort]])), ").")
      next
    }
    
    design_formula_cohort <- as.formula(paste("~ SEX + AGE +", tri_high_colname_cohort))
    design(dds_sub_cohort) <- design_formula_cohort
    
    dds_cohort_file <- paste0("output/", visit, "_cohort_", tri_variable, "_dds.rds")
    if (file.exists(dds_cohort_file)) {
      dds_sub_cohort <- readRDS(dds_cohort_file)
      message("    Loaded existing COHORT DESeq2 results for ", visit, ".")
    } else {
      message("    Running COHORT DESeq2 for ", visit, ".")
      dds_sub_cohort <- DESeq(dds_sub_cohort, parallel = TRUE, BPPARAM = MulticoreParam(2))
      saveRDS(dds_sub_cohort, file = dds_cohort_file)
    }
    
    dgea_results_tri_cohort[[tri_variable]][[visit]] <- run_dgea(
      dds         = dds_sub_cohort,
      conditions  = "High",
      method      = "contrast",
      factor_name = tri_high_colname_cohort,
      ref         = "Low"
    )
    
    nes_cohort_file <- paste0("output/", visit, "_cohort_", tri_variable, "_nes.rds")
    if (file.exists(nes_cohort_file)) {
      gsea_results_tri_cohort[[tri_variable]][[visit]] <- readRDS(nes_cohort_file)
      message("    Loaded existing COHORT fgsea results for ", visit, ".")
    } else {
      message("    Running COHORT fgsea for ", visit, ".")
      gsea_results_tri_cohort[[tri_variable]][[visit]] <- run_gsea(
        dgea_result = dgea_results_tri_cohort[[tri_variable]][[visit]],
        pathways    = PATHWAYS,
        conditions  = "High"
      )
      saveRDS(gsea_results_tri_cohort[[tri_variable]][[visit]], file = nes_cohort_file)
    }
  }
  
  # combine GSEA results for plotting multiple diagnoses side by side for each visit, including COHORT
  for (visit in VISITS[1:5]) {
    message("Plotting combined dotmatrix for tri variable: ", tri_variable, " at visit: ", visit)
    
    combined_gsea_results <- list()
    # Add diagnosis-specific results
    for (diag in DIAGNOSES) {
      if (!is.null(gsea_results_tri[[tri_variable]][[diag]][[visit]])) {
        gsea_res_temp <- gsea_results_tri[[tri_variable]][[diag]][[visit]]
        colnames(gsea_res_temp) <- gsub("NES.High", paste0("NES.", diag), colnames(gsea_res_temp))
        colnames(gsea_res_temp) <- gsub("padj.High", paste0("padj.", diag), colnames(gsea_res_temp))
        combined_gsea_results[[diag]] <- gsea_res_temp
      }
    }
    
    # Add COHORT results
    if (!is.null(gsea_results_tri_cohort[[tri_variable]][[visit]])) {
      gsea_res_cohort_temp <- gsea_results_tri_cohort[[tri_variable]][[visit]]
      colnames(gsea_res_cohort_temp) <- gsub("NES.High", "NES.COHORT", colnames(gsea_res_cohort_temp))
      colnames(gsea_res_cohort_temp) <- gsub("padj.High", "padj.COHORT", colnames(gsea_res_cohort_temp))
      combined_gsea_results[["COHORT"]] <- gsea_res_cohort_temp
    }
    
    # if there are results to combine, merge them into a single data frame
    if (length(combined_gsea_results) > 0) {
      merged_gsea_df <- combined_gsea_results[[1]]
      
      if (length(combined_gsea_results) > 1) {
        for (i in 2:length(combined_gsea_results)) {
          merged_gsea_df <- full_join(merged_gsea_df, combined_gsea_results[[i]], by = "PATHWAY")
        }
      }
      
      # ensure conditions are ordered as in DIAGNOSES, with "COHORT" at the end
      conditions_for_plot <- c(DIAGNOSES[DIAGNOSES %in% sub("NES.", "", grep("^NES\\.", colnames(merged_gsea_df), value = TRUE))], "COHORT")
      conditions_for_plot <- unique(conditions_for_plot[conditions_for_plot %in% sub("NES.", "", grep("^NES\\.", colnames(merged_gsea_df), value = TRUE))]) # Ensure only existing conditions are kept and unique
      
      if (length(conditions_for_plot) > 0 && nrow(merged_gsea_df) > 0) {
        dm <- plot_gsea_dotmatrix(
          gsea_result = merged_gsea_df,
          conditions = conditions_for_plot,
          padj_threshold  = 0.01,
          pathway2genes = PATHWAYS,
          cluster = TRUE,
          cluster_height  = 0.99
        )
        if (!is.null(dm)) {
          dm <- dm + force_panelsizes(cols = unit(1.5, "cm")) +
            ggtitle(paste0(tri_variable, " High vs Low by Group (Diagnosis + COHORT) at ", visit))
          
          ggsave(
            filename = paste0("figs/analysis/gsea/tri/", tri_variable, "_", visit, "_combined_and_cohort_dotmatrix.png"),
            plot = dm,
            width = 15, height = 6, dpi = 300 # adjusted width for 5 columns
          )
        }
      } else {
        message("Skipping combined dotmatrix for ", tri_variable, " at ", visit, ": No valid data to plot.")
      }
    }
  }
}
# DESeq2 TRI_ELISA lower quartile Low vs High across all diagnoses + the cohort & visits ~ SEX + AGE + TRI_ELISA_LOW ---------------

tri_variables = c("TRI_ELISA_V3", "TRI_ELISA_V4", "TRI_ELISA_V5")

# Assuming 'dds' and 'colData(dds)' are already defined from your main script
metadata <- as.data.frame(colData(dds))

for (tri_variable in tri_variables) {
  # dynamically create the label for the new column based on the current tri_variable
  tri_low_colname <- paste0(tri_variable, "_LOW")
  
  dgea_results_tri_low <- list()
  gsea_results_tri_low <- list()
  
  dgea_results_tri_cohort_low <- list()
  gsea_results_tri_cohort_low <- list()
  
  # initialize the dynamic TRI_ELISA_VX_LOW column
  colData(dds)[[tri_low_colname]] <- NA_character_
  
  for (diag in DIAGNOSES) {
    message("Processing diagnosis: ", diag, " for TRI variable: ", tri_variable, " (Lower Quartile)")
    
    # 1) compute threshold for tri within this diagnosis at V1 (to not count any subjects twice)
    idx_threshold_calc <- colData(dds)$DIAGNOSIS == diag & colData(dds)$VISIT == "V1"
    
    # calculate lower quartile threshold for Low responders
    lower_quantile_threshold <- quantile(colData(dds)[[tri_variable]][idx_threshold_calc], probs = 0.25, na.rm = TRUE)
    message("Using lower quartile (", round(lower_quantile_threshold, 3), ") for ", tri_variable, " in ", diag)
    
    # Assign Low/High for V1 samples in this diag:
    # All values in the lower quartile get "Low", otherwise "High".
    colData(dds)[[tri_low_colname]][idx_threshold_calc] <- ifelse(
      colData(dds)[[tri_variable]][idx_threshold_calc] <= lower_quantile_threshold, "Low", "High"
    )
    
    # 3) propagate each STUDY_ID’s V1 label to all visits of that diag
    v1_labels <- unique(data.frame(
      STUDY_ID = colData(dds)$STUDY_ID[idx_threshold_calc],
      dynamic_tri_label = colData(dds)[[tri_low_colname]][idx_threshold_calc],
      stringsAsFactors = FALSE
    ))
    
    idx_all <- colData(dds)$DIAGNOSIS == diag
    colData(dds)[[tri_low_colname]][idx_all] <-
      v1_labels$dynamic_tri_label[match(
        colData(dds)$STUDY_ID[idx_all],
        v1_labels$STUDY_ID
      )]
    
    # 4) set factor (High reference, Low as contrast)
    if (! tri_low_colname %in% names(colData(dds))) {
      stop(paste0("Column ", tri_low_colname, " was not created properly."))
    }
    colData(dds)[[tri_low_colname]] <- factor(
      colData(dds)[[tri_low_colname]],
      levels = c("High", "Low") # High as reference, Low as the group of interest
    )
    
    # store results per diagnosis for later combined plotting
    if (! (tri_variable %in% names(dgea_results_tri_low))) {
      dgea_results_tri_low[[tri_variable]] <- list()
      gsea_results_tri_low[[tri_variable]] <- list()
    }
    dgea_results_tri_low[[tri_variable]][[diag]] <- list()
    gsea_results_tri_low[[tri_variable]][[diag]] <- list()
    
    # 5) loop over visits for DGEA and GSEA
    for (visit in VISITS[1:5]) {
      message(" ", diag, " | Visit = ", visit)
      
      dds_sub <- subset_dds(dds, condition = paste0("DIAGNOSIS == '", diag, "'"))
      dds_sub <- subset_dds(dds_sub, condition = paste0("VISIT == '", visit, "'"))
      # colData(dds_sub) <- droplevels(colData(dds_sub)) # should already be dropped by subset_dds function!
      
      tab <- colData(dds)[colData(dds)$DIAGNOSIS == diag & colData(dds)$VISIT == visit,
                          c("SEX", "AGE", tri_low_colname)]
      if (any(!complete.cases(tab)) && nrow(tab) > 0) {
        message("The following sample has an NA value in at least one of the design variables:")
        print(tab[!complete.cases(tab), ])
      }
      
      vars_for_design <- c("SEX", "AGE", tri_low_colname)
      complete_samples <- complete.cases(as.data.frame(colData(dds_sub)[, vars_for_design]))
      
      if (sum(!complete_samples) > 0) {
        message("Removing ", sum(!complete_samples), " sample(s) with missing design variables")
        dds_sub <- dds_sub[, complete_samples]
        colData(dds_sub) <- droplevels(colData(dds_sub))
      }
      
      if (ncol(dds_sub) < 2 || length(unique(na.omit(colData(dds_sub)[[tri_low_colname]]))) < 2) {
        message("Skipping: not enough valid samples or no contrast possible for ", tri_low_colname, ".")
        next
      }
      
      design_formula <- as.formula(paste("~ SEX + AGE +", tri_low_colname))
      design(dds_sub) <- design_formula
      
      dds_file <- paste0("output/", visit, "_", diag, "_", tri_variable, "_low_dds.rds") # Modified filename
      if (file.exists(dds_file)) {
        dds_sub <- readRDS(dds_file)
      } else {
        dds_sub <- DESeq(dds_sub, parallel = TRUE, BPPARAM = MulticoreParam(2))
        saveRDS(dds_sub, file = dds_file)
      }
      
      dgea_results_tri_low[[tri_variable]][[diag]][[visit]] <- run_dgea(
        dds = dds_sub,
        conditions = "Low", # Changed to Low
        method = "contrast",
        factor_name = tri_low_colname,
        ref = "High" # Changed to High
      )
      
      nes_file <- paste0("output/", visit, "_", diag, "_", tri_variable, "_low_nes.rds") # Modified filename
      if (file.exists(nes_file)) {
        gsea_results_tri_low[[tri_variable]][[diag]][[visit]] <- readRDS(nes_file)
      } else {
        gsea_results_tri_low[[tri_variable]][[diag]][[visit]] <- run_gsea(
          dgea_result = dgea_results_tri_low[[tri_variable]][[diag]][[visit]],
          pathways = PATHWAYS,
          conditions = "Low" # Changed to Low
        )
        saveRDS(gsea_results_tri_low[[tri_variable]][[diag]][[visit]], file = nes_file)
      }
    }
  }
  
  # --- Add COHORT Analysis to the combined dotmatrix ---
  message("\n--- Running COHORT analysis for TRI variable: ", tri_variable, " (Lower Quartile) ---")
  tri_low_colname_cohort <- paste0(tri_variable, "_LOW")
  
  # step 1: Calculate COHORT threshold for Low/High classification
  idx_cohort_threshold_calc <- colData(dds)$VISIT == "V1" # All V1 samples, regardless of diagnosis
  
  if (sum(idx_cohort_threshold_calc, na.rm = TRUE) == 0) {
    message("  No V1 samples found for COHORT threshold calculation for ", tri_variable, ". Skipping COHORT analysis.")
    next
  }
  
  cohort_tri_values_v1 <- colData(dds)[[tri_variable]][idx_cohort_threshold_calc]
  
  # Calculate COHORT lower quartile threshold for Low responders
  cohort_lower_quantile_threshold <- quantile(cohort_tri_values_v1, probs = 0.25, na.rm = TRUE)
  message("  COHORT lower quartile (", round(cohort_lower_quantile_threshold, 3), ") used for ", tri_variable, " (across all V1 samples).")
  
  # step 2: Assign COHORT Low/High labels to ALL relevant samples in colData(dds)
  colData(dds)[[tri_low_colname_cohort]] <- NA_character_ # clear previous per-diagnosis labels
  
  # Assign based on V1 values to the V1 samples:
  colData(dds)[[tri_low_colname_cohort]][idx_cohort_threshold_calc] <- ifelse(
    colData(dds)[[tri_variable]][idx_cohort_threshold_calc] <= cohort_lower_quantile_threshold, "Low", "High"
  )
  message("  COHORT classification: Values in lower quartile are labeled 'Low'.")
  
  # propagate these COHORT V1 labels to all visits for each STUDY_ID
  cohort_v1_labels <- unique(data.frame(
    STUDY_ID = colData(dds)$STUDY_ID[idx_cohort_threshold_calc],
    dynamic_tri_label = colData(dds)[[tri_low_colname_cohort]][idx_cohort_threshold_calc],
    stringsAsFactors = FALSE
  ))
  
  colData(dds)[[tri_low_colname_cohort]] <- cohort_v1_labels$dynamic_tri_label[match(
    colData(dds)$STUDY_ID,
    cohort_v1_labels$STUDY_ID
  )]
  
  # set factor (High reference)
  if (! all(c("Low", "High") %in% unique(na.omit(colData(dds)[[tri_low_colname_cohort]])))) {
    message("  Warning: After COHORT classification, not both 'Low' and 'High' levels exist for ", tri_variable, ". Skipping COHORT analysis for this variable.")
    next
  }
  colData(dds)[[tri_low_colname_cohort]] <- factor(
    colData(dds)[[tri_low_colname_cohort]],
    levels = c("High", "Low") # High as reference, Low as the group of interest
  )
  
  # Loop over visits for COHORT DGEA and GSEA
  for (visit in VISITS[1:5]) {
    message("  Analyzing COHORT | Visit = ", visit)
    
    dds_sub_cohort <- dds[, colData(dds)$VISIT == visit]
    colData(dds_sub_cohort) <- droplevels(colData(dds_sub_cohort))
    
    # one could experiment with adjusting for diagnosis here as well, but lets keep designs consistent for now
    vars_for_design_cohort <- c("SEX", "AGE", tri_low_colname_cohort)
    tab_cohort <- colData(dds_sub_cohort)[, vars_for_design_cohort]
    
    if (any(!complete.cases(tab_cohort)) && nrow(tab_cohort) > 0) {
      message("    The following COHORT sample(s) have an NA value in at least one of the design variables:")
      print(tab_cohort[!complete.cases(tab_cohort), ])
    }
    
    complete_samples_cohort <- complete.cases(as.data.frame(colData(dds_sub_cohort)[, vars_for_design_cohort]))
    
    if (sum(!complete_samples_cohort) > 0) {
      message("    Removing ", sum(!complete_samples_cohort), " COHORT sample(s) with missing design variables")
      dds_sub_cohort <- dds_sub_cohort[, complete_samples_cohort]
      colData(dds_sub_cohort) <- droplevels(colData(dds_sub_cohort))
    }
    
    if (ncol(dds_sub_cohort) < 2 || length(unique(na.omit(colData(dds_sub_cohort)[[tri_low_colname_cohort]]))) < 2) {
      message("    Skipping COHORT DGEA/GSEA for ", visit, ": Not enough valid samples (", ncol(dds_sub_cohort), ") or not enough levels for contrast (", length(unique(na.omit(colData(dds_sub_cohort)[[tri_low_colname_cohort]]))), ").")
      next
    }
    
    design_formula_cohort <- as.formula(paste("~ SEX + AGE +", tri_low_colname_cohort))
    design(dds_sub_cohort) <- design_formula_cohort
    
    dds_cohort_file <- paste0("output/", visit, "_cohort_", tri_variable, "_low_dds.rds") # Modified filename
    if (file.exists(dds_cohort_file)) {
      dds_sub_cohort <- readRDS(dds_cohort_file)
      message("    Loaded existing COHORT DESeq2 results for ", visit, ".")
    } else {
      message("    Running COHORT DESeq2 for ", visit, ".")
      dds_sub_cohort <- DESeq(dds_sub_cohort, parallel = TRUE, BPPARAM = MulticoreParam(2))
      saveRDS(dds_sub_cohort, file = dds_cohort_file)
    }
    
    dgea_results_tri_cohort_low[[tri_variable]][[visit]] <- run_dgea(
      dds         = dds_sub_cohort,
      conditions  = "Low", # Changed to Low
      method      = "contrast",
      factor_name = tri_low_colname_cohort,
      ref         = "High" # Changed to High
    )
    
    nes_cohort_file <- paste0("output/", visit, "_cohort_", tri_variable, "_low_nes.rds") # Modified filename
    if (file.exists(nes_cohort_file)) {
      gsea_results_tri_cohort_low[[tri_variable]][[visit]] <- readRDS(nes_cohort_file)
      message("    Loaded existing COHORT fgsea results for ", visit, ".")
    } else {
      message("    Running COHORT fgsea for ", visit, ".")
      gsea_results_tri_cohort_low[[tri_variable]][[visit]] <- run_gsea(
        dgea_result = dgea_results_tri_cohort_low[[tri_variable]][[visit]],
        pathways    = PATHWAYS,
        conditions  = "Low" # Changed to Low
      )
      saveRDS(gsea_results_tri_cohort_low[[tri_variable]][[visit]], file = nes_cohort_file)
    }
  }
  
  # combine GSEA results for plotting multiple diagnoses side by side for each visit, including COHORT
  for (visit in VISITS[1:5]) {
    message("Plotting combined dotmatrix for tri variable: ", tri_variable, " (Lower Quartile) at visit: ", visit)
    
    combined_gsea_results_low <- list()
    # Add diagnosis-specific results
    for (diag in DIAGNOSES) {
      if (!is.null(gsea_results_tri_low[[tri_variable]][[diag]][[visit]])) {
        gsea_res_temp <- gsea_results_tri_low[[tri_variable]][[diag]][[visit]]
        colnames(gsea_res_temp) <- gsub("NES.Low", paste0("NES.", diag), colnames(gsea_res_temp)) # Changed to Low
        colnames(gsea_res_temp) <- gsub("padj.Low", paste0("padj.", diag), colnames(gsea_res_temp)) # Changed to Low
        combined_gsea_results_low[[diag]] <- gsea_res_temp
      }
    }
    
    # Add COHORT results
    if (!is.null(gsea_results_tri_cohort_low[[tri_variable]][[visit]])) {
      gsea_res_cohort_temp <- gsea_results_tri_cohort_low[[tri_variable]][[visit]]
      colnames(gsea_res_cohort_temp) <- gsub("NES.Low", "NES.COHORT", colnames(gsea_res_cohort_temp)) # Changed to Low
      colnames(gsea_res_cohort_temp) <- gsub("padj.Low", "padj.COHORT", colnames(gsea_res_cohort_temp)) # Changed to Low
      combined_gsea_results_low[["COHORT"]] <- gsea_res_cohort_temp
    }
    
    # if there are results to combine, merge them into a single data frame
    if (length(combined_gsea_results_low) > 0) {
      merged_gsea_df_low <- combined_gsea_results_low[[1]]
      
      if (length(combined_gsea_results_low) > 1) {
        for (i in 2:length(combined_gsea_results_low)) {
          merged_gsea_df_low <- full_join(merged_gsea_df_low, combined_gsea_results_low[[i]], by = "PATHWAY")
        }
      }
      
      # ensure conditions are ordered as in DIAGNOSES, with "COHORT" at the end
      conditions_for_plot <- c(DIAGNOSES[DIAGNOSES %in% sub("NES.", "", grep("^NES\\.", colnames(merged_gsea_df_low), value = TRUE))], "COHORT")
      conditions_for_plot <- unique(conditions_for_plot[conditions_for_plot %in% sub("NES.", "", grep("^NES\\.", colnames(merged_gsea_df_low), value = TRUE))]) # Ensure only existing conditions are kept and unique
      
      if (length(conditions_for_plot) > 0 && nrow(merged_gsea_df_low) > 0) {
        dm <- plot_gsea_dotmatrix(
          gsea_result = merged_gsea_df_low,
          conditions = conditions_for_plot,
          padj_threshold  = 0.01,
          pathway2genes = PATHWAYS,
          cluster = TRUE,
          cluster_height  = 0.99
        )
        if (!is.null(dm)) {
          dm <- dm + force_panelsizes(cols = unit(1.5, "cm")) +
            ggtitle(paste0(tri_variable, " Low vs High by Group (Diagnosis + COHORT) at ", visit)) # Changed title
          
          ggsave(
            filename = paste0("figs/analysis/gsea/tri/", tri_variable, "_", visit, "_low_combined_and_cohort_dotmatrix.png"), # Modified filename
            plot = dm,
            width = 15, height = 6, dpi = 300 # adjusted width for 5 columns
          )
        }
      } else {
        message("Skipping combined dotmatrix for ", tri_variable, " (Lower Quartile) at ", visit, ": No valid data to plot.")
      }
    }
  }
}

# DESeq2 analysis for NSCLC ~ CHEMOTHERAPY_AT_V1 --------------------------


nsclc_dds <- subset_dds(dds, condition = "DIAGNOSIS == 'NSCLC'")
design(nsclc_dds) <- ~ CHEMOTHERAPY_AT_V1

dgea_results_NSCLC <- list()
gsea_results_NSCLC <- list()

for (visit in VISITS) {
  message(" NSCLC | Visit = ", visit)
  
  dds_v <- nsclc_dds[, nsclc_dds$VISIT == visit]
  colData(dds_v) <- droplevels(colData(dds_v))
  
  # DESeq2
  dds_file <- paste0("output/", visit, "_NSCLC_CHEMO_dds.rds")
  if (file.exists(dds_file)) {
    dds_v <- readRDS(dds_file)
  } else {
    dds_v <- DESeq(dds_v, parallel=TRUE, BPPARAM=MulticoreParam(2))
    saveRDS(dds_v, dds_file)
  }
  
  # run_dgea: only "Y" vs "N"
  dgea_results_NSCLC[[visit]] <- run_dgea(
    dds         = dds_v,
    conditions  = "Y",
    method      = "contrast",
    factor_name = "CHEMOTHERAPY_AT_V1",
    ref         = "N"
  )
  
  # run_gsea: only "Y"
  nes_file <- paste0("output/", visit, "_NSCLC_CHEMO_nes.rds")
  if (file.exists(nes_file)) {
    gsea_results_NSCLC[[visit]] <- readRDS(nes_file)
  } else {
    gsea_results_NSCLC[[visit]] <- run_gsea(
      dgea_result = dgea_results_NSCLC[[visit]],
      pathways    = PATHWAYS,
      conditions  = "Y"
    )
    saveRDS(gsea_results_NSCLC[[visit]], nes_file)
  }
  
  # # heatmap for NES.Y
  # ht <- plot_gsea_heatmap(
  #   gsea_results_NSCLC[[visit]],
  #   conditions = "Y"
  # )
  # png(
  #   filename = paste0("figs/analysis/gsea/cancer_treatment/NSCLC_CHEMO_", visit, "_heatmap.png"),
  #   width    = 1400, height = 1000
  # )
  # draw(ht,
  #      padding         = unit(c(10,80,10,20), "mm"),
  #      column_title    = paste("NSCLC CHEMO GSEA at", visit),
  #      column_title_gp = gpar(fontsize=16, fontface="bold"))
  # dev.off()
  
  # dotmatrix for padj.Y / NES.Y
  dm <- plot_gsea_dotmatrix(
    gsea_result   = gsea_results_NSCLC[[visit]],
    conditions    = "Y",
    pathway2genes = PATHWAYS,
    cluster       = TRUE,
    cluster_height= 0.96
  ) + force_panelsizes(cols = unit(1.5, "cm")) +
    ggtitle(paste("NSCLC CHEMO at", visit, "vs no chemo"))
  
  ggsave(
    filename = paste0("figs/analysis/gsea/cancer_treatment/NSCLC_CHEMO_", visit, "_dotmatrix.png"),
    plot     = dm,
    width    = 8, height = 5, dpi = 300
  )
}

# DESeq2 analysis for MM ~ LENALIDOMIDE_AT_V1 -----------------------------


# subset data to MM
mm_dds <- subset_dds(dds, condition = "DIAGNOSIS == 'MM'")

# set design
design(mm_dds) <- ~ LENALIDOMIDE_AT_V1

# prepare result containers
dgea_results_MM_LENA <- list()
gsea_results_MM_LENA <- list()

# loop over visits
for (visit in VISITS) {
  message(" MM | Visit =", visit)
  
  # subset by visit and drop unused levels
  dds_v <- mm_dds[, mm_dds$VISIT == visit]
  colData(dds_v) <- droplevels(colData(dds_v))
  
  # DESeq2
  dds_file <- paste0("output/", visit, "_MM_LENA_dds.rds")
  if (file.exists(dds_file)) {
    dds_v <- readRDS(dds_file)
  } else {
    dds_v <- DESeq(dds_v, parallel = TRUE, BPPARAM = MulticoreParam(2))
    saveRDS(dds_v, file = dds_file)
  }
  
  # differential expression (Y vs N)
  dgea_results_MM_LENA[[visit]] <- run_dgea(
    dds          = dds_v,
    conditions   = "Y",                   # only Y vs N
    method       = "contrast",
    factor_name  = "LENALIDOMIDE_AT_V1",
    ref          = "N"
  )
  
  # gene set enrichment (only for Y)
  nes_file <- paste0("output/", visit, "_MM_LENA_nes.rds")
  if (file.exists(nes_file)) {
    gsea_results_MM_LENA[[visit]] <- readRDS(nes_file)
  } else {
    gsea_results_MM_LENA[[visit]] <- run_gsea(
      dgea_result = dgea_results_MM_LENA[[visit]],
      pathways    = PATHWAYS,
      conditions  = "Y"
    )
    saveRDS(gsea_results_MM_LENA[[visit]], file = nes_file)
  }
  
  # # plot GSEA heatmap for NES.Y
  # ht <- plot_gsea_heatmap(
  #   gsea_results_MM_LENA[[visit]],
  #   conditions = "Y"
  # )
  # png(
  #   filename = paste0("figs/analysis/gsea/cancer_treatment/MM_LENA_", visit, "_heatmap.png"),
  #   width    = 1400,
  #   height   = 1000
  # )
  # draw(
  #   ht,
  #   padding         = unit(c(10, 80, 10, 20), "mm"),
  #   column_title    = paste("MM LENALIDOMIDE GSEA at", visit),
  #   column_title_gp = gpar(fontsize = 16, fontface = "bold")
  # )
  # dev.off()
  
  # plot GSEA dotmatrix for padj.Y / NES.Y
  dm <- plot_gsea_dotmatrix(
    gsea_result     = gsea_results_MM_LENA[[visit]],
    conditions      = "Y",
    pathway2genes   = PATHWAYS,
    cluster         = TRUE,
    cluster_height  = 0.96
  ) + force_panelsizes(cols = unit(1.5, "cm")) +
    ggtitle(paste("MM LENALIDOMIDE at", visit, "vs no treatment"))
  
  ggsave(
    filename = paste0("figs/analysis/gsea/cancer_treatment/MM_LENA_", visit, "_dotmatrix.png"),
    plot     = dm,
    width    = 8,
    height   = 5,
    dpi      = 300
  )
}


# DGEA upset plots --------------------------------------------------------------

for (visit in VISITS[-1]) {
  
  if (visit %in% c("V4", "V5", "V6")) next # these visits have no significant results
  
  upset_result <- plot_significant_degs_upset(dgea_results, visit)
  
  png(paste0("figs/analysis/dgea/", visit, "_upset.png"), width = 2800, height = 2200, res = 300)
  print(upset_result$plot)
  grid.text(paste0("Significant DEGs UpsetPlot ", visit), x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
  dev.off()
  
}

# DGEA volcano plots ----------------------------------------

for (list in names(dgea_results)) {
  for (sublist in names(dgea_results[[list]])) {
    v <- plot_degs_volcano(dgea_results[[list]][[sublist]][["results"]], y = "padj") + labs(title = paste(list, sublist))
    ggsave(paste0("figs/analysis/dgea/volcano/", list, "_", sublist, "_degs_volcano.png"), plot = v, width = 8, height = 8)
  }
}

message("DGEA and GSEA completed. All plots saved.")

# subject log2FoldChange ANOVA --------------------------------------------

# different ANOVA models are possible, let's try some and rank them via the Akaike Information Criterion
# this criterion ranks those models best that explain the greatest amount of variation using the fewest possible independent variables

# for each visit, loop over all the genes and perform ANOVA with diagnoses as groups and log2FoldChange in normalized gene counts between V1 and VX

genes <- rownames(l2fc_data_imputed$count_matrix)

anova_results <- list()

# function to visualize anova results
plot_anova_pvalues <- function(anova_df,
                               n_top = 100,
                               pval_threshold = 0.05,
                               title = "ANOVA p-values") {
  # optionally filter by p-value threshold
  if (!is.null(pval_threshold)) {
    anova_df <- subset(anova_df, P_VALUE <= pval_threshold)
  }
  
  # sort by P_VALUE ascending, then take top n_top
  top_df <- anova_df[order(anova_df$P_VALUE), , drop=FALSE]
  if (nrow(top_df) > n_top) {
    top_df <- top_df[1:n_top, ]
  }
  
  # build the barplot of –log10(p)
  p <- ggplot(top_df, aes(x = reorder(GENE, P_VALUE), y = -log10(P_VALUE))) +
    geom_col(fill = "steelblue") +
    xlab("Gene") +
    ylab(expression(-log[10](P-value))) +
    ggtitle(title) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # return both plot and the table of genes/p-values
  return(list(
    plot  = p,
    genes = top_df
  ))
}

make_anova_entry <- function(g, l2fc_data_subset) {
  
  # extract expression
  expr <- as.numeric(l2fc_data_subset$count_matrix[g, ])
  
  # build a temporary df with expression and diagnosis as input for ANOVA
  temp_df <- data.frame(
    EXPRESSION = expr,
    DIAGNOSIS = l2fc_data_subset$metadata[, "DIAGNOSIS"],
    AGE = l2fc_data_subset$metadata[, "AGE"],
    STUDY = l2fc_data_subset$metadata[, "STUDY"]
  )
  
  one_way_model <- aov(EXPRESSION ~ DIAGNOSIS, data = temp_df)       
  two_way_model <- aov(EXPRESSION ~ DIAGNOSIS + AGE, data = temp_df)
  interaction_model <- aov(EXPRESSION ~ DIAGNOSIS*AGE, data = temp_df)
  blocking_model <- aov(EXPRESSION ~ DIAGNOSIS + AGE + STUDY, data = temp_df)
  
  model_set <- list(one_way_model, two_way_model, interaction_model, blocking_model)
  model_names <- c("one_way", "two_way", "interaction", "blocking")
  
  aic  <- aictab(model_set, modnames = model_names) 
  
  list(
    one_way_model = one_way_model,
    two_way_model = two_way_model,
    interaction_model = interaction_model,
    blocking_model = blocking_model,
    aic  = aic
  )
}

if (COMPARE_ANOVA_MODELS == TRUE) {
  
  # calculate 4 gene-level models
  # takes around 15 - 20 minutes, requires several gigabytes of R memory
  for (visit in VISITS[-1]) {
    
    # subset data
    l2fc_data_subset <- subset_l2fc_data(l2fc_data_imputed, condition = paste0("VISIT ==", "'", visit, "'"))
    
    res_list <- lapply(genes, make_anova_entry, l2fc_data_subset = l2fc_data_subset)
    names(res_list) <- genes
    
    anova_results[[visit]] <- res_list
  }
  
  # analyse which model came out best
  # create a histogram of AIC scores for each model and gene at a visit
  for (visit in VISITS[-1]) {
    res_list <- anova_results[[visit]]
    
    # build a long data.frame: one row per (gene, model)
    df <- do.call(rbind, lapply(names(res_list), function(g) {
      aic_obj <- res_list[[g]]$aic
      
      # take the raw AICc numbers
      aicc_vals <- aic_obj$AICc
      
      # and name them from Modnames
      names(aicc_vals) <- aic_obj$Modnames
      
      # build your data.frame
      data.frame(
        gene  = g,
        model = names(aicc_vals),
        AICc  = as.numeric(aicc_vals),
        stringsAsFactors = FALSE
      )
    }))
    
    
    # compute mean AICc per model, dropping -Inf
    mean_aicc <- df %>%
      filter(is.finite(AICc)) %>%
      group_by(model) %>%
      dplyr::summarize(mean_AICc = mean(AICc), .groups="drop") %>%
      mutate(
        # prepare label text
        label = paste0("mean = ", round(mean_AICc, 2))
      )
    
    # histogram + mean label
    ggplot(df %>% filter(is.finite(AICc)), aes(x = AICc)) +
      geom_histogram(bins = 250, fill = "grey80", color = "black") +
      facet_wrap(~ model, scales = "free") +
      # add text labels: x at right edge (Inf), y at top (Inf)
      geom_text(
        data    = mean_aicc,
        aes(x = Inf, y = Inf, label = label),
        hjust   = 1.1,        # push slightly left from right edge
        vjust   = 1.5,        # push slightly down from top
        size    = 3           # text size
      ) +
      labs(
        title    = paste("AIC score distribution by model at", visit),
        x        = "AIC score",
        y        = "Number of genes"
      ) +
      theme_bw()
    
    ggsave(paste0("figs/analysis/subject_l2fc_anova/", visit, "_model_aic_scores_histogram.png"))
  }
  
  # if we look at all genes and compare the models, they perform very similarly
  # it seems like the simple one-way model has a slight edge in the mean - we will therefore stick with it for now and retrieve results
  for (visit in VISITS[-1]) {
    
    # subset your data for this visit
    l2fc_data_subset <- subset_l2fc_data(
      l2fc_data_imputed,
      condition = paste0("VISIT == '", visit, "'")
    )
    
    # pull out the list of one_way models for each gene
    res_list <- anova_results[[visit]]
    
    # build a data.frame of Gene, F_value and P_VALUE from the one_way model
    anova_res_df <- do.call(
      rbind,
      lapply(names(res_list), function(gene) {
        # extract the aov object
        aov_res     <- res_list[[gene]]$one_way_model
        aov_summary <- summary(aov_res)
        
        # grab the row corresponding to 'DIAGNOSIS'
        fval <- aov_summary[[1]]["DIAGNOSIS", "F value"]
        pval <- aov_summary[[1]]["DIAGNOSIS", "Pr(>F)"]
        
        data.frame(
          GENE    = gene,
          F_value = fval,
          P_VALUE = pval,
          stringsAsFactors = FALSE
        )
      })
    )
    
    # now plot top 100 by p-value
    top_anova <- plot_anova_pvalues(anova_res_df, n_top = 100)
    top_anova$plot +
      ggtitle(paste0(
        "ANOVA p-values by gene, comparison of log2FoldChange in normalized counts by diagnoses at ",
        visit
      ))
    ggsave(
      filename = paste0("figs/analysis/subject_l2fc_anova/", visit, "_top_anova_plot.png"),
      height   = 6,
      width    = 10
    )
    
    # GO-term analysis on those top genes
    ego <- dotplot_enrichment(top_anova$genes$GENE, gmt_file = BTM_FILE)
    if (!is.null(ego$plot)) {
      ggsave(
        file = paste0("figs/analysis/subject_l2fc_anova/", visit, "_top_anova_goterms_dotplot.png"),
        plot = ego$plot
      )
    }
    
    # boxplots for the top 10 genes
    for (top_gene in top_anova$genes$GENE[1:10]) {
      p_box <- boxplot_gene_count_with_test(
        l2fc_data_subset$count_matrix,
        l2fc_data_subset$metadata,
        top_gene,
        "DIAGNOSIS",
        comparisons = list(
          c("HC", "NSCLC"),
          c("HC", "MM"),
          c("HC", "SMM")
        )
      ) +
        ggtitle(paste0(
          top_gene,
          " l2fc of gene counts from V1 to ",
          visit
        )) +
        xlab("Diagnosis") +
        ylab(paste0("V1 to ", visit, " l2fc of gene counts"))
      
      ggsave(
        filename = paste0(
          "figs/analysis/subject_l2fc_anova/",
          visit, "_", top_gene, "_boxplot.png"
        ),
        plot = p_box,
        height = 6,
        width  = 8
      )
    }
  }
} else {
  
  # compute anova for just the simple one–way model
  for (visit in VISITS[-1]) {
    
    # subset data for this visit
    l2fc_data_subset <- subset_l2fc_data(
      l2fc_data_imputed,
      condition = paste0("VISIT == '", visit, "'")
    )
    
    # initialize empty results dataframe for this visit
    anova_res_df <- data.frame(
      GENE    = character(0),
      F_value = numeric(0),
      P_VALUE = numeric(0),
      stringsAsFactors = FALSE
    )
    
    # loop over genes and row-bind each result
    for (g in genes) {
      
      # extract expression vector
      expr <- as.numeric(l2fc_data_subset$count_matrix[g, ])
      
      # build df for aov
      temp_df <- data.frame(
        EXPRESSION = expr,
        DIAGNOSIS  = l2fc_data_subset$metadata[, "DIAGNOSIS"]
      )
      
      # one-way ANOVA
      one_way_model <- aov(EXPRESSION ~ DIAGNOSIS, data = temp_df)
      aov_summary   <- summary(one_way_model)[[1]]
      
      # pull F and p
      fval <- aov_summary["DIAGNOSIS", "F value"]
      pval <- aov_summary["DIAGNOSIS", "Pr(>F)"]
      
      # append to this visit’s results
      anova_res_df <- rbind(
        anova_res_df,
        data.frame(
          GENE    = g,
          F_value = fval,
          P_VALUE = pval,
          stringsAsFactors = FALSE
        )
      )
    }
    
    # sort by P_VALUE ascending
    anova_res_df <- anova_res_df[order(anova_res_df$P_VALUE), , drop=FALSE]
    
    # store in list
    anova_results[[visit]] <- anova_res_df
    
    # save data as .csv
    write.csv(
      anova_res_df,
      file      = paste0("cleaned_output/subject_l2fc_anova/", visit, "_anova_results.csv"),
      row.names = FALSE
    )
  }
}

# loop over visits to plot anova results
for (visit in VISITS[-1]) {
  
  anova_res_df <- anova_results[[visit]]
  
  # now do your plotting/enrichment using anova_res_df
  top_anova <- plot_anova_pvalues(anova_res_df,
                                  n_top = 50,
                                  pval_threshold = 0.001,
                                  title = "Top 100 Genes with p ≤ 0.001")
  top_anova$plot +
    ggtitle(paste0(
      "ANOVA p-values by gene, comparison of log2FoldChange in normalized counts by diagnoses at ",
      visit
    ))
  ggsave(
    filename = paste0("figs/analysis/subject_l2fc_anova/", visit, "_top_anova_plot.png"),
    height   = 6,
    width    = 10,
    dpi = 300
  )
  
  
  # perform enrichment just based on p value and do not limit number of genes
  top_genes <- anova_res_df[anova_res_df$P_VALUE < 0.001, ]$GENE
  
  ego <- dotplot_enrichment(top_genes)
  if (!is.null(ego$plot)) {
    ggsave(
      file = paste0("figs/analysis/subject_l2fc_anova/", visit, "_top_anova_enrichment_dotplot.png"),
      plot = ego$plot,
      height = 6,
      width = 6,
      dpi = 300
    )
  }
  
  # # boxplots for the top 10 genes
  # for (top_gene in top_anova$genes$GENE[1:10]) {
  #   p_box <- boxplot_gene_count_with_test(
  #     l2fc_data_subset$count_matrix,
  #     l2fc_data_subset$metadata,
  #     top_gene,
  #     "DIAGNOSIS",
  #     comparisons = list(
  #       c("HC", "NSCLC"),
  #       c("HC", "MM"),
  #       c("HC", "SMM")
  #     )
  #   ) +
  #     ggtitle(paste0(
  #       top_gene,
  #       " l2fc in gene counts from V1 to ",
  #       visit
  #     )) +
  #     xlab("Diagnosis") +
  #     ylab(paste0("V1 to ", visit, " l2fc in gene counts"))
  #   
  #   ggsave(
  #     filename = paste0(
  #       "figs/analysis/subject_l2fc_anova/",
  #       visit, "_", top_gene, "_boxplot.png"
  #     ),
  #     plot = p_box,
  #     height = 6,
  #     width  = 8
  #   )
  # }
}


# print top genes
#cat(head(V3_anova_results$GENE[order(V3_anova_results$P_VALUE)], 200), sep = "\n")


#### Correlation ####



# Spearman correlation of pathway GSVA with TRI_ELISA_V4 ---------------------------
#–– 1) Define pathways of interest ––#
pathways_of_interest <- names(PATHWAYS)

#–– 2) Extract GSVA scores and metadata ––#
gsva_mat <- V1_pathway_data$gsva_matrix
meta     <- V1_pathway_data$metadata

#–– 3) Filter to valid pathways in the GSVA matrix ––#
valid_pathways   <- intersect(pathways_of_interest, rownames(gsva_mat))
missing_pathways <- setdiff(pathways_of_interest, rownames(gsva_mat))

if (length(missing_pathways) > 0) {
  warning(
    "The following pathways were not found in gsva_matrix and will be skipped:\n",
    paste(missing_pathways, collapse = ", ")
  )
}

#–– 4) Reshape GSVA to long format ––#
gsva_sub <- gsva_mat[valid_pathways, , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Pathway") %>%
  pivot_longer(
    cols      = -Pathway,
    names_to  = "CEGAT_ID",
    values_to = "GSVA"
  )

#–– 5) Merge with metadata (using existing CEGAT_ID) ––#
dat <- gsva_sub %>%
  left_join(
    meta %>% select(CEGAT_ID, TRI_ELISA_V4, DIAGNOSIS),
    by = "CEGAT_ID"
  ) %>%
  na.omit()

#–– 6) Get diagnosis groups ––#
diagnoses <- sort(unique(dat$DIAGNOSIS))

#–– 7) Compute Spearman correlations per diagnosis ––#
corr_list <- list()
p_list    <- list()

for (diag in diagnoses) {
  tmp  <- filter(dat, DIAGNOSIS == diag)
  cors <- pvs <- numeric(length(valid_pathways))
  names(cors) <- names(pvs) <- valid_pathways
  
  for (pw in valid_pathways) {
    sel <- tmp$Pathway == pw
    x   <- tmp$GSVA[sel]
    y   <- tmp$TRI_ELISA_V4[sel]
    
    if (sum(!is.na(x) & !is.na(y)) >= 2) {
      rc        <- rcorr(x, y, type = "spearman")
      cors[pw]  <- rc$r[1,2]
      pvs [pw]  <- rc$P[1,2]
    } else {
      cors[pw] <- NA
      pvs [pw] <- NA
    }
  }
  
  corr_list[[diag]] <- cors
  p_list   [[diag]] <- pvs
}

#–– 8) Compute cohort-wide Spearman correlations ––#
cors <- pvs <- numeric(length(valid_pathways))
names(cors) <- names(pvs) <- valid_pathways

for (pw in valid_pathways) {
  sel <- dat$Pathway == pw
  x   <- dat$GSVA[sel]
  y   <- dat$TRI_ELISA_V4[sel]
  
  if (sum(!is.na(x) & !is.na(y)) >= 2) {
    rc        <- rcorr(x, y, type = "spearman")
    cors[pw]  <- rc$r[1,2]
    pvs [pw]  <- rc$P[1,2]
  } else {
    cors[pw] <- NA
    pvs [pw] <- NA
  }
}

corr_list[["COHORT"]] <- cors
p_list   [["COHORT"]] <- pvs

#–– 9) Assemble correlation & p-value matrices ––#
corr_mat <- do.call(cbind, corr_list)
p_mat    <- do.call(cbind, p_list)

rownames(corr_mat) <- valid_pathways
colnames(corr_mat) <- names(corr_list)
rownames(p_mat)    <- valid_pathways
colnames(p_mat)    <- names(p_list)

#–– 10) (Optional) Filter to significant & strong effects ––#
sig_pw <- intersect(
  rownames(p_mat)[apply(p_mat, 1, function(r) any(r < 0.05, na.rm = TRUE))],
  rownames(corr_mat)[apply(corr_mat, 1, function(r) any(abs(r) > 0.2, na.rm = TRUE))]
)
if (length(sig_pw) > 0) {
  corr_mat <- corr_mat[sig_pw, , drop = FALSE]
  p_mat    <- p_mat   [sig_pw, , drop = FALSE]
}

#–– 11) Define colour ramp and draw the heatmap ––#
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

ht <- Heatmap(
  corr_mat,
  name              = "Spearman ρ",
  col               = col_fun,
  cluster_rows      = TRUE,
  cluster_columns   = TRUE,
  show_row_names    = TRUE,
  show_column_names = TRUE,
  column_title      = "GSVA Pathway vs. TRI_ELISA_V4\nby Diagnosis & Cohort",
  column_title_gp   = gpar(fontsize = 14, fontface = "bold"),
  row_names_gp      = gpar(fontsize = 10),
  column_names_gp   = gpar(fontsize = 10),
  cell_fun          = function(j, i, x, y, width, height, fill) {
    pval <- p_mat[i, j]
    if (!is.na(pval) && pval < 0.05) {
      stars <- if      (pval < 0.001) "***"
      else if (pval < 0.01 ) "**"
      else                    "*"
      grid.text(stars, x, y, gp = gpar(fontsize = 12))
    }
  }
)

#–– 12) Save as PNG ––#
png(
  filename = "figs/correlation/pathway_gsva_spearman_corr_heatmap.png",
  width    = 1800,
  height   = 2000,
  res      = 300
)
draw(ht, padding = unit(c(10, 40, 10, 10), "mm"))
dev.off()



