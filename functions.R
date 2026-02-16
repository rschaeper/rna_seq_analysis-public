# Define global analysis parameters (move this elsewhere at some point)

# fgsea uses Monte Carlo Simulations and is inherently stochastic, seed the code
set.seed(42)

# CONSTANTS ---------------------------------------------------------------

# define DIAGNOSIS via constant list of strings
DIAGNOSES <- c("HC", "SMM", "MM", "NSCLC")

# define VISITS to contrast as constant list of strings
VISITS <- c("V1", "V2", "V3", "V4", "V5", "V6")

# whether to calculate 4 ANOVA models in 05_kinetics.R (takes ~ 20 min)
COMPARE_ANOVA_MODELS = FALSE

HALLMARK_FILE <- "data/gene_sets/h.all.v2024.1.Hs.symbols.gmt"

# define pathways for gene set enrichment analysis
BTM_FILE <- "data/gene_sets/BTM_for_GSEA_20131008.gmt"

REACTOME_FILE <- "data/gene_sets/c2.cp.reactome.v2024.1.Hs.symbols.gmt"

KEGG_FILE <- "data/gene_sets/c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt"

C7_IMMUNO_FILE <- "data/gene_sets/c7.all.v2024.1.Hs.symbols.gmt"

VAX_FILE <- "data/gene_sets/c7.vax.v2024.1.Hs.symbols.gmt"

# load pathways
library(fgsea)
PATHWAYS <- fgsea::gmtPathways(BTM_FILE)

# define diagnosis colors
DIAGNOSIS_COLKEY <- c(
  "HC" = "#FFC700",
  "SMM" = "#00C7E1",
  "MM" = "#00559E",
  "NSCLC" = "#FE1F04"
)

# desaturated version of colors can be created in code with colorspace package

# define visit colors
VISIT_COLKEY <- c(V1 = "#e41a1c", V2 = "#377eb8", V3 = "#4daf4a", V4 = "#984ea3", V5 = "#ff7f00", V6 = "#ffff33")

# functions ---------------------------------------------------------------

# function to easily subset dds

subset_dds <- function(dds, condition = "STUDY == 'SYS01' && (DIAGNOSIS == ('HC' || 'MM'))") {
  # evaluate the condition string
  keep <- eval(parse(text = condition), envir = colData(dds))
  
  # subset the dds data using the 'keep' logical vector
  dds <- dds[, keep]
  
  # make sure to reset factor levels
  colData(dds) <- droplevels(colData(dds))
  
  return(dds)
}

# example usage:

# dds <- subset_dds(dds, condition = "STUDY == 'SYS01'")


subset_l2fc_data <- function(l2fc_data, condition = "STUDY == 'SYS01' && (DIAGNOSIS == ('HC' || 'MM'))") {
  
  l2fc_metadata <- l2fc_data$metadata
  l2fc_count_matrix <- l2fc_data$count_matrix
  
  keep <- eval(parse(text = condition), envir = l2fc_metadata)
  
  l2fc_metadata_filtered <- l2fc_metadata[keep, ]
  
  # relevel metadata
  l2fc_metadata_filtered <- droplevels(l2fc_metadata_filtered)
  
  # reorder the columns of our count matrix to follow the order of our metadata
  idx <- match(l2fc_metadata_filtered$CEGAT_ID, colnames(l2fc_count_matrix))
  l2fc_count_matrix_filtered <- l2fc_count_matrix[, idx]
  
  l2fc_data$metadata <- l2fc_metadata_filtered
  l2fc_data$count_matrix <- l2fc_count_matrix_filtered
  
  return(l2fc_data)
}

# example usage:

# l2fc_data_imputed <- subset_l2fc_data(l2fc_data_imputed, condition = "STUDY == 'SYS01'")


subset_euclid_dists <- function(euclid_dists,
                                condition = "STUDY == 'SYS01' & DIAGNOSIS %in% c('HC','MM')") {
  
  # pull apart
  md   <- euclid_dists$metadata
  mat  <- euclid_dists$matrix
  
  # evaluate your condition inside the metadata
  keep <- with(md, eval(parse(text = condition)))
  
  # subset metadata
  md_f <- droplevels(md[keep, ])
  
  # get matching indices in the matrix
  idx <- match(md_f$CEGAT_ID, colnames(mat))
  
  # subset *both* rows and columns of the (symmetric) distance matrix
  mat_f <- mat[idx, idx, drop = FALSE]
  
  list(
    metadata = md_f,
    matrix   = mat_f
  )
}

subset_pathway_data <- function(pathway_data, condition = "STUDY == 'SYS01' && (DIAGNOSIS == ('HC' || 'MM'))") {
  
  pathway_metadata <- pathway_data$metadata
  pathway_gsva_matrix <- pathway_data$gsva_matrix
  
  keep <- eval(parse(text = condition), envir = pathway_metadata)
  
  pathway_metadata_filtered <- pathway_metadata[keep, ]
  
  # relevel metadata
  pathway_metadata_filtered <- droplevels(pathway_metadata_filtered)
  
  # reorder the columns of our pathway matrix to follow the order of our metadata
  idx <- match(pathway_metadata_filtered$CEGAT_ID, colnames(pathway_gsva_matrix))
  pathway_gsva_matrix_filtered <- pathway_gsva_matrix[, idx]
  
  pathway_data$metadata <- pathway_metadata_filtered
  pathway_data$gsva_matrix <- pathway_gsva_matrix_filtered
  
  return(pathway_data)
}

# simple wrapper for enrichment using GO or custom GMT
# if gmt_file is NULL, perform GO enrichment; otherwise use GMT for custom gene sets
dotplot_enrichment <- function(genes, gmt_file = BTM_FILE, organism_db = org.Hs.eg.db,
                               keyType = "SYMBOL", pAdjustMethod = "BH",
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                               showCategory = 10) {
  library(clusterProfiler)
  
  if (!is.null(gmt_file)) {
    # read GMT file into term2gene data.frame
    # uses clusterProfiler::read.gmt
    term2gene_df <- read.gmt(gmt_file)
    # term2gene_df has columns: term, gene
    # run enricher with custom gene sets
    enr <- enricher(
      gene = genes,
      TERM2GENE = term2gene_df,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff
    )
    res_df <- as.data.frame(enr)
    if (nrow(res_df) == 0) {
      warning("No enriched terms found in GMT for provided gene list.")
      return(list(plot = NULL, data = enr))
    }
    dp <- dotplot(enr, showCategory = showCategory)
    return(list(plot = dp, data = enr))
  } else {
    # GO enrichment
    ego <- enrichGO(
      gene          = genes,
      OrgDb         = organism_db,
      keyType       = keyType,
      ont           = "BP",
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff  = pvalueCutoff,
      qvalueCutoff  = qvalueCutoff
    )
    ego_df <- as.data.frame(ego)
    if (nrow(ego_df) == 0) {
      warning("No enriched GO terms found for the provided gene list.")
      return(list(plot = NULL, data = ego))
    }
    dp <- dotplot(ego, showCategory = showCategory)
    return(list(plot = dp, data = ego))
  }
}

# function to boxplot values from two columns of dataframe
# make sure to use wilcoxon for e.g. log2FoldChange in Ab levels - this is not normally distributed!
boxplot_with_test <- function(data, x_col, y_col, comparisons = NULL, test_method = "t.test") {
  # creates a boxplot with optional pairwise comparisons.
  # 
  # args:
  #   data: The data frame.
  #   x_col: The column name for the x-axis (grouping).
  #   y_col: The column name for the y-axis (values).
  #   comparisons: A list of comparisons (pairs of group names) or NULL for no comparisons.
  
  # generate boxplot
  p <- ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[x_col]])) +
    geom_boxplot(width = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.15, aes(color = factor(.data[[x_col]]))) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title   = element_text(hjust = 0.5),
      axis.text.x  = element_text(angle = 45, hjust = 1)
    )
  if (!is.null(comparisons)) {
    # add comparisons if they are provided
    p <- p + stat_compare_means(comparisons = comparisons,
                                method = test_method,
                                label = "p.signif")
  } else {
    # consider giving the option for Kruskal-Wallis as non-parametric alt. to ANOVA
    p <- p + stat_compare_means(method = "anova", label.x = 1.5)
  }
  
  return(p)
}

# create a boxplot for normalized counts and group by metadata
boxplot_gene_count_with_test <- function(count_matrix, metadata, gene, group_col, comparisons = NULL, test_method = "t.test") {
  library(stringr)
  
  # check if the gene exists in the count matrix
  if (!gene %in% rownames(count_matrix)) {
    stop(paste("Gene", gene, "not found in count_matrix."))
  }
  
  # create a data frame of gene counts for the specified gene.
  df_counts <- data.frame(
    CEGAT_ID = colnames(count_matrix),
    COUNT = as.numeric(count_matrix[gene, ])
  )
  
  # merge the counts with metadata on sample_id
  df <- merge(df_counts, metadata, by = "CEGAT_ID")
  
  # generate the boxplot using the grouping specified by group_col
  p <- ggplot(df, aes(x = .data[[group_col]], y = COUNT, color = .data[[group_col]])) +
    geom_boxplot(width = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.15) + # color is already in the main aes
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title   = element_text(hjust = 0.5),
      axis.text.x  = element_text(angle = 45, hjust = 1)
    ) +
    labs(title = paste(gene, "Gene Count"), fill = str_to_title(group_col))
  
  # add statistical comparisons
  if (!is.null(comparisons)) {
    p <- p + stat_compare_means(comparisons = comparisons,
                                method = test_method,
                                label = "p.signif")
  } else {
    # consider giving the option for Kruskal-Wallis as non-parametric alt. to ANOVA
    # read "compare means" as “statistically compare central tendencies”
    p <- p + stat_compare_means(method = "kruskal.test", label.x = 0.8, label.y = max(df$COUNT), show.legend = FALSE)
  }
  
  return(p)
}

# example usage:
#boxplot_gene_count_with_test(count_matrix, metadata, "MZB1", "DIAGNOSIS", comparisons = list(c("HC", "SMM"), c("HC", "NSCLC"), c("HC", "MM")))




