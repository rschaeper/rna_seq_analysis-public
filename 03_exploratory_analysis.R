# script to uncover structured parts of data and whether there could be a link to metadata variables

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
# library(Biobase)
# library(plyr)
# library(foreach)
# library(xtable)
# library(biomaRt)
# library(GOstats)
# library(cluster)
# library(marray)
# library(mclust)
# library(RColorBrewer)
# library(igraph)
# library(Rgraphviz)
# library(graph)
# library(colorspace)
# library(annotate)
# library(scales)
# library(gtools)
# 
# library(MineICA)

# prepare data ---------------------------------------------------------------

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

# ideally set to FALSE, but changes won't be drastic and takes much longer
V1_vsd <- as.data.frame(assay(vst(V1_dds, blind = TRUE)))

# load l2fc data
l2fc_data_imputed <- readRDS("output/l2fc_data_imputed_qc.rds")

# remove V1 samples from l2fc count data for this script
l2fc_data_imputed <- subset_l2fc_data(l2fc_data_imputed, condition = "VISIT != 'V1'")

# load euclidean distances for dispersion study
euclid_dists <- readRDS("output/euclid_dists_qc.rds")

V1_euclid_dists <- subset_euclid_dists(euclid_dists, condition = "VISIT == 'V1'")

# visualize TRI -----------------------------------------------------------

boxplot_with_cohort_and_quartiles <- function(data, x_col, y_col, comparisons = NULL,
                                              test_method = "t.test", show_quartiles = TRUE) {
  df <- data %>%
    mutate(!!x_col := as.character(.data[[x_col]])) %>%
    bind_rows(
      transmute(data,
                !!x_col := "COHORT",
                !!y_col := .data[[y_col]])
    ) %>%
    mutate(!!x_col := factor(.data[[x_col]], levels = c("HC", "MM", "NSCLC", "SMM", "COHORT")))
  
  if (show_quartiles) {
    # compute Q1 and Q3
    quantiles <- df %>%
      group_by(across(all_of(x_col))) %>%
      summarise(
        q25 = quantile(.data[[y_col]], 0.25, na.rm = TRUE),
        q75 = quantile(.data[[y_col]], 0.75, na.rm = TRUE),
        .groups = "drop"
      )
    
    # join and annotate
    df <- df %>%
      left_join(quantiles, by = x_col) %>%
      mutate(
        is_upper = .data[[y_col]] >= q75,
        is_lower = .data[[y_col]] <= q25,
        quartile_flag = case_when(
          is_upper ~ "upper",
          is_lower ~ "lower",
          TRUE     ~ "middle"
        ),
        quartile_flag = factor(quartile_flag, levels = c("upper", "middle", "lower"))
      )
  }
  
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_boxplot(aes(color = .data[[x_col]]), width = 0.8, outlier.shape = NA)
  
  if (show_quartiles) {
    p <- p +
      geom_jitter(aes(color = .data[[x_col]], fill = quartile_flag),
                  shape = 21, width = 0.15, size = 2, stroke = 0.6, alpha = 0.8) +
      scale_fill_manual(
        name = "Quartile",
        values = c("upper" = "red", "middle" = "white", "lower" = "blue"),
        labels = c("Upper quartile", "Middle 50%", "Lower quartile")
      )
  } else {
    p <- p +
      geom_jitter(aes(color = .data[[x_col]]),
                  shape = 21, fill = NA, width = 0.15, size = 2, stroke = 0.6, alpha = 0.8)
  }
  
  p <- p +
    { if (!is.null(comparisons)) {
      stat_compare_means(comparisons = comparisons, method = test_method, label = "p.signif")
    } else {
      stat_compare_means(method = "anova", label.x = 1.5)
    }
    } +
    scale_color_manual(values = c(DIAGNOSIS_COLKEY, COHORT = "grey")) +
    theme_bw() +
    theme(
      legend.position = if (show_quartiles) "right" else "none",
      plot.title      = element_text(hjust = 0.5),
      axis.text.x     = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}


for (ab_var in c("TRI_ELISA_V3", "TRI_ELISA_V4", "TRI_ELISA_V5",
                 "IgG_V1",       "IgG_V3",       "IgG_V4")) {
  
  # pull out the "1", "3", "4", or "5" after "_V"
  v_num <- sub(".*_V", "", ab_var)
  
  # inline map to d0/d7/d30/d90
  day_label <- if (v_num == "1") {
    "d0"
  } else if (v_num == "3") {
    "d7"
  } else if (v_num == "4") {
    "d30"
  } else if (v_num == "5") {
    "d90"
  } else {
    stop("Unknown V number: ", v_num)
  }
  
  # pick the assay‑specific part and build full y‑label
  y_lab <- paste0(
    day_label, " ",
    if (grepl("^IgG", ab_var)) {
      "influenza binding IgG titer"
    } else {
      "titer response index"
    }
  )
  
  # make the boxplot
  boxp <- boxplot_with_cohort_and_quartiles(
    data           = V1_metadata,
    x_col          = "DIAGNOSIS",
    y_col          = ab_var,
    comparisons    = list(c("HC","SMM"), c("HC","NSCLC"), c("HC","MM")),
    test_method    = "wilcox.test",
    show_quartiles = FALSE
  ) +
    labs(
      title = ab_var,
      x     = NULL,
      y     = y_lab
    )
  
  # save it
  ggsave(
    filename = sprintf(
      "figs/exploration/COHORT_%s_%s_colby_DIAGNOSIS.png",
      day_label, ab_var
    ),
    plot   = boxp,
    width  = 6,
    height = 6,
    dpi    = 300
  )
}

# visualize AGE -----------------------------------------------------------

# get the quartiles
quartiles_age <- quantile(V1_metadata$AGE, probs = c(0.25, 0.75))

# extract the lower bound (Q1)
lower_bound_age <- quartiles_age[1]

# extract the upper bound (Q3)
upper_bound_age <- quartiles_age[2]

message("The interquartile range of age for the cohort is:")
print(paste("Lower Bound (Q1):", lower_bound_age))
print(paste("Upper Bound (Q3):", upper_bound_age))

boxp <- boxplot_with_test(V1_metadata,
                          "DIAGNOSIS",
                          "AGE",
                          comparisons = list(c("HC", "SMM"), c("HC", "NSCLC"), c("HC", "MM")),
                          test_method = "t.test") +
  labs(
    title = "Age",
    x     = "",
    y     = "Age in years"
  ) +
  scale_color_manual(values = DIAGNOSIS_COLKEY)
ggsave(paste0("figs/exploration/COHORT_V1_AGE_colby_DIAGNOSIS.png"), plot = boxp, width = 6, height = 6, dpi = 300)



# normalized gene count sample dispersion ------------------------------------

# compute pairwise‐SD within each group
group <- V1_euclid_dists$metadata$DIAGNOSIS
pairwise_sd <- sapply(levels(group), function(g) {
  idx  <- which(group == g)
  subm <- V1_euclid_dists$matrix[idx, idx]               # sub‑matrix for group g
  sd(subm[lower.tri(subm)])         # SD of all i<j distances (lower triangle, no duplicates or self-pairs)
})
print(pairwise_sd)

# boxplot of all within‐group pairwise distances
# (this gives the full distribution, not just the single SD summary)

# melt distance matrix to long form
dist_df <- melt(V1_euclid_dists$matrix, varnames = c("S1","S2"), value.name="Distance")
# keep each pair once, and only within‐group pairs
dist_df <- subset(dist_df, as.character(S1) < as.character(S2))
dist_df$G1 <- V1_euclid_dists$metadata[as.character(dist_df$S1),"DIAGNOSIS"]
dist_df$G2 <- V1_euclid_dists$metadata[as.character(dist_df$S2),"DIAGNOSIS"]
within_pairs <- subset(dist_df, G1 == G2)

ggplot(within_pairs, aes(x = G1, y = Distance, colour = G1)) +
  geom_boxplot() +
  geom_jitter(width=0.3, alpha=0.5) +
  theme_bw() +
  labs(
    x = "Diagnosis",
    y = "Pairwise euclidean distance",
    title = "Distribution of within‑group pairwise distances for V1 samples"
  ) +
  scale_colour_manual(
    name = "Diagnosis",
    values = DIAGNOSIS_COLKEY) # + stat_compare_means(method = "kruskal.test", label.y = max(within_pairs$Distance) * 1.05)
  

ggsave("figs/exploration/within-group_distances_boxplot.png", width = 8, height = 6)

# test whether those distributions differ across diagnoses
kruskal_res <- kruskal.test(Distance ~ G1, data = within_pairs)
print(kruskal_res)



# keep each pair once, and only cross‐group pairs
cross_pairs <- subset(dist_df, 
                      as.character(S1) < as.character(S2))

# create a unified 'Pair' factor (e.g. "CTRL–AD" rather than "AD–CTRL" and "CTRL–AD")
cross_pairs$Pair <- with(cross_pairs, {
  g1_g2 <- paste(G1, G2, sep = "–")
  g2_g1 <- paste(G2, G1, sep = "–")
  # if g1_g2 is lexicographically ≤ g2_g1, keep that; otherwise swap
  ifelse(g1_g2 <= g2_g1, g1_g2, g2_g1)
})
cross_pairs$Pair <- factor(cross_pairs$Pair,
                           levels = sort(unique(cross_pairs$Pair)))

# ─────────────────────────────────────────────────────────────────────────────
# compute mean±SD for each diagnosis-pair (fixed)
# ─────────────────────────────────────────────────────────────────────────────

cross_stats <- cross_pairs %>%
  group_by(Pair) %>%
  summarise(
    MeanDist = mean(Distance),
    SdDist   = sd(Distance),
    Npairs   = length(Distance),    # <-- replace n() with length(Distance)
    .groups  = "drop"
  )

print(cross_stats)

# ─────────────────────────────────────────────────────────────────────────────
# visualize between‐group distances
# ─────────────────────────────────────────────────────────────────────────────

# boxplot + jitter for all between‐group pairs
ggplot(cross_pairs, aes(x = Pair, y = Distance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  theme_bw(base_size = 14) +
  labs(
    x = "Diagnosis Pair",
    y = "Pairwise Euclidean Distance",
    title = "Between‐Group Pairwise Distance Distributions"
  ) +
  #stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") + # shows the mean as well
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) # + stat_compare_means(method = "kruskal.test", label.y = max(cross_pairs$Distance) * 1.05)
  

# save figure
ggsave("figs/exploration/between-group_distances_boxplot.png",
       width = 8, height = 5, dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# statistical tests (e.g. compare each pair to the overall distribution)
# ─────────────────────────────────────────────────────────────────────────────

# Kruskal–Wallis test across all cross‐group pairs
kruskal_between <- kruskal.test(Distance ~ Pair, data = cross_pairs)
print(kruskal_between)

# Post‐hoc Dunn’s test (requires FSA or dunn.test package)
library(FSA)
dunn_between <- dunnTest(Distance ~ Pair, data = cross_pairs, method = "bh")
print(dunn_between)

# ─────────────────────────────────────────────────────────────────────────────
# plot heatmap
# ─────────────────────────────────────────────────────────────────────────────

# get your diagnosis levels
groups <- sort(unique(V1_euclid_dists$metadata$DIAGNOSIS))

# melt distance matrix (you already did this)
# dist_df has columns S1, S2, Distance, G1, G2, and contains only S1 < S2 pairs

# initialize an empty matrix
mat <- matrix(NA,
              nrow = length(groups),
              ncol = length(groups),
              dimnames = list(groups, groups))

# fill in medians for each pair of groups
for (i in seq_along(groups)) {
  for (j in seq_along(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]
    if (g1 == g2) {
      # within‐group: filter where G1==G2==g1
      vals <- dist_df$Distance[ dist_df$G1 == g1 & dist_df$G2 == g1 ]
    } else {
      # between‐group: filter both orientations
      vals <- dist_df$Distance[
        (dist_df$G1 == g1 & dist_df$G2 == g2) |
          (dist_df$G1 == g2 & dist_df$G2 == g1)
      ]
    }
    mat[g1, g2] <- median(vals, na.rm = TRUE)
  }
}

# slightly cheeky variant:
col_fun <- colorRamp2(c(min(mat), max(mat)), c("white", "red"))

# now plot with pheatmap
png("figs/exploration/median_distances_heatmap.png", width = 2400, height = 2400, res = 300)
Heatmap(
  mat,
  name = "Median\nEuclidean\ndistance",
  column_title = "Median of pairwise Euclidean distances",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 18),
  column_names_gp = gpar(fontsize = 18),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.0f", mat[i, j]),
      x = x,
      y = y,
      gp = gpar(fontsize = 18)
    )
  }
)
dev.off()


# subject log2FoldChange gene count sample dispersion --------------

# ─────────────────────────────────────────────────────────────────────────────
# 1) Compute or load Euclidean distances between samples (log2FC)
# ─────────────────────────────────────────────────────────────────────────────
l2fc_euclid_dists_file <- "output/l2fc_euclid_dists.rds"
if (file.exists(l2fc_euclid_dists_file)) {
  euclid_dists <- readRDS(l2fc_euclid_dists_file)
  message("Loaded existing sample l2fc euclidean distances from file.")
} else {
  message("Calculating l2fc euclidean distances between samples")
  euclid_dists <- list()
  # dist() expects rows = samples, so transpose count_matrix
  euclid_dists$matrix <- as.matrix(dist(t(l2fc_data_imputed$count_matrix), method = "euclidean"))
  euclid_dists$metadata <- l2fc_data_imputed$metadata
  saveRDS(euclid_dists, file = l2fc_euclid_dists_file)
}

# ─────────────────────────────────────────────────────────────────────────────
# 2) Boxplot: Within‐group pairwise distances across visits, faceted by diagnosis
# ─────────────────────────────────────────────────────────────────────────────

# Helper function to subset euclid_dists by a condition on metadata.
# Assumes euclid_dists is a list with:
#   $matrix   : sample × sample numeric matrix
#   $metadata : data.frame with rownames = sample IDs; has columns including "VISIT" and "DIAGNOSIS"
subset_euclid_dists <- function(euclid_dists, condition) {
  # 'condition' is a string to pass to subset(), e.g. "VISIT == 'V2'"
  md <- euclid_dists$metadata
  keep_samples <- rownames(md)[eval(parse(text = condition), envir = md)]
  sub_mat <- euclid_dists$matrix[keep_samples, keep_samples, drop = FALSE]
  sub_md  <- md[keep_samples, , drop = FALSE]
  return(list(matrix = sub_mat, metadata = sub_md))
}

# 2.1) Build a single data‐frame of within‐group distances for all visits (excluding baseline if VISITS[1])
all_within <- lapply(VISITS[2:4], function(visit) {
  # subset by this visit
  ed <- subset_euclid_dists(euclid_dists, condition = paste0("VISIT == '", visit, "'"))
  
  # melt the sub‐matrix into long format
  dist_df <- melt(ed$matrix, varnames = c("S1","S2"), value.name = "Distance")
  dist_df <- dist_df %>%
    filter(as.character(S1) < as.character(S2)) %>%             # keep unique pairs
    mutate(
      G1    = ed$metadata[as.character(S1), "DIAGNOSIS"],
      G2    = ed$metadata[as.character(S2), "DIAGNOSIS"],
      Visit = visit
    ) %>%
    filter(G1 == G2) %>%                                         # only within‐diagnosis
    dplyr::select(Visit, Distance, Diagnosis = G1)
  
  return(dist_df)
}) %>%
  bind_rows() %>%
  mutate(
    Visit     = factor(Visit, levels = VISITS[2:3]),
    Diagnosis = factor(Diagnosis, levels = names(DIAGNOSIS_COLKEY))
  )

# 2.2) Plot boxplot + jitter
p <- ggplot(all_within, aes(x = Visit, y = Distance, fill = Diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour = Diagnosis), width = 0.2, size = 0.8, alpha = 0.4) +
  facet_grid(Diagnosis ~ ., scales = "fixed", switch = "y") +
  scale_fill_manual(values = DIAGNOSIS_COLKEY) +
  scale_colour_manual(values = DIAGNOSIS_COLKEY) +
  theme_bw() +
  theme(
    strip.background  = element_blank(),
    strip.placement   = "outside",
    panel.spacing.y   = unit(0.5, "lines"),
    axis.text.x       = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x     = "Visit",
    y     = "Within‐group pairwise Euclidean distance",
    title = "Log2FoldChange dispersion across visits, by diagnosis"
  )

ggsave(
  filename = "figs/exploration/all_visits_within-group_distances.png",
  plot    = p,
  width   = 6,
  height  = 8,
  dpi     = 300
)

# ─────────────────────────────────────────────────────────────────────────────
# 3) Complex Heatmap: Median pairwise distances for every (Visit × Diagnosis)
# ─────────────────────────────────────────────────────────────────────────────

# 3.1) Melt the full euclidean distance matrix into a long data.frame
dist_df_all <- melt(
  euclid_dists$matrix,
  varnames   = c("S1", "S2"),
  value.name = "Distance"
) %>%
  filter(as.character(S1) < as.character(S2)) %>%
  mutate(
    Visit1  = euclid_dists$metadata[as.character(S1), "VISIT"],
    Diag1   = euclid_dists$metadata[as.character(S1), "DIAGNOSIS"],
    Visit2  = euclid_dists$metadata[as.character(S2), "VISIT"],
    Diag2   = euclid_dists$metadata[as.character(S2), "DIAGNOSIS"],
    Group1  = paste(Visit1, Diag1, sep = "_"),
    Group2  = paste(Visit2, Diag2, sep = "_")
  )

# 3.2) Enumerate all (Visit × Diagnosis) combos in desired order
visits    <- VISITS[2:3]                     # exclude baseline if needed
diagnoses <- names(DIAGNOSIS_COLKEY)        # diagnosis levels
all_groups <- as.vector(outer(visits, diagnoses, paste, sep = "_"))

# 3.3) Initialize empty square matrix of medians
n_grp      <- length(all_groups)
median_mat <- matrix(NA_real_,
                     nrow = n_grp,
                     ncol = n_grp,
                     dimnames = list(all_groups, all_groups))

# 3.4) Fill in medians for each pair of groups
for (i in seq_len(n_grp)) {
  for (j in seq_len(n_grp)) {
    g1 <- all_groups[i]
    g2 <- all_groups[j]
    if (g1 == g2) {
      # within‐group
      vals <- dist_df_all$Distance[
        dist_df_all$Group1 == g1 & dist_df_all$Group2 == g1
      ]
    } else {
      # between‐group (both orientations)
      vals <- dist_df_all$Distance[
        (dist_df_all$Group1 == g1 & dist_df_all$Group2 == g2) |
          (dist_df_all$Group1 == g2 & dist_df_all$Group2 == g1)
      ]
    }
    median_mat[i, j] <- median(vals, na.rm = TRUE)
  }
}

# 3.5) Define color ramp from min → max of median_mat
col_fun <- colorRamp2(
  c(min(median_mat, na.rm = TRUE),
    max(median_mat, na.rm = TRUE)),
  c("white", "red")# c("#F7F7F7", "#B2182B")
)

# 3.6) Build a small annotation: split each “Visit_Diagnosis” back out
split_df <- data.frame(
  Group     = all_groups,
  Visit     = rep(visits, each = length(diagnoses)),
  Diagnosis = rep(diagnoses, times = length(visits)),
  stringsAsFactors = FALSE
)

row_anno <- rowAnnotation(
  Diagnosis           = split_df$Diagnosis,
  col                 = list(Diagnosis = DIAGNOSIS_COLKEY),
  show_annotation_name = FALSE,
  width               = unit(4, "mm")
)

# 3.7) Draw and save the ComplexHeatmap
png(
  filename = "figs/exploration/median_pairwise_distances_visit_diag_heatmap.png",
  width    = 3000,
  height   = 3000,
  res      = 300
)
Heatmap(
  median_mat,
  name               = "Median\nEuclidean\ndistance",
  col                = col_fun,
  cluster_rows       = FALSE,
  cluster_columns    = FALSE,
  row_names_gp       = gpar(fontsize = 10),
  column_names_gp    = gpar(fontsize = 10),
  #left_annotation    = row_anno,
  column_title       = "Median pairwise Euclidean distance\n(Visit × Diagnosis vs Visit × Diagnosis)",
  row_title          = NULL,
  show_row_dend      = FALSE,
  show_column_dend   = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.0f", median_mat[i, j]),
      x     = x,
      y     = y,
      gp    = gpar(fontsize = 12)
    )
  }
)
dev.off()

# COHORT exploration ------------------------------------------------------
##### 1. correlate metadata variables #####
form <- ~ STUDY + SEX + AGE + DIAGNOSIS + IgG_V1 + IgG_V3 + IgG_V4 +
  INFLUENZA_VACCINATIONS_LAST_5_YEARS + COVID19_VACCINATIONS +
  HISTORY_OF_CANCER + REGULAR_MEDICATION

# compute CCA
C <- canCorPairs(form, V1_metadata)

# adjust rownames too long
old_rows <- rownames(C)
new_rows <- ifelse(
  old_rows == "INFLUENZA_VACCINATIONS_LAST_5_YEARS",
  "INFLUENZA_VACCINATIONS_\nLAST_5_YEARS",
  old_rows
)
rownames(C) <- new_rows

# adjust columnnames too long
old_columns <- colnames(C)
new_columns <- ifelse(
  old_columns == "INFLUENZA_VACCINATIONS_LAST_5_YEARS",
  "INFLUENZA_VACCINATIONS_\nLAST_5_YEARS",
  old_columns
)
colnames(C) <- new_columns


# define a diverging color scale
col_fun <- colorRamp2(c(0, 1), c("white", "red"))

# open graphics device if saving to file
png("figs/exploration/COHORT_cca_heatmap.png",
    width  = 2400,   # pixels
    height = 2400,
    res    = 300)

# build the heatmap object
ht <- Heatmap(
  C,
  name                   = "Correlation",    # legend title (must be character)
  col                    = col_fun,
  cluster_rows           = TRUE,
  cluster_columns        = TRUE,
  show_row_dend          = TRUE,
  show_column_dend       = TRUE,
  # control dendrogram space
  row_dend_width         = unit(2, "cm"),
  row_dend_side          = "right",
  column_dend_height     = unit(2, "cm"),
  # label styling
  row_names_side         = "left", 
  row_names_gp           = gpar(fontsize = 12),
  column_names_gp        = gpar(fontsize = 12),
  column_names_rot       = 50,
  # cell size
  width                  = unit(10, "cm"),
  height                 = unit(10, "cm"),
  # add a heatmap title
  column_title           = "Cohort CCA Heatmap",
  column_title_gp        = gpar(fontsize = 14, fontface = "bold"),
  column_title_side      = "top"
)

# draw the heatmap (titles & legends will render correctly)
draw(
  ht,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right"
)

dev.off()

##### 2. normalized gene count variance partitioning analysis #####

# Let us investigate how well the metadata variables we assume to explain a lot of differences in our data actually perform

# We assume inter-subject variability to be high and thus that STUDY_ID can explain most of the variance
# Residual variance will likely be technical and due to the temporal axis (i.e. vaccine induced or subject's life circumstance's effects)

# let us model STUDY_ID + VISIT
form <- ~ (1 | STUDY_ID) + (1 | VISIT)

# the varPart model enforces us to model either all or none of the categorical variables as random effects
# it is recommended by the developers to model categorical variables with many levels, like inidividual (STUDY_ID) as random effects

varPart <- fitExtractVarPartModel(vsd, form, metadata)

# create violin plot of contribution of each variable to total variance
plotVarPart(varPart) + ggtitle("Variance explained by subject and visit")

# save the plot
ggsave("figs/exploration/normalized_counts_STUDY_ID_VISIT_variance_violins.png", width = 6, height = 4, dpi = 300)

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
    enriched_plot <- ego$plot + ggtitle(paste0("COHORT - BTM enrichment for top genes with most variance explained by ", var, " in normalized count space"))
    ggsave(
      filename = paste0("figs/exploration/normalized_counts_variance_violins_", var, "_top_genes_enrichment.png"),
      plot = enriched_plot,
      height = 8,
      width = 10
    )
  } else {
    message(paste("No enrichment found for var", var))
  }
  
}

# attempt to break down variance explained by STUDY_ID
form <- ~ (1 | SEX) + (1 | DIAGNOSIS) + AGE + (1 | VISIT)

varPart <- fitExtractVarPartModel(vsd, form, metadata)

# create violin plot of contribution of each variable to total variance
plotVarPart(varPart) + ggtitle("Variance explained by diagnosis, sex, visit and age")

# save the plot
ggsave("figs/exploration/normalized_counts_SEX_DIAGNOSIS_AGE_VISIT_variance_violins.png", width = 6, height = 4, dpi = 300)

message("Completed exploratory analysis.")
