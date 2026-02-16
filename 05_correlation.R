# script to correlate transcriptomic signatures with antibodies

library(dplyr)
library(tidyr)
library(PCAtools)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(DESeq2)
library(Hmisc)
library(fgsea)

plot_gsea_dotmatrix <- function(gsea_result,
                                conditions,
                                padj_prefix      = "padj.",
                                nes_prefix       = "NES.",
                                padj_threshold   = 0.001,
                                selected_pathways = NULL,
                                cluster           = FALSE,
                                pathway2genes     = NULL,
                                cluster_method    = "average",
                                cluster_height    = 0.5) {
  
  padj_cols <- paste0(padj_prefix, conditions)
  nes_cols  <- paste0(nes_prefix,  conditions)
  
  # 1) FILTER or SELECT
  if (!is.null(selected_pathways)) {
    df <- gsea_result %>% filter(PATHWAY %in% selected_pathways)
    missing <- setdiff(selected_pathways, df$PATHWAY)
    if (length(missing))
      warning("These selected_pathways were not found: ",
              paste(missing, collapse = ", "))
  } else {
    df <- gsea_result %>%
      filter_at(vars(all_of(padj_cols)), any_vars(. < padj_threshold)) %>%
      filter(!grepl("TBA", PATHWAY, ignore.case = TRUE))
  }
  
  if (nrow(df) < 2) {
    message("  Skipping dotmatrix: Not enough pathways (", nrow(df), ") after initial filtering.")
    return(invisible(NULL))
  }
  
  # 2) OPTIONAL JACCARD CLUSTERING + COLLAPSE
  if (cluster) {
    if (is.null(pathway2genes))
      stop("To cluster, supply `pathway2genes` (named list of gene vectors).")
    
    # build binary matrix genes x pathways
    pgs      <- pathway2genes[df$PATHWAY]
    all_genes<- unique(unlist(pgs))
    mat      <- sapply(pgs, function(g) all_genes %in% g)
    rownames(mat) <- all_genes
    
    if (ncol(mat) < 2 || nrow(mat) == 0) { # ncol(mat) corresponds to number of pathways
      message("  Skipping Jaccard clustering: Not enough pathways (", ncol(mat), ") or genes (", nrow(mat), ").")
      # You might want to handle this by setting cluster to FALSE and continuing without Jaccard clustering
      # For now, let's make it return NULL, as the primary intent is clustering.
      return(invisible(NULL))
    }
    
    # Jaccard similarity
    inter  <- crossprod(mat)
    sizes  <- colSums(mat)
    union  <- outer(sizes, sizes, "+") - inter
    sim    <- inter / union
    dist_j <- as.dist(1 - sim)
    hc     <- hclust(dist_j, method = cluster_method)
    
    # assign clusters
    cls    <- cutree(hc, h = cluster_height)
    
    # compute each pathway's min padj
    min_p <- df %>%
      rowwise() %>%
      mutate(min_padj = min(c_across(all_of(padj_cols)), na.rm = TRUE)) %>%
      ungroup() %>%
      select(PATHWAY, min_padj)
    
    # pick one rep per cluster
    rep_paths <- tibble(PATHWAY = names(cls), cluster = cls) %>%
      inner_join(min_p, by = "PATHWAY") %>%
      group_by(cluster) %>%
      slice_min(min_padj, with_ties = FALSE) %>%
      pull(PATHWAY)
    
    # subset df to reps only
    df <- df %>% filter(PATHWAY %in% rep_paths)
    # reorder by dendrogram order of clusters
    ordered_labels <- hc$labels[hc$order]
    # keep only those that survived rep selection
    row_order <- intersect(ordered_labels, rep_paths)
    if (length(row_order) < 2) {
      message("  Skipping dotmatrix: Not enough pathways (", length(row_order), ") after Jaccard clustering and representative selection.")
      return(invisible(NULL))
    }
  } else {
    row_order <- unique(df$PATHWAY)
  }
  
  # 3) PIVOT LONG
  df_nes  <- df %>%
    tidyr::pivot_longer(all_of(nes_cols), names_to="Condition", values_to="NES") %>%
    mutate(Condition = sub(paste0("^", nes_prefix), "", Condition))
  df_padj <- df %>%
    tidyr::pivot_longer(all_of(padj_cols), names_to="Condition", values_to="padj") %>%
    mutate(Condition = sub(paste0("^", padj_prefix), "", Condition))
  df_long <- left_join(df_nes, df_padj, by=c("PATHWAY","Condition"))
  
  # this ensures the facet order matches the color order.
  df_long$Condition <- factor(df_long$Condition, levels = conditions)
  
  # 4) SIGNIFICANCE LEVEL
  df_long <- df_long %>%
    mutate(sig_level = case_when(
      padj < padj_threshold ~ paste0("padj<", padj_threshold),
      # padj < 0.01  ~ "padj<0.01",
      # padj < 0.05  ~ "padj<0.05",
      TRUE         ~ "ns"
    ))
  
  # 5) HIERARCHICAL CLUSTERING OF ROWS BASED ON NES
  nes_mat <- df %>%
    select(all_of(nes_cols)) %>%
    as.matrix()
  rownames(nes_mat) <- df$PATHWAY
  
  if (nrow(nes_mat) < 2) {
    message("  Skipping dotmatrix: NES matrix has less than 2 rows (", nrow(nes_mat), ") for final clustering.")
    return(invisible(NULL))
  }
  
  dist_nes <- dist(nes_mat)
  hc_nes   <- hclust(dist_nes, method = "ward.D2")
  row_order <- hc_nes$labels[hc_nes$order]
  
  # 6) APPLY ROW ORDER
  df_long$PATHWAY <- factor(df_long$PATHWAY, levels = rev(row_order))
  
  # Build strip fill vector in order of provided conditions
  strip_fills <- unname(DIAGNOSIS_COLKEY[conditions])
  
  # 7) PLOT
  p <- ggplot(df_long, aes(x = Condition, y = PATHWAY,
                           color = NES,
                           size  = -log10(padj),
                           shape = sig_level)) +
    geom_point() +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, name = "NES") +
    scale_size_continuous(name = expression(-log[10](padj))) +
    scale_shape_manual(
      values = setNames(c(1,16), c("ns", paste0("padj<", padj_threshold))),
      breaks = c("ns", paste0("padj<", padj_threshold)),
      name   = "Significance"
    ) +
    labs(x = NULL, y = NULL, title = "GSEA Dotmatrix") +
    theme_bw() +
    theme(
      axis.text.y   = element_text(size = 10),
      axis.text.x   = element_blank(),
      axis.ticks.x  = element_blank(),
      plot.title    = element_text(hjust = 0.5)
    ) +
    ggh4x::facet_grid2(
      . ~ Condition,
      scales = "free_x",
      space  = "free_x",
      switch = "x",
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = strip_fills
        )
      )
    )
  
  print(p)
}


# plot a single-column gsea dotmatrix for a given condition using the
# top_n up‑regulated (highest NES) and top_n down‑regulated (lowest NES)
# pathways (excluding any that contain "TBA"). the function re‑uses
# `plot_gsea_dotmatrix()` but supplies a reduced pathway set and a single
# condition vector. jaccard clustering is disabled, yet hierarchical
# clustering on the NES values (within `plot_gsea_dotmatrix()`) is still
# performed to keep the rows readable.

plot_gsea_top_bottom_single <- function(gsea_result,
                                        condition,
                                        top_n          = 10,
                                        padj_prefix    = "padj.",
                                        nes_prefix     = "NES.",
                                        padj_threshold = 0.001) {
  nes_col  <- paste0(nes_prefix,  condition)
  padj_col <- paste0(padj_prefix, condition)
  
  # guard: skip if the required columns are missing
  if (!all(c(nes_col, padj_col) %in% colnames(gsea_result))) {
    warning(sprintf("columns %s / %s not found in gsea_result — skipping", nes_col, padj_col))
    return(invisible(NULL))
  }
  
  df <- gsea_result[, c("PATHWAY", nes_col, padj_col)]
  colnames(df) <- c("PATHWAY", "NES", "padj")
  
  # remove missing & TBA entries ------------------------------------------------
  df <- df[!is.na(df$NES) & !grepl("TBA", df$PATHWAY, ignore.case = TRUE), ]
  if (nrow(df) == 0) {
    warning(sprintf("no pathways remain for %s after filtering", condition))
    return(invisible(NULL))
  }
  
  # select top_n up & down pathways -------------------------------------------
  df_up   <- head(df[order(-df$NES), ], top_n)
  df_down <- head(df[order(df$NES),  ], top_n)
  sel_pw  <- unique(c(df_up$PATHWAY, df_down$PATHWAY))
  
  # rebuild a wide df compatible with plot_gsea_dotmatrix ----------------------
  df_wide <- df[df$PATHWAY %in% sel_pw, ]
  df_wide <- df_wide %>%
    dplyr::rename(!!paste0("NES.",  condition) := NES,
           !!paste0("padj.", condition) := padj)
  
  # call existing plotting routine --------------------------------------------
  p <- plot_gsea_dotmatrix(
    gsea_result     = df_wide,
    conditions      = condition,
    padj_threshold  = padj_threshold,
    selected_pathways = sel_pw,
    cluster         = FALSE,           # disable jaccard clustering
    cluster_height  = 1                # ignored because cluster = FALSE
  )
  
  return(p)
}

plot_gsea_combined_top_bottom <- function(gsea_result,
                                          conditions,
                                          top_n          = 10,
                                          padj_prefix    = "padj.",
                                          nes_prefix     = "NES.",
                                          padj_threshold = 0.001) {
  sel_pw <- character()
  for (cond in conditions) {
    nes_col <- paste0(nes_prefix, cond)
    if (!nes_col %in% colnames(gsea_result)) next
    tmp <- gsea_result[, c("PATHWAY", nes_col)]
    colnames(tmp) <- c("PATHWAY", "NES")
    tmp <- tmp[!is.na(tmp$NES) & !grepl("TBA", tmp$PATHWAY, ignore.case = TRUE), ]
    if (nrow(tmp) == 0) next
    sel_pw <- c(sel_pw,
                head(tmp[order(-tmp$NES), ], top_n)$PATHWAY,
                head(tmp[order(tmp$NES),  ], top_n)$PATHWAY)
  }
  sel_pw <- unique(sel_pw)
  if (length(sel_pw) == 0) return(invisible(NULL))
  
  p <- plot_gsea_dotmatrix(
    gsea_result       = gsea_result,
    conditions        = conditions,
    padj_threshold    = padj_threshold,
    selected_pathways = sel_pw,
    cluster           = FALSE
  )
  return(p)
}

plot_gene_count_metadata_cor_heatmap <- function(count_matrix, metadata, genes_of_interest, metadata_columns = NULL, method = "spearman") {
  # Extract GSVA scores for pathways of interest
  gene_counts <- count_matrix[genes_of_interest, , drop = FALSE]
  
  # Transpose to match sample rows
  gene_counts_t <- t(gene_counts)
  
  # Extract metadata and ensure sample IDs align
  metadata <- metadata[match(rownames(gene_counts_t), rownames(metadata)), , drop = FALSE]
  
  # Filter metadata to numeric columns only or user-specified
  if (is.null(metadata_columns)) {
    metadata_num <- metadata[, sapply(metadata, is.numeric), drop = FALSE]
  } else {
    metadata_num <- metadata[, metadata_columns, drop = FALSE]
  }
  
  # Calculate correlation
  cor_matrix <- matrix(NA, nrow = length(genes_of_interest), ncol = ncol(metadata_num))
  rownames(cor_matrix) <- genes_of_interest
  colnames(cor_matrix) <- colnames(metadata_num)
  
  for (i in seq_along(genes_of_interest)) {
    for (j in seq_len(ncol(metadata_num))) {
      cor_matrix[i, j] <- suppressWarnings(cor(gene_counts_t[, i], metadata_num[, j], method = method, use = "pairwise.complete.obs"))
    }
  }
  
  # neutral variant:
  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  # cheeky variant:
  #col_fun <- colorRamp2(c(min(cor_matrix), 0, max(cor_matrix)), c("#2166AC", "#F7F7F7", "#B2182B"))
  #col_fun <- colorRamp2(c(min(cor_matrix), 0, max(cor_matrix)), c("blue", "white", "red"))
  
  ht <- Heatmap(cor_matrix,
                name = "Spearman\ncorrelation",  # legend title
                col = col_fun,
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                row_names_gp = gpar(fontsize = 10),
                row_names_side = "left",
                row_dend_side = "right",
                column_names_gp = gpar(fontsize = 10))
  
  return(ht)
}

# prepare data ------------------------------------------------------------
# load dds
dds <- readRDS("output/dds_qc.rds")

# define V1 subset
V1_dds <- subset_dds(dds, condition = "VISIT == 'V1'")

V1_count_matrix <- as.data.frame(counts(V1_dds, normalized = TRUE))

V1_metadata <- as.data.frame(colData(V1_dds))

# set to TRUE for debugging (much faster)
V1_vsd <- as.data.frame(assay(vst(V1_dds, blind = TRUE)))

# load l2fc data
l2fc_data_imputed <- readRDS("output/l2fc_data_imputed_qc.rds")

# load gsva scores
pathway_data <- readRDS("output/pathway_data.rds")
# create V1 subset
V1_pathway_data <- subset_pathway_data(pathway_data, condition = "VISIT == 'V1'")

# correlating changes in gene expression with changes in antibody levels
l2fc_pathway_data <- readRDS("output/l2fc_pathway_data.rds")


# Count corr TRI_ELISA enrichment --------------------------------

tri_variables = c("TRI_ELISA_V3", "TRI_ELISA_V4", "TRI_ELISA_V5")

gsea_results_count_tri_corr_file <- "output/gsea_results_count_tri_corr.rds"
gsea_results_count_tri_corr_cohort_file <- "output/gsea_results_count_tri_corr_cohort.rds"

if (!((file.exists(gsea_results_count_tri_corr_file)) && (file.exists(gsea_results_count_tri_corr_cohort_file)))) {
  
  message("Computing gsea for count TRI corr")
  
  # make sure dds is loaded (at the start)
  
  # initialize containers for correlation-based GSEA
  gsea_results_count_tri_corr        <- list()
  gsea_results_count_tri_corr_cohort <- list()
  
  # loop over each TRI_ELISA variable
  for (tri_variable in tri_variables) {
    message(sprintf("\n--- TRI corr enrichment for: %s ---", tri_variable))
    gsea_results_count_tri_corr[[tri_variable]]        <- list()
    gsea_results_count_tri_corr_cohort[[tri_variable]] <- list()
    
    # loop over visits of interest
    for (visit in VISITS[1:5]) {
      message(sprintf(" Processing visit: %s", visit))
      
      # --- Cohort-wide correlation ---
      dds_cohort <- subset_dds(dds, paste0("VISIT == '", visit, "'"))
      tri_vec   <- colData(dds_cohort)[[tri_variable]]
      expr_mat  <- counts(dds_cohort, normalized = TRUE)
      
      # ensure dimensions align
      if (length(tri_vec) != ncol(expr_mat)) stop("Length of TRI vector does not match number of samples.")
      
      # compute Spearman or Pearson correlation and p-values per gene
      cors_cohort  <- apply(expr_mat, 1, function(g) cor(g, tri_vec, method = "spearman", use = "pairwise.complete.obs"))
      pvals_cohort <- apply(expr_mat, 1, function(g) cor.test(g, tri_vec, method = "spearman")$p.value)
      
      # combine effect size (correlation) and significance: signed -log10(p-value), beware of the sign, but should be pos. within pval range (0, 1)
      ranking_metric_cohort <- sign(cors_cohort) * -log10(pvals_cohort)
      ranking_metric_cohort[!is.finite(ranking_metric_cohort)] <- NA
      geneList_cohort <- sort(na.omit(ranking_metric_cohort), decreasing = TRUE)
      
      # run fgsea enrichment using the combined ranking metric
      fgsea_cohort <- fgsea(pathways = PATHWAYS, stats = geneList_cohort, nperm = 10000)
      gsea_results_count_tri_corr_cohort[[tri_variable]][[visit]] <- fgsea_cohort
      
      # --- Diagnosis-specific correlations ---
      gsea_results_count_tri_corr[[tri_variable]][[visit]] <- list()
      for (diag in DIAGNOSES) {
        message(sprintf("  - Diagnosis: %s", diag))
        dds_diag  <- subset_dds(
          dds,
          paste0("(VISIT == '", visit, "') & (DIAGNOSIS == '", diag, "')")
        )
        tri_vec_d  <- colData(dds_diag)[[tri_variable]]
        expr_mat_d <- counts(dds_diag, normalized = TRUE)
        if (length(tri_vec_d) < 3) {
          message(sprintf("    Skipping %s (n < 3 samples)", diag))
          next
        }
        cors_diag  <- apply(expr_mat_d, 1, function(g) cor(g, tri_vec_d, method = "spearman", use = "pairwise.complete.obs"))
        pvals_diag <- apply(expr_mat_d, 1, function(g) cor.test(g, tri_vec_d, method = "spearman")$p.value)
        
        # combine into ranking metric
        ranking_metric_diag <- sign(cors_diag) * -log10(pvals_diag)
        ranking_metric_diag[!is.finite(ranking_metric_diag)] <- NA
        geneList_diag <- sort(na.omit(ranking_metric_diag), decreasing = TRUE)
        
        # run fgsea
        fgsea_diag <- fgsea(pathways = PATHWAYS, stats = geneList_diag, nperm = 10000)
        gsea_results_count_tri_corr[[tri_variable]][[visit]][[diag]] <- fgsea_diag
      }
    }
  }
  
  saveRDS(gsea_results_count_tri_corr, gsea_results_count_tri_corr_file)
  saveRDS(gsea_results_count_tri_corr_cohort, gsea_results_count_tri_corr_cohort_file)
  
} else {
  message("Loaded existing gsea results for count TRI corr")
  gsea_results_count_tri_corr <- readRDS(gsea_results_count_tri_corr_file)
  gsea_results_count_tri_corr_cohort <- readRDS(gsea_results_count_tri_corr_cohort_file)
}

# combine and plot dotmatrices for correlation-based enrichment
for (tri_variable in tri_variables) {
  for (visit in VISITS[1:5]) {
    message(sprintf("Plotting corr dotmatrix: %s at %s", tri_variable, visit))
    combined_gsea <- list()
    
    # diagnosis‑specific results --------------------------------------------
    for (diag in DIAGNOSES) {
      res_diag <- gsea_results_count_tri_corr[[tri_variable]][[visit]][[diag]]
      if (is.null(res_diag)) next
      res_df <- as.data.frame(res_diag)[, c("pathway", "NES", "padj")]
      colnames(res_df) <- c("PATHWAY", paste0("NES.", diag), paste0("padj.", diag))
      combined_gsea[[diag]] <- res_df
      
      # single‑column plot ---------------------------------------------------
      p_single <- plot_gsea_top_bottom_single(res_df, diag, padj_threshold = 0.05)
      if (!is.null(p_single)) {
        ggsave(sprintf("figs/correlation/tri/%s_%s_%s_single_dotmatrix.png",
                       tri_variable, visit, diag),
               p_single + ggtitle(sprintf("%s – %s (visit %s)", tri_variable, diag, visit)),
               width = 7, height = 6, dpi = 300)
      }
    }
    
    # cohort result ---------------------------------------------------------
    res_cohort <- gsea_results_count_tri_corr_cohort[[tri_variable]][[visit]]
    if (!is.null(res_cohort)) {
      res_cohort <- as.data.frame(res_cohort)[, c("pathway", "NES", "padj")]
      colnames(res_cohort) <- c("PATHWAY", "NES.COHORT", "padj.COHORT")
      combined_gsea[["COHORT"]] <- res_cohort
      
      # cohort single plot ---------------------------------------------------
      p_co <- plot_gsea_top_bottom_single(res_cohort, "COHORT", padj_threshold = 0.05)
      if (!is.null(p_co)) {
        ggsave(sprintf("figs/correlation/tri/%s_%s_COHORT_single_dotmatrix.png",
                       tri_variable, visit),
               p_co + ggtitle(sprintf("%s – COHORT (visit %s)", tri_variable, visit)),
               width = 7, height = 6, dpi = 300)
      }
    }
    
    # merged df -------------------------------------------------------------
    merged_df <- Reduce(function(x,y) full_join(x,y,by="PATHWAY"), combined_gsea)
    conditions <- c(DIAGNOSES, "COHORT")[c(DIAGNOSES, "COHORT") %in%
                                           sub(".*\\.(.*)$", "\\1", grep("^NES\\.", colnames(merged_df), value = TRUE))]
    
    # full dotmatrix --------------------------------------------------------
    dm <- plot_gsea_dotmatrix(merged_df, conditions, padj_threshold = 0.05,
                              pathway2genes = PATHWAYS, cluster = TRUE, cluster_height = 0.99)
    if (!is.null(dm)) {
      ggsave(sprintf("figs/correlation/tri/%s_%s_count_corr_dotmatrix_combined.png",
                     tri_variable, visit),
             dm + ggtitle(sprintf("%s: TRI corr enrichment with %s", tri_variable, visit)),
             width = 10, height = 6, dpi = 300)
    }
    
    # combined ±10 per group -----------------------------------------------
    dm_comb2 <- plot_gsea_combined_top_bottom(merged_df, conditions, top_n = 10, padj_threshold = 0.05)
    if (!is.null(dm_comb2)) {
      ggsave(sprintf("figs/correlation/tri/%s_%s_count_corr_dotmatrix_combined_2.png",
                     tri_variable, visit),
             dm_comb2 + ggtitle(sprintf("%s: TRI corr enrichment %s (±10 each)", tri_variable, visit)),
             width = 10, height = 10, dpi = 300)
    }
  }
}

# l2fc count corr TRI_ELISA enrichment --------------------------------

tri_variables = c("TRI_ELISA_V3", "TRI_ELISA_V4", "TRI_ELISA_V5")

gsea_results_l2fc_count_tri_corr_file <- "output/gsea_results_l2fc_count_tri_corr.rds"
gsea_results_l2fc_count_tri_corr_cohort_file <- "output/gsea_results_l2fc_tri_corr_cohort.rds"

if (!((file.exists(gsea_results_l2fc_count_tri_corr_file)) && (file.exists(gsea_results_l2fc_count_tri_corr_cohort_file)))) {
  
  message("Computing gsea for l2fc count TRI corr")
  
  # make sure imputed log2FC data is loaded (at the start)
  
  # initialize containers for correlation-based GSEA
  gsea_results_l2fc_count_tri_corr        <- list()
  gsea_results_l2fc_count_tri_corr_cohort  <- list()
  
  # loop over each TRI_ELISA variable
  for (tri_variable in tri_variables) {
    message(sprintf("\n--- TRI corr enrichment for: %s ---", tri_variable))
    gsea_results_l2fc_count_tri_corr[[tri_variable]]        <- list()
    gsea_results_l2fc_count_tri_corr_cohort[[tri_variable]] <- list()
    
    # loop over visits of interest
    for (visit in VISITS[2:5]) {
      message(sprintf(" Processing visit: %s", visit))
      
      # --- Cohort-wide correlation ---
      df_cohort <- subset_l2fc_data(l2fc_data_imputed, paste0("VISIT == '", visit, "'"))
      tri_vec   <- df_cohort$metadata[[tri_variable]]
      expr_mat  <- df_cohort$count_matrix
      
      # ensure dimensions align
      if (length(tri_vec) != ncol(expr_mat)) stop("Length of TRI vector does not match number of samples.")
      
      # compute Spearman or Pearson correlation and p-values per gene
      cors_cohort  <- apply(expr_mat, 1, function(g) cor(g, tri_vec, method = "spearman", use = "pairwise.complete.obs"))
      pvals_cohort <- apply(expr_mat, 1, function(g) cor.test(g, tri_vec, method = "spearman")$p.value)
      
      # combine effect size (correlation) and significance: signed -log10(p-value), beware of sign, but should be pos. in pval range (0,1)
      ranking_metric_cohort <- sign(cors_cohort) * -log10(pvals_cohort)
      ranking_metric_cohort[!is.finite(ranking_metric_cohort)] <- NA
      geneList_cohort <- sort(na.omit(ranking_metric_cohort), decreasing = TRUE)
      
      # run fgsea enrichment using the combined ranking metric
      fgsea_cohort <- fgsea(pathways = PATHWAYS, stats = geneList_cohort, nperm = 10000)
      gsea_results_l2fc_count_tri_corr_cohort[[tri_variable]][[visit]] <- fgsea_cohort
      
      # --- Diagnosis-specific correlations ---
      gsea_results_l2fc_count_tri_corr[[tri_variable]][[visit]] <- list()
      for (diag in DIAGNOSES) {
        message(sprintf("  - Diagnosis: %s", diag))
        df_diag  <- subset_l2fc_data(
          l2fc_data_imputed,
          paste0("(VISIT == '", visit, "') & (DIAGNOSIS == '", diag, "')")
        )
        tri_vec_d  <- df_diag$metadata[[tri_variable]]
        expr_mat_d <- df_diag$count_matrix
        if (length(tri_vec_d) < 3) {
          message(sprintf("    Skipping %s (n < 3 samples)", diag))
          next
        }
        cors_diag  <- apply(expr_mat_d, 1, function(g) cor(g, tri_vec_d, method = "spearman", use = "pairwise.complete.obs"))
        pvals_diag <- apply(expr_mat_d, 1, function(g) cor.test(g, tri_vec_d, method = "spearman")$p.value)
        
        # combine into ranking metric
        ranking_metric_diag <- sign(cors_diag) * -log10(pvals_diag)
        ranking_metric_diag[!is.finite(ranking_metric_diag)] <- NA
        geneList_diag <- sort(na.omit(ranking_metric_diag), decreasing = TRUE)
        
        # run fgsea
        fgsea_diag <- fgsea(pathways = PATHWAYS, stats = geneList_diag, nperm = 10000)
        gsea_results_l2fc_count_tri_corr[[tri_variable]][[visit]][[diag]] <- fgsea_diag
      }
    }
  }
  
  saveRDS(gsea_results_l2fc_count_tri_corr, gsea_results_l2fc_count_tri_corr_file)
  saveRDS(gsea_results_l2fc_count_tri_corr_cohort, gsea_results_l2fc_count_tri_corr_cohort_file)
  
} else {
  message("Loaded existing gsea results for l2fc count TRI corr")
  gsea_results_l2fc_count_tri_corr <- readRDS(gsea_results_l2fc_count_tri_corr_file)
  gsea_results_l2fc_count_tri_corr_cohort <- readRDS(gsea_results_l2fc_count_tri_corr_cohort_file)
}

# combine and plot dotmatrices for correlation-based enrichment
for (tri_variable in tri_variables) {
  for (visit in VISITS[2:5]) {
    message(sprintf("Plotting l2fc corr dotmatrix: %s %svsV1", tri_variable, visit))
    combined_gsea <- list()
    
    # diagnosis‑specific results --------------------------------------------
    for (diag in DIAGNOSES) {
      res_diag <- gsea_results_l2fc_count_tri_corr[[tri_variable]][[visit]][[diag]]
      if (is.null(res_diag)) next
      res_df <- as.data.frame(res_diag)[, c("pathway", "NES", "padj")]
      colnames(res_df) <- c("PATHWAY", paste0("NES.", diag), paste0("padj.", diag))
      combined_gsea[[diag]] <- res_df
      
      p_single <- plot_gsea_top_bottom_single(res_df, diag, padj_threshold = 0.05)
      if (!is.null(p_single)) {
        ggsave(sprintf("figs/correlation/tri/%s_%s_%s_l2fc_single_dotmatrix.png",
                       tri_variable, visit, diag),
               p_single + ggtitle(sprintf("%s – %s (%svsV1)", tri_variable, diag, visit)),
               width = 7, height = 6, dpi = 300)
      }
    }
    
    # cohort result ---------------------------------------------------------
    res_cohort <- gsea_results_l2fc_count_tri_corr_cohort[[tri_variable]][[visit]]
    if (!is.null(res_cohort)) {
      res_cohort <- as.data.frame(res_cohort)[, c("pathway", "NES", "padj")]
      colnames(res_cohort) <- c("PATHWAY", "NES.COHORT", "padj.COHORT")
      combined_gsea[["COHORT"]] <- res_cohort
      
      p_co <- plot_gsea_top_bottom_single(res_cohort, "COHORT", padj_threshold = 0.05)
      if (!is.null(p_co)) {
        ggsave(sprintf("figs/correlation/tri/%s_%s_COHORT_l2fc_single_dotmatrix.png",
                       tri_variable, visit),
               p_co + ggtitle(sprintf("%s – COHORT (%svsV1)", tri_variable, visit)),
               width = 7, height = 6, dpi = 300)
      }
    }
    
    merged_df <- Reduce(function(x,y) full_join(x,y,by="PATHWAY"), combined_gsea)
    conditions <- c(DIAGNOSES, "COHORT")[c(DIAGNOSES, "COHORT") %in%
                                           sub(".*\\.(.*)$", "\\1", grep("^NES\\.", colnames(merged_df), value = TRUE))]
    
    dm <- plot_gsea_dotmatrix(merged_df, conditions, padj_threshold = 0.05,
                              pathway2genes = PATHWAYS, cluster = TRUE, cluster_height = 0.98)
    if (!is.null(dm)) {
      ggsave(sprintf("figs/correlation/tri/%s_%s_l2fc_count_corr_dotmatrix.png",
                     tri_variable, visit),
             dm + ggtitle(sprintf("%s: TRI corr enrichment %svsV1 l2fc", tri_variable, visit)),
             width = 10, height = 6, dpi = 300)
    }
    
    dm_comb2 <- plot_gsea_combined_top_bottom(merged_df, conditions, top_n = 10, padj_threshold = 0.05)
    if (!is.null(dm_comb2)) {
      ggsave(sprintf("figs/correlation/tri/%s_%s_l2fc_count_corr_dotmatrix_combined_2.png",
                     tri_variable, visit),
             dm_comb2 + ggtitle(sprintf("%s: TRI corr enrichment %svsV1 l2fc (±10 each)", tri_variable, visit)),
             width = 10, height = 10, dpi = 300)
    }
  }
}

# Count corr IgG enrichment --------------------------------

igg_variables = c("IgG_V1", "IgG_V3", "IgG_V4")

gsea_results_count_igg_corr_file <- "output/gsea_results_count_igg_corr.rds"
gsea_results_count_igg_corr_cohort_file <- "output/gsea_results_count_igg_corr_cohort.rds"

if (!((file.exists(gsea_results_count_igg_corr_file)) && (file.exists(gsea_results_count_igg_corr_cohort_file)))) {
  
  message("Computing gsea for count IgG corr")
  
  # make sure dds is loaded (at the start)
  
  # initialize containers for correlation-based GSEA
  gsea_results_count_igg_corr         <- list()
  gsea_results_count_igg_corr_cohort <- list()
  
  # loop over each IgG variable
  for (igg_variable in igg_variables) {
    message(sprintf("\n--- IgG corr enrichment for: %s ---", igg_variable))
    gsea_results_count_igg_corr[[igg_variable]]         <- list()
    gsea_results_count_igg_corr_cohort[[igg_variable]] <- list()
    
    # loop over visits of interest
    for (visit in VISITS[1:4]) {
      message(sprintf(" Processing visit: %s", visit))
      
      # --- Cohort-wide correlation ---
      dds_cohort <- subset_dds(dds, paste0("VISIT == '", visit, "'"))
      igg_vec   <- colData(dds_cohort)[[igg_variable]]
      expr_mat  <- counts(dds_cohort, normalized = TRUE)
      
      # ensure dimensions align
      if (length(igg_vec) != ncol(expr_mat)) stop("Length of IgG vector does not match number of samples.")
      
      # compute Spearman or Pearson correlation and p-values per gene
      cors_cohort  <- apply(expr_mat, 1, function(g) cor(g, igg_vec, method = "spearman", use = "pairwise.complete.obs"))
      pvals_cohort <- apply(expr_mat, 1, function(g) cor.test(g, igg_vec, method = "spearman")$p.value)
      
      # combine effect size (correlation) and significance: signed -log10(p-value), beware of the sign, but should be pos. within pval range (0, 1)
      ranking_metric_cohort <- sign(cors_cohort) * -log10(pvals_cohort)
      ranking_metric_cohort[!is.finite(ranking_metric_cohort)] <- NA
      geneList_cohort <- sort(na.omit(ranking_metric_cohort), decreasing = TRUE)
      
      # run fgsea enrichment using the combined ranking metric
      fgsea_cohort <- fgsea(pathways = PATHWAYS, stats = geneList_cohort, nperm = 10000)
      gsea_results_count_igg_corr_cohort[[igg_variable]][[visit]] <- fgsea_cohort
      
      # --- Diagnosis-specific correlations ---
      gsea_results_count_igg_corr[[igg_variable]][[visit]] <- list()
      for (diag in DIAGNOSES) {
        message(sprintf("  - Diagnosis: %s", diag))
        dds_diag  <- subset_dds(
          dds,
          paste0("(VISIT == '", visit, "') & (DIAGNOSIS == '", diag, "')")
        )
        igg_vec_d  <- colData(dds_diag)[[igg_variable]]
        expr_mat_d <- counts(dds_diag, normalized = TRUE)
        if (length(igg_vec_d) < 3) {
          message(sprintf("    Skipping %s (n < 3 samples)", diag))
          next
        }
        cors_diag  <- apply(expr_mat_d, 1, function(g) cor(g, igg_vec_d, method = "spearman", use = "pairwise.complete.obs"))
        pvals_diag <- apply(expr_mat_d, 1, function(g) cor.test(g, igg_vec_d, method = "spearman")$p.value)
        
        # combine into ranking metric
        ranking_metric_diag <- sign(cors_diag) * -log10(pvals_diag)
        ranking_metric_diag[!is.finite(ranking_metric_diag)] <- NA
        geneList_diag <- sort(na.omit(ranking_metric_diag), decreasing = TRUE)
        
        # run fgsea
        fgsea_diag <- fgsea(pathways = PATHWAYS, stats = geneList_diag, nperm = 10000)
        gsea_results_count_igg_corr[[igg_variable]][[visit]][[diag]] <- fgsea_diag
      }
    }
  }
  
  saveRDS(gsea_results_count_igg_corr, gsea_results_count_igg_corr_file)
  saveRDS(gsea_results_count_igg_corr_cohort, gsea_results_count_igg_corr_cohort_file)
  
} else {
  message("Loaded existing gsea results for count IgG corr")
  gsea_results_count_igg_corr <- readRDS(gsea_results_count_igg_corr_file)
  gsea_results_count_igg_corr_cohort <- readRDS(gsea_results_count_igg_corr_cohort_file)
}

# combine and plot dotmatrices for correlation-based enrichment
for (igg_variable in igg_variables) {
  for (visit in VISITS[1:4]) {
    message(sprintf("Plotting corr dotmatrix: %s at %s", igg_variable, visit))
    combined_gsea <- list()
    
    # diagnosis-specific results --------------------------------------------
    for (diag in DIAGNOSES) {
      res_diag <- gsea_results_count_igg_corr[[igg_variable]][[visit]][[diag]]
      if (is.null(res_diag)) next
      res_df <- as.data.frame(res_diag)[, c("pathway", "NES", "padj")]
      colnames(res_df) <- c("PATHWAY", paste0("NES.", diag), paste0("padj.", diag))
      combined_gsea[[diag]] <- res_df
      
      # single-column plot ---------------------------------------------------
      p_single <- plot_gsea_top_bottom_single(res_df, diag, padj_threshold = 0.05)
      if (!is.null(p_single)) {
        ggsave(sprintf("figs/correlation/igg/%s_%s_%s_single_dotmatrix.png",
                       igg_variable, visit, diag),
               p_single + ggtitle(sprintf("%s – %s (visit %s)", igg_variable, diag, visit)),
               width = 7, height = 6, dpi = 300)
      }
    }
    
    # cohort result ---------------------------------------------------------
    res_cohort <- gsea_results_count_igg_corr_cohort[[igg_variable]][[visit]]
    if (!is.null(res_cohort)) {
      res_cohort <- as.data.frame(res_cohort)[, c("pathway", "NES", "padj")]
      colnames(res_cohort) <- c("PATHWAY", "NES.COHORT", "padj.COHORT")
      combined_gsea[["COHORT"]] <- res_cohort
      
      # cohort single plot ---------------------------------------------------
      p_co <- plot_gsea_top_bottom_single(res_cohort, "COHORT", padj_threshold = 0.05)
      if (!is.null(p_co)) {
        ggsave(sprintf("figs/correlation/igg/%s_%s_COHORT_single_dotmatrix.png",
                       igg_variable, visit),
               p_co + ggtitle(sprintf("%s – COHORT (visit %s)", igg_variable, visit)),
               width = 7, height = 6, dpi = 300)
      }
    }
    
    # merged df -------------------------------------------------------------
    merged_df <- Reduce(function(x,y) full_join(x,y,by="PATHWAY"), combined_gsea)
    conditions <- c(DIAGNOSES, "COHORT")[c(DIAGNOSES, "COHORT") %in%
                                           sub(".*\\.(.*)$", "\\1", grep("^NES\\.", colnames(merged_df), value = TRUE))]
    
    # full dotmatrix --------------------------------------------------------
    dm <- plot_gsea_dotmatrix(merged_df, conditions, padj_threshold = 0.05,
                              pathway2genes = PATHWAYS, cluster = TRUE, cluster_height = 0.99)
    if (!is.null(dm)) {
      ggsave(sprintf("figs/correlation/igg/%s_%s_count_corr_dotmatrix_combined.png",
                     igg_variable, visit),
             dm + ggtitle(sprintf("%s: IgG corr enrichment with %s", igg_variable, visit)),
             width = 10, height = 6, dpi = 300)
    }
    
    # combined ±10 per group -----------------------------------------------
    dm_comb2 <- plot_gsea_combined_top_bottom(merged_df, conditions, top_n = 10, padj_threshold = 0.05)
    if (!is.null(dm_comb2)) {
      ggsave(sprintf("figs/correlation/igg/%s_%s_count_corr_dotmatrix_combined_2.png",
                     igg_variable, visit),
             dm_comb2 + ggtitle(sprintf("%s: IgG corr enrichment %s (±10 each)", igg_variable, visit)),
             width = 10, height = 10, dpi = 300)
    }
  }
}
# l2fc count corr IgG enrichment --------------------------------

igg_variables = c("IgG_V3", "IgG_V4")

gsea_results_l2fc_count_igg_corr_file <- "output/gsea_results_l2fc_count_igg_corr.rds"
gsea_results_l2fc_count_igg_corr_cohort_file <- "output/gsea_results_l2fc_igg_corr_cohort.rds"

if (!((file.exists(gsea_results_l2fc_count_igg_corr_file)) && (file.exists(gsea_results_l2fc_count_igg_corr_cohort_file)))) {
  
  message("Computing gsea for l2fc count IgG corr")
  
  # make sure imputed log2FC data is loaded (at the start)
  
  # initialize containers for correlation-based GSEA
  gsea_results_l2fc_count_igg_corr         <- list()
  gsea_results_l2fc_count_igg_corr_cohort  <- list()
  
  # loop over each IgG variable
  for (igg_variable in igg_variables) {
    message(sprintf("\n--- IgG corr enrichment for: %s ---", igg_variable))
    gsea_results_l2fc_count_igg_corr[[igg_variable]]         <- list()
    gsea_results_l2fc_count_igg_corr_cohort[[igg_variable]] <- list()
    
    # loop over visits of interest
    for (visit in VISITS[2:4]) {
      message(sprintf(" Processing visit: %s", visit))
      
      # --- Cohort-wide correlation ---
      df_cohort <- subset_l2fc_data(l2fc_data_imputed, paste0("VISIT == '", visit, "'"))
      igg_vec   <- df_cohort$metadata[[igg_variable]]
      expr_mat  <- df_cohort$count_matrix
      
      # ensure dimensions align
      if (length(igg_vec) != ncol(expr_mat)) stop("Length of IgG vector does not match number of samples.")
      
      # compute Spearman or Pearson correlation and p-values per gene
      cors_cohort  <- apply(expr_mat, 1, function(g) cor(g, igg_vec, method = "spearman", use = "pairwise.complete.obs"))
      pvals_cohort <- apply(expr_mat, 1, function(g) cor.test(g, igg_vec, method = "spearman")$p.value)
      
      # combine effect size (correlation) and significance: signed -log10(p-value), beware of sign, but should be pos. in pval range (0,1)
      ranking_metric_cohort <- sign(cors_cohort) * -log10(pvals_cohort)
      ranking_metric_cohort[!is.finite(ranking_metric_cohort)] <- NA
      geneList_cohort <- sort(na.omit(ranking_metric_cohort), decreasing = TRUE)
      
      # run fgsea enrichment using the combined ranking metric
      fgsea_cohort <- fgsea(pathways = PATHWAYS, stats = geneList_cohort, nperm = 10000)
      gsea_results_l2fc_count_igg_corr_cohort[[igg_variable]][[visit]] <- fgsea_cohort
      
      # --- Diagnosis-specific correlations ---
      gsea_results_l2fc_count_igg_corr[[igg_variable]][[visit]] <- list()
      for (diag in DIAGNOSES) {
        message(sprintf("  - Diagnosis: %s", diag))
        df_diag  <- subset_l2fc_data(
          l2fc_data_imputed,
          paste0("(VISIT == '", visit, "') & (DIAGNOSIS == '", diag, "')")
        )
        igg_vec_d  <- df_diag$metadata[[igg_variable]]
        expr_mat_d <- df_diag$count_matrix
        if (length(igg_vec_d) < 3) {
          message(sprintf("    Skipping %s (n < 3 samples)", diag))
          next
        }
        cors_diag  <- apply(expr_mat_d, 1, function(g) cor(g, igg_vec_d, method = "spearman", use = "pairwise.complete.obs"))
        pvals_diag <- apply(expr_mat_d, 1, function(g) cor.test(g, igg_vec_d, method = "spearman")$p.value)
        
        # combine into ranking metric
        ranking_metric_diag <- sign(cors_diag) * -log10(pvals_diag)
        ranking_metric_diag[!is.finite(ranking_metric_diag)] <- NA
        geneList_diag <- sort(na.omit(ranking_metric_diag), decreasing = TRUE)
        
        # run fgsea
        fgsea_diag <- fgsea(pathways = PATHWAYS, stats = geneList_diag, nperm = 10000)
        gsea_results_l2fc_count_igg_corr[[igg_variable]][[visit]][[diag]] <- fgsea_diag
      }
    }
  }
  
  saveRDS(gsea_results_l2fc_count_igg_corr, gsea_results_l2fc_count_igg_corr_file)
  saveRDS(gsea_results_l2fc_count_igg_corr_cohort, gsea_results_l2fc_count_igg_corr_cohort_file)
  
} else {
  message("Loaded existing gsea results for l2fc count IgG corr")
  gsea_results_l2fc_count_igg_corr <- readRDS(gsea_results_l2fc_count_igg_corr_file)
  gsea_results_l2fc_count_igg_corr_cohort <- readRDS(gsea_results_l2fc_count_igg_corr_cohort_file)
}

# combine and plot dotmatrices for correlation-based enrichment
for (igg_variable in igg_variables) {
  for (visit in VISITS[2:4]) {
    message(sprintf("Plotting l2fc corr dotmatrix: %s %svsV1", igg_variable, visit))
    combined_gsea <- list()
    
    # diagnosis-specific results --------------------------------------------
    for (diag in DIAGNOSES) {
      res_diag <- gsea_results_l2fc_count_igg_corr[[igg_variable]][[visit]][[diag]]
      if (is.null(res_diag)) next
      res_df <- as.data.frame(res_diag)[, c("pathway", "NES", "padj")]
      colnames(res_df) <- c("PATHWAY", paste0("NES.", diag), paste0("padj.", diag))
      combined_gsea[[diag]] <- res_df
      
      p_single <- plot_gsea_top_bottom_single(res_df, diag, padj_threshold = 0.05)
      if (!is.null(p_single)) {
        ggsave(sprintf("figs/correlation/igg/%s_%s_%s_l2fc_single_dotmatrix.png",
                       igg_variable, visit, diag),
               p_single + ggtitle(sprintf("%s – %s (%svsV1)", igg_variable, diag, visit)),
               width = 7, height = 6, dpi = 300)
      }
    }
    
    # cohort result ---------------------------------------------------------
    res_cohort <- gsea_results_l2fc_count_igg_corr_cohort[[igg_variable]][[visit]]
    if (!is.null(res_cohort)) {
      res_cohort <- as.data.frame(res_cohort)[, c("pathway", "NES", "padj")]
      colnames(res_cohort) <- c("PATHWAY", "NES.COHORT", "padj.COHORT")
      combined_gsea[["COHORT"]] <- res_cohort
      
      p_co <- plot_gsea_top_bottom_single(res_cohort, "COHORT", padj_threshold = 0.05)
      if (!is.null(p_co)) {
        ggsave(sprintf("figs/correlation/igg/%s_%s_COHORT_l2fc_single_dotmatrix.png",
                       igg_variable, visit),
               p_co + ggtitle(sprintf("%s – COHORT (%svsV1)", igg_variable, visit)),
               width = 7, height = 6, dpi = 300)
      }
    }
    
    merged_df <- Reduce(function(x,y) full_join(x,y,by="PATHWAY"), combined_gsea)
    conditions <- c(DIAGNOSES, "COHORT")[c(DIAGNOSES, "COHORT") %in%
                                           sub(".*\\.(.*)$", "\\1", grep("^NES\\.", colnames(merged_df), value = TRUE))]
    
    dm <- plot_gsea_dotmatrix(merged_df, conditions, padj_threshold = 0.05,
                              pathway2genes = PATHWAYS, cluster = TRUE, cluster_height = 0.98)
    if (!is.null(dm)) {
      ggsave(sprintf("figs/correlation/igg/%s_%s_l2fc_count_corr_dotmatrix.png",
                     igg_variable, visit),
             dm + ggtitle(sprintf("%s: IgG corr enrichment %svsV1 l2fc", igg_variable, visit)),
             width = 10, height = 6, dpi = 300)
    }
    
    dm_comb2 <- plot_gsea_combined_top_bottom(merged_df, conditions, top_n = 10, padj_threshold = 0.05)
    if (!is.null(dm_comb2)) {
      ggsave(sprintf("figs/correlation/igg/%s_%s_l2fc_count_corr_dotmatrix_combined_2.png",
                     igg_variable, visit),
             dm_comb2 + ggtitle(sprintf("%s: IgG corr enrichment %svsV1 l2fc (±10 each)", igg_variable, visit)),
             width = 10, height = 10, dpi = 300)
    }
  }
}
# Spearman correlation of V1 gene counts scores with TRI_ELISA within diagnosis -------------------------------------------


pathways_of_interest <- c("enriched in B cells (I) (M47.0)",
                          "T cell activation (I) (M7.1)",
                          "enriched in NK cells (I) (M7.2)",
                          "enriched for ubiquitination (M138)",
                          "proteasome (M226)",
                          "cell cycle and transcription (M4.0)",
                          "transcription regulation in cell development (M49)",
                          "immune activation - generic cluster (M37.0)",
                          "respiratory electron transport chain (mitochondrion) (M219)",
                          "platelet activation (II) (M32.1)",
                          "enriched in monocytes (II) (M11.0)",
                          "enriched in antigen presentation (I) (M71)")


# Create gene heatmaps
for (pathway_name in pathways_of_interest) {
  # Extract the M ID from the pathway name
  m_id_match <- regmatches(pathway_name, regexec("\\((M[0-9]+\\.?[0-9]*)\\)", pathway_name))
  pathway_m_id <- if (length(m_id_match) > 0 && length(m_id_match[[1]]) > 1) m_id_match[[1]][2] else gsub("[^A-Za-z0-9_.-]", "", pathway_name) # Fallback if M ID not found
  
  genes_of_interest <- PATHWAYS[[pathway_name]]
  genes_of_interest <- intersect(genes_of_interest, rownames(V1_count_matrix))
  
  if (length(genes_of_interest) == 0) {
    warning("No genes found for pathway: ", pathway_name, ". Skipping heatmap generation.")
    next
  }
  
  # 1. prepare data
  V1_count_matrix_t <- as.data.frame(t(V1_count_matrix[genes_of_interest, ]))
  V1_count_matrix_t$CEGAT_ID <- rownames(V1_count_matrix_t)
  
  data_merged <- left_join(
    V1_count_matrix_t,
    V1_metadata %>% select(CEGAT_ID, TRI_ELISA_V4, DIAGNOSIS),
    by = "CEGAT_ID"
  ) %>% na.omit()
  
  diagnosis_groups <- unique(data_merged$DIAGNOSIS)
  
  # 2. perform per-group Spearman correlations
  correlation_results_list <- list()
  p_value_results_list <- list()
  
  for (diag in diagnosis_groups) {
    diag_data <- filter(data_merged, DIAGNOSIS == diag)
    if (nrow(diag_data) < 2) {
      warning("Skipping group '", diag, "' for pathway '", pathway_name, "': fewer than 2 samples.")
      correlation_results_list[[diag]] <- rep(NA, length(genes_of_interest))
      p_value_results_list[[diag]] <- rep(NA, length(genes_of_interest))
      next
    }
    cors <- pvals <- numeric(length(genes_of_interest))
    names(cors) <- names(pvals) <- genes_of_interest
    
    for (gene in genes_of_interest) {
      valid_idx <- !is.na(diag_data[[gene]]) & !is.na(diag_data$TRI_ELISA_V4)
      if (sum(valid_idx) >= 2) {
        res <- rcorr(diag_data[[gene]][valid_idx], diag_data$TRI_ELISA_V4[valid_idx], type = "spearman")
        cors[gene] <- res$r[1,2]
        pvals[gene] <- res$P[1,2]
      } else {
        cors[gene] <- NA
        pvals[gene] <- NA
      }
    }
    correlation_results_list[[diag]] <- cors
    p_value_results_list[[diag]] <- pvals
  }
  
  # 3. compute COHORT (all-samples) correlation
  cohort_cors <- pvals_all <- numeric(length(genes_of_interest))
  names(cohort_cors) <- names(pvals_all) <- genes_of_interest
  
  for (gene in genes_of_interest) {
    valid_idx <- !is.na(data_merged[[gene]]) & !is.na(data_merged$TRI_ELISA_V4)
    if (sum(valid_idx) >= 2) {
      res_all <- rcorr(data_merged[[gene]][valid_idx], data_merged$TRI_ELISA_V4[valid_idx], type = "spearman")
      cohort_cors[gene] <- res_all$r[1,2]
      pvals_all[gene] <- res_all$P[1,2]
    } else {
      cohort_cors[gene] <- NA
      pvals_all[gene] <- NA
    }
  }
  
  # add COHORT to lists
  correlation_results_list[["COHORT"]] <- cohort_cors
  p_value_results_list[["COHORT"]] <- pvals_all
  
  # 4. build matrices for heatmap
  correlation_matrix <- do.call(cbind, correlation_results_list)
  p_value_matrix <- do.call(cbind, p_value_results_list)
  
  rownames(correlation_matrix) <- genes_of_interest
  colnames(correlation_matrix) <- names(correlation_results_list)
  rownames(p_value_matrix) <- genes_of_interest
  colnames(p_value_matrix) <- names(p_value_results_list)
  
  # # - 4.1. filter to genes with >=1 significant correlation (p < 0.05) -
  # sig_genes_p_value <- rownames(p_value_matrix)[
  #   apply(p_value_matrix, 1, function(pv) any(pv < 0.05, na.rm = TRUE))
  # ]
  # 
  # # - 4.2. filter to genes with >=1 absolute correlation > 0.2 -
  # sig_genes_correlation_strength <- rownames(correlation_matrix)[
  #   apply(correlation_matrix, 1, function(cor) any(abs(cor) > 0.2, na.rm = TRUE))
  # ]
  # 
  # # Append/Combine: keep only genes that are in *both* lists (satisfy both conditions)
  # sig_genes <- intersect(sig_genes_p_value, sig_genes_correlation_strength)
  # 
  # 
  # # warn when nothing is significant:
  # if (length(sig_genes) == 0) {
  #   warning("No genes pass the significance filter for pathway: ", pathway_name, "; plotting all genes.")
  #   sig_genes <- rownames(p_value_matrix)
  # }
  # 
  # # subset both matrices
  # correlation_matrix <- correlation_matrix[sig_genes, , drop = FALSE]
  # p_value_matrix <- p_value_matrix[sig_genes, , drop = FALSE]
  # 
  # # check if there are still genes to plot after filtering
  # if (nrow(correlation_matrix) == 0) {
  #   warning("No genes to plot for pathway: ", pathway_name, " after filtering. Skipping heatmap generation.")
  #   next
  # }
  
  # 5. plot heatmap with significance stars
  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  ht <- Heatmap(correlation_matrix,
                name = "Spearman ρ",
                col = col_fun,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                column_title = paste0("Spearman Correlation of V1 Gene Counts with TRI_ELISA_V4\nPathway: ", pathway_name),
                column_title_gp = gpar(fontsize = 16, fontface = "bold"),
                row_names_gp = gpar(fontsize = 5),
                column_names_gp = gpar(fontsize = 10) #,
                # cell_fun = function(j, i, x, y, width, height, fill) {
                #   p_val <- p_value_matrix[i, j]
                #   if (!is.na(p_val) && p_val < 0.05) {
                #     stars <- if (p_val < 0.001) "***"
                #     else if (p_val < 0.01) "**"
                #     else "*"
                #     grid.text(stars, x, y, gp = gpar(fontsize = 10, col = "black"))
                #   }
                # }
                )
  
  # define filename using the M ID
  filename <- paste0("figs/correlation/pathway_genes_spearman_corr/", pathway_m_id, "_correlation_heatmap.png")
  
  png(filename, width = 1400, height = 1400, res = 300)
  draw(ht, padding = unit(c(10, 40, 10, 10), "mm"))
  dev.off()
}


# Spearman correlation of l2fc gene counts scores with TRI_ELISA within diagnosis -------------------------------------------

pathways_of_interest <- c(# V2 related
                          #"cell cycle and transcription (M4.0)" # also for V3, V4
                          #"enriched in T cells (I) (M7.0)", # also V4
                          #"enriched in B cells (I) (M47.0)", # also for V3, V4
                          #"enriched in activated dendritic cells (II) (M165)" # also V4
                          # V3 related
                          #"phosphatidylinositol signaling system (M101)",
                          #"golgi membrane (II) (M237)",
                          #"antiviral IFN signature (M75)",
                          #"immune activation - generic cluster (M37.0)",
                          #"Plasma cell surface signature (S3)"
                          # V4 related
                          # "heme biosynthesis (I) (M171)",
                          #"RA, WNT, CSF receptors network (monocyte) (M23)",
                          #"enriched in activated dendritic cells/monocytes (M64)",
                          #"inositol phosphate metabolism (M129)",
                          #"golgi membrane (II) (M237)"
                          )

for (visit in VISITS[2:4]) {
  
  l2fc_data_imputed_subset <- subset_l2fc_data(l2fc_data_imputed, condition = paste0("VISIT ==", "'", visit, "'"))
  
  for (pathway_name in pathways_of_interest) {
    # Extract the ID from the pathway name
    id_match <- regmatches(pathway_name, regexec("\\(([MS][0-9]+\\.?[0-9]*)\\)", pathway_name))
    pathway_id <- if (length(id_match) > 0 && length(id_match[[1]]) > 1) id_match[[1]][2] else gsub("[^A-Za-z0-9_.-]", "", pathway_name) # Fallback if ID not found
    
    genes_of_interest <- PATHWAYS[[pathway_name]]
    genes_of_interest <- intersect(genes_of_interest, rownames(l2fc_data_imputed_subset$count_matrix))
    
    if (length(genes_of_interest) == 0) {
      warning("No genes found for pathway: ", pathway_name, ". Skipping heatmap generation.")
      next
    }
    
    # 1. prepare data
    l2fc_data_imputed_subset$count_matrix_t <- as.data.frame(t(l2fc_data_imputed_subset$count_matrix[genes_of_interest, ]))
    l2fc_data_imputed_subset$count_matrix_t$CEGAT_ID <- rownames(l2fc_data_imputed_subset$count_matrix_t)
    
    data_merged <- left_join(
      l2fc_data_imputed_subset$count_matrix_t,
      l2fc_data_imputed_subset$metadata %>% select(CEGAT_ID, TRI_ELISA_V4, DIAGNOSIS),
      by = "CEGAT_ID"
    ) %>% na.omit()
    
    diagnosis_groups <- unique(data_merged$DIAGNOSIS)
    
    # 2. perform per-group Spearman correlations
    correlation_results_list <- list()
    p_value_results_list <- list()
    
    for (diag in diagnosis_groups) {
      diag_data <- filter(data_merged, DIAGNOSIS == diag)
      if (nrow(diag_data) < 2) {
        warning("Skipping group '", diag, "' for pathway '", pathway_name, "': fewer than 2 samples.")
        correlation_results_list[[diag]] <- rep(NA, length(genes_of_interest))
        p_value_results_list[[diag]] <- rep(NA, length(genes_of_interest))
        next
      }
      cors <- pvals <- numeric(length(genes_of_interest))
      names(cors) <- names(pvals) <- genes_of_interest
      
      for (gene in genes_of_interest) {
        valid_idx <- !is.na(diag_data[[gene]]) & !is.na(diag_data$TRI_ELISA_V4)
        if (sum(valid_idx) >= 2) {
          res <- rcorr(diag_data[[gene]][valid_idx], diag_data$TRI_ELISA_V4[valid_idx], type = "spearman")
          cors[gene] <- res$r[1,2]
          pvals[gene] <- res$P[1,2]
        } else {
          cors[gene] <- NA
          pvals[gene] <- NA
        }
      }
      correlation_results_list[[diag]] <- cors
      p_value_results_list[[diag]] <- pvals
    }
    
    # 3. compute COHORT (all-samples) correlation
    cohort_cors <- pvals_all <- numeric(length(genes_of_interest))
    names(cohort_cors) <- names(pvals_all) <- genes_of_interest
    
    for (gene in genes_of_interest) {
      valid_idx <- !is.na(data_merged[[gene]]) & !is.na(data_merged$TRI_ELISA_V4)
      if (sum(valid_idx) >= 2) {
        res_all <- rcorr(data_merged[[gene]][valid_idx], data_merged$TRI_ELISA_V4[valid_idx], type = "spearman")
        cohort_cors[gene] <- res_all$r[1,2]
        pvals_all[gene] <- res_all$P[1,2]
      } else {
        cohort_cors[gene] <- NA
        pvals_all[gene] <- NA
      }
    }
    
    # add COHORT to lists
    correlation_results_list[["COHORT"]] <- cohort_cors
    p_value_results_list[["COHORT"]] <- pvals_all
    
    # 4. build matrices for heatmap
    correlation_matrix <- do.call(cbind, correlation_results_list)
    p_value_matrix <- do.call(cbind, p_value_results_list)
    
    rownames(correlation_matrix) <- genes_of_interest
    colnames(correlation_matrix) <- names(correlation_results_list)
    rownames(p_value_matrix) <- genes_of_interest
    colnames(p_value_matrix) <- names(p_value_results_list)
    
    # # - 4.1. filter to genes with >=1 significant correlation (p < 0.05) -
    # sig_genes_p_value <- rownames(p_value_matrix)[
    #   apply(p_value_matrix, 1, function(pv) any(pv < 0.05, na.rm = TRUE))
    # ]
    # 
    # # - 4.2. filter to genes with >=1 absolute correlation > 0.2 -
    # sig_genes_correlation_strength <- rownames(correlation_matrix)[
    #   apply(correlation_matrix, 1, function(cor) any(abs(cor) > 0.2, na.rm = TRUE))
    # ]
    # 
    # # Append/Combine: keep only genes that are in *both* lists (satisfy both conditions)
    # sig_genes <- intersect(sig_genes_p_value, sig_genes_correlation_strength)
    # 
    # 
    # # warn when nothing is significant:
    # if (length(sig_genes) == 0) {
    #   warning("No genes pass the significance filter for pathway: ", pathway_name, "; plotting all genes.")
    #   sig_genes <- rownames(p_value_matrix)
    # }
    # 
    # # subset both matrices
    # correlation_matrix <- correlation_matrix[sig_genes, , drop = FALSE]
    # p_value_matrix <- p_value_matrix[sig_genes, , drop = FALSE]
    
    # check if there are still genes to plot after filtering
    if (nrow(correlation_matrix) == 0) {
      warning("No genes to plot for pathway: ", pathway_name, " after filtering. Skipping heatmap generation.")
      next
    }
    
    # 5. plot heatmap with significance stars
    col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    
    ht <- Heatmap(correlation_matrix,
                  name = "Spearman ρ",
                  col = col_fun,
                  cluster_rows = TRUE,
                  cluster_columns = TRUE,
                  show_row_names = TRUE,
                  show_column_names = TRUE,
                  column_title = paste0("Spearman Correlation of ", visit, "vsV1 l2fc in", "Gene Counts with TRI_ELISA_V4\nPathway: ", pathway_name),
                  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
                  row_names_gp = gpar(fontsize = 5),
                  column_names_gp = gpar(fontsize = 10) #,
                  # cell_fun = function(j, i, x, y, width, height, fill) {
                  #   p_val <- p_value_matrix[i, j]
                  #   if (!is.na(p_val) && p_val < 0.05) {
                  #     stars <- if (p_val < 0.001) "***"
                  #     else if (p_val < 0.01) "**"
                  #     else "*"
                  #     grid.text(stars, x, y, gp = gpar(fontsize = 10, col = "black"))
                  #   }
                  # }
                  )
    
    # define filename using the M ID
    filename <- paste0("figs/correlation/l2fc_pathway_genes_spearman_corr/", visit, "_", pathway_id, "_correlation_heatmap.png")
    
    png(filename, width = 1400, height = 1600, res = 300)
    draw(ht, padding = unit(c(10, 40, 10, 10), "mm"))
    dev.off()
  }
}



