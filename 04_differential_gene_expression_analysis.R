# script to explore the transcriptomic response over time

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
library(GSVA)        # Subject-wise gene set enrichment

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

# differential gene expression analysis (DGEA) -------------------------------------

# generic function to run differential gene expression analysis (essentially a wrapper for results())
run_dgea <- function(dds,
                     conditions,
                     method = c("contrast", "contrast_list", "name"),
                     factor_name = NULL,
                     ref = NULL,
                     contrast_list = NULL,
                     name_list = NULL) {
  
  method <- match.arg(method)
  
  # initialize nested list for visits
  nested_list <- list()
  
  # loop through each condition (visit)
  for (cond in conditions) {
    
    # Get number of genes before filtering
    num_genes_before <- length(rownames(counts(dds)))
    
    # Determine the parameters for results() based on the chosen method
    if (method == "contrast") {
      if (is.null(factor_name) || is.null(ref)) {
        stop("When using method 'contrast', please supply both 'factor_name' and 'ref'.")
      }
      contrast_param <- c(factor_name, cond, ref)
      res <- results(dds, contrast = contrast_param, tidy = TRUE, independentFiltering = TRUE)
    } else if (method == "contrast_list") {
      if (is.null(contrast_list) || is.null(contrast_list[[cond]])) {
        stop("When using method 'contrast_list', please supply a contrast_list with an element for each condition.")
      }
      res <- results(dds, contrast = contrast_list[[cond]], tidy = TRUE, independentFiltering = TRUE)
    } else if (method == "name") {
      if (is.null(name_list) || is.null(name_list[[cond]])) {
        stop("When using method 'name', please supply a name_list with an element for each condition.")
      }
      res <- results(dds, name = name_list[[cond]], tidy = TRUE, independentFiltering = TRUE)
    }
    
    # Remove rows with NA p-values, these only get set to NA if all samples have 0 counts or an extreme outlier was detected via Cook's distance
    res <- res[!is.na(res$pvalue), ]
    num_genes_after <- length(res$row)
    
    # Print how many genes were removed due to being extreme outliers
    removed_genes <- num_genes_before - num_genes_after
    message(paste("Applying results() function has removed", removed_genes, "gene(s) due to being extreme outliers for", cond))
    
    # Subset significant genes using padj < 0.05 and |log2FoldChange| > 0.2
    sig <- subset(res, (padj < 0.05) & (abs(log2FoldChange) > 0.2))
    sig <- sig[order(sig$padj, decreasing = FALSE), ]
    
    # Store both full results and significant results under the visit key
    nested_list[[cond]] <- list(
      results = res,
      significant_results = sig
    )
  }
  
  return(nested_list)
}


# function to compare No. of DEGs by group and condition (e.g. Diagnosis and Visit)
plot_bar_by_group_and_condition <- function(dgea_results,
                                            groups,         # e.g., c("HC", "NSCLC", "MM", "SMM")
                                            conditions,     # e.g., VISITS[-1] or a single condition like "V1"
                                            ref_condition = NULL,  # reference condition (e.g., "V1")
                                            group_label = "Diagnosis",
                                            condition_label = "Visit",
                                            legend_colors) {  # named vector: names like "HC_Downregulated", "HC_Upregulated", etc.
  # Prepare data frame to store gene counts
  plot_data <- data.frame(GROUP = character(),
                          CONDITION = character(),
                          REGULATION = character(),
                          GENE_COUNT = integer(),
                          stringsAsFactors = FALSE)
  
  # Loop over conditions and groups to count genes
  for (cond in conditions) {
    for (group in groups) {
      if (!is.null(dgea_results[[group]][[cond]][["significant_results"]])) {
        up_count <- sum(dgea_results[[group]][[cond]][["significant_results"]]$log2FoldChange > 0)
        down_count <- sum(dgea_results[[group]][[cond]][["significant_results"]]$log2FoldChange < 0)
      } else {
        up_count <- 0
        down_count <- 0
      }
      # Append data for up- and downregulated genes
      plot_data <- rbind(plot_data, data.frame(GROUP = group,
                                               CONDITION = cond,
                                               REGULATION = "Upregulated",
                                               GENE_COUNT = up_count,
                                               stringsAsFactors = FALSE))
      plot_data <- rbind(plot_data, data.frame(GROUP = group,
                                               CONDITION = cond,
                                               REGULATION = "Downregulated",
                                               GENE_COUNT = down_count,
                                               stringsAsFactors = FALSE))
    }
  }
  
  # Create a combined factor for the fill aesthetic.
  # The factor levels are built so that for each group the order is Downregulated then Upregulated.
  plot_data$group_reg <- factor(paste0(plot_data$GROUP, "_", plot_data$REGULATION),
                                levels = unlist(lapply(groups, function(g) {
                                  c(paste0(g, "_Downregulated"), paste0(g, "_Upregulated"))
                                })))
  
  # Create the main barplot.
  p <- ggplot(plot_data, aes(x = GROUP, y = GENE_COUNT, fill = group_reg)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "Significantly Differentially Expressed Genes",
         x = group_label,
         y = "Number of Genes") +
    theme_minimal() +
    scale_fill_manual(values = legend_colors) +
    theme(legend.position = "none",    # hide the default legend
          panel.grid = element_blank(), # hide the grid
          # X-axis tick labels (e.g., "HC", "NSCLC", "MM", "SMM")
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Adjust 'size' here
          
          # Y-axis tick labels (the gene counts)
          axis.text.y = element_text(size = 12), # Add/Adjust 'size' here
          
          # X-axis title (e.g., "Diagnosis")
          axis.title.x = element_text(size = 14, face = "plain"), # Add/Adjust 'size' and 'face' here
          
          # Y-axis title (e.g., "Number of Genes")
          axis.title.y = element_text(size = 14, face = "plain"), # Add/Adjust 'size' and 'face' here
          
          # Plot Title
          plot.title = element_text(size = 16, hjust = 0.5, face = "plain")) # Adjust 'size', 'hjust', and 'face' here
  
  # Add facets if more than one condition is provided.
  if (length(conditions) > 1) {
    # If a reference condition is specified, use that; otherwise, default to the first condition.
    ref_cond <- if (!is.null(ref_condition)) ref_condition else conditions[1]
    p <- p + facet_wrap(~ CONDITION, nrow = 1,
                        labeller = as_labeller(function(x) paste0(x, " vs ", ref_cond))) +
      theme(
        strip.text = element_text(size = 12, face = "plain") # Adjust 'size' and 'face' here
      )
  }
  
  # --- Build the custom legend as a tile plot ---
  # Create a data frame similar to your original legend_df.
  legend_df <- data.frame(
    DIAGNOSIS = rep(groups, each = 2),
    REGULATION = rep(c("Down", "Up"), times = length(groups)),
    stringsAsFactors = FALSE
  )
  # Add the corresponding fill color for each combination.
  legend_df$fill <- NA
  for (i in seq_len(nrow(legend_df))) {
    grp <- legend_df$DIAGNOSIS[i]
    reg <- legend_df$REGULATION[i]
    if (reg == "Down") {
      legend_df$fill[i] <- legend_colors[paste0(grp, "_Downregulated")]
    } else {
      legend_df$fill[i] <- legend_colors[paste0(grp, "_Upregulated")]
    }
  }
  
  # Set the factor levels for DIAGNOSIS; we want the groups in the order provided.
  legend_df$DIAGNOSIS <- factor(legend_df$DIAGNOSIS, levels = groups)
  
  # Create the custom legend plot: smaller tiles with DIAGNOSIS on the y-axis in reversed order.
  legend_plot <- ggplot(legend_df, aes(x = REGULATION, y = DIAGNOSIS, fill = fill)) +
    geom_tile(color = "black", width = 0.6, height = 0.6) +
    scale_fill_identity() +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev(levels(legend_df$DIAGNOSIS))) +  # reverse y-axis order so that the first group is on top
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(face = "plain", size = 12),
          plot.margin = margin(1, 1, 1, 1),
          aspect.ratio = 0.5)
  
  # --- Combine the main plot and the custom legend ---
  final_plot <- cowplot::plot_grid(p, legend_plot, ncol = 2, rel_widths = c(1, 0.3))
  
  return(final_plot)
}



plot_dgea_variable_distribution <- function(dgea_results,
                                       group,          # e.g., "HC", "NSCLC", etc.
                                       condition,      # e.g., "V1", "V2", etc.
                                       variable_to_plot = c("padj", "l2fc"),
                                       hist_breaks = 400) {
  # Ensure a valid choice among "padj" and "l2fc"
  variable_to_plot <- match.arg(variable_to_plot)
  
  # Check if the specified group is present
  if (!group %in% names(dgea_results)) {
    stop("Group not found in dgea_results")
  }
  
  # Extract the data for the selected group and condition.
  # Here we assume that each group has a list "results" with the conditions.
  diagnosis_data <- dgea_results[[group]][[condition]][["results"]]
  
  if (variable_to_plot == "padj") {
    # Verify the column exists
    if (!"padj" %in% colnames(diagnosis_data)) {
      stop("No 'padj' column found in the data for the specified group and condition.")
    }
    vals <- diagnosis_data$padj
    main_title <- paste("padj Distribution for", group, "at", condition)
    x_label <- "Adjusted p-value (padj)"
    x_limits <- c(0.0, 0.1)
  } else if (variable_to_plot == "l2fc") {
    if (!"log2FoldChange" %in% colnames(diagnosis_data)) {
      stop("No 'log2FoldChange' column found in the data for the specified group and condition.")
    }
    vals <- diagnosis_data$log2FoldChange
    main_title <- paste("log2FoldChange Distribution for", group, "at", condition)
    x_label <- "log2FoldChange"
    x_limits <- c(-4, 4)
  }
  
  # Create the histogram
  p<- hist(vals,
       breaks = hist_breaks,
       col = "skyblue",
       border = "white",
       main = main_title,
       xlab = x_label,
       xlim = x_limits,
       ylab = "Frequency")
  
  # Add vertical threshold lines depending on the plotted variable.
  if (variable_to_plot == "padj") {
    abline(v = 0.05, col = "red", lty = 2, lwd = 2)
  } else if (variable_to_plot == "l2fc") {
    abline(v = -0.2, col = "red", lty = 2, lwd = 2)
    abline(v =  0.2, col = "red", lty = 2, lwd = 2)
  }
}

# to plot an adjusted p-value histogram for group "HC" at visit "V2":
#plot_dgea_variable_distribution(dgea_results, group = "HC", condition = "V2", variable_to_plot = "padj")


# gene set enrichment analysis (GSEA) -------------------------------------

# generic function to run gene set enrichment
run_gsea <- function(dgea_result, pathways, conditions) {
  
  # initialize list to store NES values for each condition
  nes_list <- list()
  
  # loop through each condition (e.g. visits or diagnoses)
  for (cond in conditions) {
    # access the corresponding dgea_result
    res <- as.data.frame(dgea_result[[cond]][["results"]])
    
    # rename column "row" containing gene_ids to GENE_ID
    res <- res %>%
      dplyr::rename(GENE_ID = row)
    
    # prepare results: select gene id and stat, remove missing values,
    # and for each gene take the most extreme stat value
    res2 <- res %>%
      dplyr::select(GENE_ID, stat) %>%
      na.omit() %>%
      distinct() %>%
      group_by(GENE_ID) %>%
      slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
      ungroup()
    
    # convert tibble into a named numeric vector for ranking
    ranks <- deframe(res2)
    
    # run fgsea
    fgseaRes <- fgsea(pathways = pathways, stats = ranks, nperm = 10000)
    
    # tidy the fgsea results: select and rename columns using the current condition name
    fgseaResTidy <- as_tibble(fgseaRes) %>%
      dplyr::select(pathway, NES, pval, padj) %>%
      dplyr::rename(!!sym(paste0("pval.", cond)) := pval) %>%
      dplyr::rename(!!sym(paste0("padj.", cond)) := padj) %>%
      filter(!is.na(NES)) %>%
      arrange(desc(NES))
    
    # rename NES column using the condition and store the result
    nes_list[[cond]] <- fgseaResTidy %>% 
      dplyr::rename(!!sym(paste0("NES.", cond)) := NES)
  }
  
  # merge all results to create a single data frame and rename the pathway column
  nes_df <- purrr::reduce(nes_list, full_join, by = "pathway") %>%
    dplyr::rename(PATHWAY = pathway)
  
  return(nes_df)  # Return the final cleaned NES matrix
}

# generic function to plot heatmap of gene set enrichment result
plot_gsea_heatmap <- function(gsea_result, 
                              conditions, 
                              padj_prefix = "padj.", 
                              nes_prefix = "NES.", 
                              padj_filter_threshold = 0.05,
                              selected_pathways = NULL) {
  # construct full column names for padj and NES using the condition names
  padj_cols <- paste0(padj_prefix, conditions)
  nes_cols  <- paste0(nes_prefix, conditions)
  
  # subset rows: either by explicit vector, or by padj threshold
  if (!is.null(selected_pathways)) {
    # keep only those pathways the user requested
    gsea_result_filtered <- gsea_result %>%
      filter(PATHWAY %in% selected_pathways)
    
    missing <- setdiff(selected_pathways, gsea_result_filtered$PATHWAY)
    if (length(missing) > 0) {
      warning("The following selected_pathways were not found: ",
              paste(missing, collapse = ", "))
    }
  } else {
    # original padj-based filter
    gsea_result_filtered <- gsea_result %>%
      filter_at(vars(one_of(padj_cols)), any_vars(. < padj_filter_threshold))
  }
  
  if (nrow(gsea_result_filtered) < 2) {
    message("  Skipping heatmap: Not enough pathways (", nrow(gsea_result_filtered), ") after filtering for clustering.")
    return(invisible(NULL)) # Return NULL or a dummy object instead of erroring
  }
  
  # Prepare NES matrix
  nes_matrix <- gsea_result_filtered %>%
    dplyr::select(PATHWAY, one_of(nes_cols)) %>%
    tidyr::pivot_longer(-PATHWAY, names_to = "Condition", values_to = "NES") %>%
    mutate(Condition = sub(paste0("^", nes_prefix), "", Condition)) %>%
    tidyr::pivot_wider(names_from = Condition, values_from = NES, values_fill = 0) %>%
    tibble::column_to_rownames("PATHWAY") %>%
    as.matrix()
  
  # Prepare padj matrix
  padj_matrix <- gsea_result_filtered %>%
    dplyr::select(PATHWAY, one_of(padj_cols)) %>%
    tidyr::pivot_longer(-PATHWAY, names_to = "Condition", values_to = "padj") %>%
    mutate(Condition = sub(paste0("^", padj_prefix), "", Condition)) %>%
    tidyr::pivot_wider(names_from = Condition, values_from = padj, values_fill = 0) %>%
    tibble::column_to_rownames("PATHWAY") %>%
    as.matrix()
  
  # Build significance annotation matrix
  sig_matrix <- apply(padj_matrix, c(1,2), function(x) {
    if (is.na(x))        return("")
    else if (x < 0.001)  return("***")
    else if (x < 0.01)   return("**")
    else if (x < 0.05)   return("*")
    else                  return("")
  })
  
  # Row clustering
  row_dist  <- dist(nes_matrix, method = "euclidean")
  row_clust <- hclust(row_dist, method = "ward.D2")
  row_order <- seriation::get_order(seriation::seriate(row_dist, method = "OLO"))
  row_dend  <- as.dendrogram(row_clust) %>% reorder(row_order)
  
  # Color ramp
  col_fun <- circlize::colorRamp2(
    c(min(nes_matrix, na.rm=TRUE), 0, max(nes_matrix, na.rm=TRUE)),
    c("blue","white","red")
  )
  
  # Draw heatmap
  ht <- ComplexHeatmap::Heatmap(
    nes_matrix,
    name = "NES",
    col = col_fun,
    cluster_rows    = row_dend,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    column_names_gp      = grid::gpar(fontsize = 14),
    show_row_names    = TRUE,
    row_names_side    = "left",
    row_dend_side     = "right",
    row_names_gp      = grid::gpar(fontsize = 14),
    row_names_max_width = grid::unit(2, "cm"),
    row_dend_width      = grid::unit(3, "cm"),
    heatmap_legend_param = list(
      legend_width = grid::unit(6, "cm"),
      title_gp     = grid::gpar(fontsize = 10, fontface = "bold"),
      label_gp     = grid::gpar(fontsize = 10),
      title_align  = "left"
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.text(sig_matrix[i, j], x, y,
                      gp = grid::gpar(fontsize = 10))
    }
  )
  
  return(ht)
}

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
    pivot_longer(all_of(nes_cols), names_to="Condition", values_to="NES") %>%
    mutate(Condition = sub(paste0("^", nes_prefix), "", Condition))
  df_padj <- df %>%
    pivot_longer(all_of(padj_cols), names_to="Condition", values_to="padj") %>%
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


plot_gsea_kinetics <- function(gsea_results,    # named list: one df per Diagnosis
                               visits,          # e.g. c("V1","V2","V3","V4")
                               padj_prefix     = "padj.",
                               nes_prefix      = "NES.",
                               padj_threshold  = 0.001,
                               cluster         = FALSE,
                               pathway2genes   = NULL,
                               cluster_method  = "average",
                               cluster_height  = 0.5) {
  
  # detect single vs multiple visits
  is_single_visit <- length(visits) == 1
  visit_label <- visits
  
  # 1) STACK across diagnoses → long table with Diagnosis, Visit, NES, padj
  df_all <- bind_rows(
    lapply(names(gsea_results), function(diag) {
      df <- gsea_results[[diag]]
      
      # pivot NES
      nes_long <- df %>%
        pivot_longer(
          cols      = matches(paste0("^", nes_prefix)),
          names_to  = "Visit",
          values_to = "NES"
        ) %>%
        mutate(
          Diagnosis = diag,
          Visit     = sub(paste0("^", nes_prefix), "", Visit)
        ) %>%
        filter(Visit %in% visits)
      
      # pivot padj
      padj_long <- df %>%
        pivot_longer(
          cols      = matches(paste0("^", padj_prefix)),
          names_to  = "Visit",
          values_to = "padj"
        ) %>%
        mutate(
          Diagnosis = diag,
          Visit     = sub(paste0("^", padj_prefix), "", Visit)
        ) %>%
        filter(Visit %in% visits)
      
      left_join(nes_long, padj_long,
                by = c("PATHWAY", "Diagnosis", "Visit"))
    })
  )
  
  # 2) FILTER: keep pathways with padj < threshold in any diag×visit, filter TBA
  keep <- df_all %>%
    group_by(PATHWAY) %>%
    dplyr::summarize(min_p = min(padj, na.rm = TRUE)) %>%
    filter(min_p < padj_threshold) %>%
    filter(!grepl("TBA", PATHWAY, ignore.case = TRUE)) %>%
    pull(PATHWAY)
  df_all <- df_all %>% filter(PATHWAY %in% keep)
  
  # 3) OPTIONAL CLUSTER + REPRESENTATIVE SELECTION (based on pathway similarity)
  if (cluster) {
    if (is.null(pathway2genes))
      stop("Need pathway2genes (named list) for clustering.")
    
    pgs       <- pathway2genes[keep]
    all_genes <- unique(unlist(pgs))
    mat       <- sapply(pgs, function(g) as.numeric(all_genes %in% g))
    rownames(mat) <- all_genes
    inter   <- crossprod(mat)
    sizes   <- colSums(mat)
    union   <- outer(sizes, sizes, "+") - inter
    sim     <- inter / union
    hc      <- hclust(as.dist(1 - sim), method = cluster_method)
    cls     <- cutree(hc, h = cluster_height)
    
    min_pth <- df_all %>%
      group_by(PATHWAY) %>%
      dplyr::summarize(min_p = min(padj, na.rm = TRUE))
    
    reps <- tibble(PATHWAY = names(cls), cluster = cls) %>%
      inner_join(min_pth, by = "PATHWAY") %>%
      group_by(cluster) %>%
      slice_min(min_p, with_ties = FALSE) %>%
      pull(PATHWAY)
    
    df_all    <- df_all %>% filter(PATHWAY %in% reps)
    ordered   <- hc$labels[hc$order]
    row_order_initial <- intersect(ordered, reps) # Store this for initial ordering
  } else {
    row_order_initial <- unique(df_all$PATHWAY) # Initial order if no clustering
  }
  
  # 4) SIGNIFICANCE BINNING + factor levels (applying initial ordering)
  df_all <- df_all %>%
    mutate(
      sig_level = case_when(
        padj < padj_threshold ~ paste0("padj<", padj_threshold),
        # padj < 0.01  ~ "padj<0.01",
        # padj < 0.05  ~ "padj<0.05",
        TRUE         ~ "ns"
      ),
      PATHWAY   = factor(PATHWAY, levels = rev(row_order_initial)), # Apply initial order here
      Visit     = factor(Visit, levels = visits),
      Diagnosis = factor(Diagnosis, levels = names(gsea_results))
    )
  
  # 5) HIERARCHICAL CLUSTERING OF ROWS BASED ON NES (for final plot ordering)
  # Reshape df_all to wide format for clustering NES values
  df_nes_wide <- df_all %>%
    unite("Diagnosis_Visit", Diagnosis, Visit, sep = "_", remove = FALSE) %>%
    dplyr::select(PATHWAY, Diagnosis_Visit, NES) %>%
    pivot_wider(names_from = Diagnosis_Visit, values_from = NES)
  
  # Prepare matrix for hclust
  # Ensure PATHWAY column is excluded and set as rownames
  nes_mat <- as.matrix(df_nes_wide[, -1]) 
  rownames(nes_mat) <- df_nes_wide$PATHWAY
  
  # Handle NA values (e.g., pathways not present in all conditions)
  # Replacing with 0 is a simple approach; consider imputation for more complex scenarios.
  nes_mat[is.na(nes_mat)] <- 0 
  
  dist_nes <- dist(nes_mat)
  hc_nes   <- hclust(dist_nes, method = "ward.D2")
  row_order_final <- hc_nes$labels[hc_nes$order] # This will be the final order
  
  # 6) APPLY FINAL ROW ORDER
  df_all$PATHWAY <- factor(df_all$PATHWAY, levels = rev(row_order_final))
  
  # Define a list of theme elements that are applied only for a single visit
  single_visit_specific_theme <- list()
  if (is_single_visit) {
    single_visit_specific_theme <- list(
      panel.spacing = unit(0.5, "lines"),      # Conditionally set panel spacing
      axis.text.x = element_blank(),         # Hide individual Visit labels
      axis.ticks.x = element_blank(),        # Hide individual Visit ticks
      axis.title.x = element_text(vjust = -1.5) # Style the overall Visit label
    )
  }
  
  # 7) PLOT: visits on x, pathways on y, dodge by Diagnosis
  p <- ggplot(df_all,
              aes(x = Visit, y = PATHWAY,
                  color = NES,
                  size  = -log10(padj),
                  shape = sig_level,
                  group = Diagnosis)) +
    geom_point(position = position_dodge(width = 0.7)) +
    scale_color_gradient2(
      low      = "blue", mid = "white", high = "red",
      midpoint = 0, name = "NES"
    ) +
    scale_size_continuous(name = expression(-log[10](padj))) +
    scale_shape_manual(
      values = setNames(
        c(1, 16),
        c("ns", paste0("padj<", padj_threshold))
      ),
      breaks = c("ns", paste0("padj<", padj_threshold)),
      name   = "Significance"
    ) +
    labs(x = if(length(visits) == 1) visits else NULL, y = NULL) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title  = element_text(hjust = 0.8)
    ) +
    # use ggh4x facet_grid2 with themed strips
    ggh4x::facet_grid2(
      . ~ Diagnosis,
      scales = "free_x",
      space  = "free_x",
      switch = "x",
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = DIAGNOSIS_COLKEY
        )
      )
    ) +
    do.call(theme, single_visit_specific_theme)
  
  print(p)
}



# prepare data ---------------------------------------------------------------

# read in dds
dds <- readRDS("output/dds_qc.rds")

metadata <- as.data.frame(colData(dds))

count_matrix <- counts(dds, normalized = TRUE)

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

# define V1 subset
V1_dds <- subset_dds(dds, condition = "VISIT == 'V1'")

# set to TRUE for debugging (much faster)
V1_vsd <- as.data.frame(assay(vst(V1_dds, blind = TRUE)))

V1_metadata <- as.data.frame(colData(V1_dds))

V1_count_matrix <- as.data.frame(counts(V1_dds, normalized = TRUE))

pathway_data <- readRDS("output/pathway_data.rds")

V1_pathway_data <- subset_pathway_data(pathway_data, condition = "VISIT == 'V1'")

l2fc_pathway_data <- readRDS("output/l2fc_pathway_data.rds")

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

# DESeq2 analysis for diagnosis subsets ~ STUDY_ID + VISIT ---------------------------------------------------------------

# set design
design(dds) <- ~ STUDY_ID + VISIT

# loop over DIAGNOSES
for (diagnosis in DIAGNOSES) {
  
  message(paste("Analysing", diagnosis, "kinetics."))
  
  dds_subset <- subset_dds(dds, condition = paste0("DIAGNOSIS == ", "'", diagnosis, "'"))
  
  # check if DESeq2 run exists
    dds_file <- paste0("output/", diagnosis, "_V1_contrast_dds.rds")
    if (file.exists(dds_file)) {
      # load existing results
      dds_subset <- readRDS(dds_file)
      message(paste("Loaded existing DESeq2 results for", diagnosis, "from file."))
    } else {
      message(paste("Running DESeq2 for", diagnosis))
      # run DESeq2
      dds_subset <- DESeq(dds_subset, parallel = TRUE, BPPARAM=MulticoreParam(2)) # two cores works quite well, more could have diminishing returns
      # save dds
      saveRDS(dds_subset, file = dds_file)
    }
  
  # get degs and store in list
  dgea_results[[diagnosis]] <- run_dgea(dds_subset, conditions = VISITS[-1], method = "contrast", factor_name = "VISIT", ref = "V1")
  
  # check if fgsea run exists
  nes_file <- paste0("output/", diagnosis, "_nes.rds")
  if (file.exists(nes_file)) {
    # load existing results
    gsea_results[[diagnosis]] <- readRDS(nes_file)
    message(paste("Loaded existing fgsea results for", diagnosis, "from file."))
  } else {
    message(paste("Running fgsea for", diagnosis))
    # get NE scores and store in list
    gsea_results[[diagnosis]] <- run_gsea(dgea_results[[diagnosis]], PATHWAYS, conditions = VISITS[-1])
    # save nes
    saveRDS(gsea_results[[diagnosis]], file = nes_file)
  }
  
  # ht <- plot_gsea_heatmap(gsea_results[[diagnosis]], conditions = VISITS[-1])
  # 
  # # save heatmap
  # png(filename = paste0("figs/analysis/gsea/", diagnosis, "_gsea_heatmap_by_VISIT.png"), width = 1000, height = 3600)
  # draw(ht,
  #      padding = unit(c(10, 80, 10, 20), "mm"),
  #      column_title = paste(diagnosis, "GSEA Heatmap by VISIT"),
  #      column_title_gp = gpar(fontsize = 16, fontface = "bold"))
  # dev.off()
  
}

# generate legend colors with desaturated downregulated versions
legend_colors <- unlist(lapply(names(DIAGNOSIS_COLKEY), function(diagnosis) {
  base_color <- DIAGNOSIS_COLKEY[[diagnosis]]
  # desaturate base color
  down_color <- hex(colorspace::desaturate(hex2RGB(base_color), 0.5))
  
  # create named vector with up/down colors
  c(
    setNames(down_color, paste0(diagnosis, "_Downregulated")),
    setNames(base_color, paste0(diagnosis, "_Upregulated"))
  )
}))

# plot differentially expressed genes by VISIT and DIAGNOSIS
deg_barplot <- plot_bar_by_group_and_condition(dgea_results,
                                               groups = DIAGNOSES,
                                               conditions = VISITS[-1],
                                               ref_condition = "V1",
                                               group_label = "Diagnosis",
                                               condition_label = "Visit",
                                               legend_colors = legend_colors)

ggsave("figs/analysis/dgea/COHORT_degs_bar_by_DIAGNOSIS_and_VISIT.png", plot = deg_barplot, width = 12, height = 8, dpi = 300)


dm <- plot_gsea_kinetics(
  gsea_results    = gsea_results,
  visits          = VISITS[2:3],
  padj_threshold  = 0.05,
  cluster         = TRUE,
  pathway2genes   = PATHWAYS,
  cluster_height  = 0.99
) + ggtitle("Pathway kinetics across visits & diagnoses")

ggsave(paste0("figs/analysis/gsea/combined_kinetics_dotmatrix.png"), plot = dm, width = 8, height = 5, dpi = 300)

for (visit in VISITS[-1]) {
  dm <- plot_gsea_kinetics(
    gsea_results    = gsea_results,
    visit          = visit,
    padj_threshold  = 0.05,
    cluster         = TRUE,
    pathway2genes   = PATHWAYS,
    cluster_height  = 0.99
  ) + force_panelsizes(cols = unit(1.5, "cm")) + ggtitle("Pathway kinetics across visits & diagnoses")
  
  ggsave(paste0("figs/analysis/gsea/", visit, "_kinetics_dotmatrix.png"), plot = dm, width = 8, height = 5, dpi = 300)
}

# # investigate DEG results with goterms
# dgea_result <- dgea_results[["SMM"]][["V4"]][["significant_results"]]
# 
# # Order data based on criterion
# ordered_data <- dgea_result[order(dgea_result[["padj"]], decreasing = FALSE), ]
# 
# # Select the top genes based on ordering
# top_genes <- head(ordered_data, 100)
# top_gene_names <- top_genes[["row"]]
# 
# ego <- dotplot_enrichment(top_gene_names)
# ego$plot

# plot DGEA variable histograms
# this does not save the plots and I have no idea why... you can however save it from R studio
for (var in c("padj", "l2fc")) {
  for (visit in VISITS[-1]) {
    for (diag in DIAGNOSES) {
      png(filename = paste0("figs/analysis/dgea/deseq_res_histogram/", visit, "_", diag,
                            "_", var, "_histogram.png"), width = 1800, height = 1800, res = 300)
      print(plot_dgea_variable_distribution(dgea_results, group = diag, condition = visit, variable_to_plot = var))
      dev.off()
    }
  }
}


# DESeq2 analysis for visit subsets ~ SEX + AGE + DIAGNOSIS ---------------------------------------------------------------

# set design
design(dds) <- ~ SEX + AGE + DIAGNOSIS

# loop over VISITS
for (visit in VISITS) {
  
  message(paste("Analysing", visit, "diagnosis differences."))
  
  # subset dds
  dds_subset <- dds[, dds$VISIT == visit]
  
  # make sure to reset the factor levels, else we run into an error with the current design
  colData(dds_subset) <- droplevels(colData(dds_subset))
  
  # check if DESeq2 run exists
  dds_file <- paste0("output/", visit, "_HC_contrast_dds.rds")
  if (file.exists(dds_file)) {
    # load existing results
    dds_subset <- readRDS(dds_file)
    message(paste("Loaded existing DESeq2 results for", visit, "from file."))
  } else {
    message(paste("Running DESeq2 for", visit))
    # run DESeq2
    dds_subset <- DESeq(dds_subset, parallel = TRUE, BPPARAM=MulticoreParam(2)) # two cores works quite well, more could have diminishing returns
    # save dds
    saveRDS(dds_subset, file = dds_file)
  }
  
  # get degs and store in list
  dgea_results[[visit]] <- run_dgea(dds_subset, conditions = c("SMM", "MM", "NSCLC"), method = "contrast", factor_name = "DIAGNOSIS", ref = "HC")
  
  # check if fgsea run exists
  nes_file <- paste0("output/", visit, "_nes.rds")
  if (file.exists(nes_file)) {
    # load existing results
    gsea_results[[visit]] <- readRDS(nes_file)
    message(paste("Loaded existing fgsea results for", visit, "from file."))
  } else {
    message(paste("Running fgsea for", visit))
    # get NE scores and store in list
    gsea_results[[visit]] <- run_gsea(dgea_results[[visit]], PATHWAYS, conditions = c("SMM", "MM", "NSCLC"))
    # save nes
    saveRDS(gsea_results[[visit]], file = nes_file)
  }
  
  # # plot GSEA heatmap
  # ht <- plot_gsea_heatmap(gsea_results[[visit]], conditions = c("SMM", "MM", "NSCLC"))
  # 
  # # save heatmap
  # png(filename = paste0("figs/analysis/gsea/", visit, "_gsea_heatmap_by_DIAGNOSIS.png"), width = 1000, height = 3600)
  # draw(ht,
  #      padding = unit(c(10, 80, 10, 20), "mm"),
  #      column_title = paste(visit, "GSEA Heatmap by DIAGNOSIS"),
  #      column_title_gp = gpar(fontsize = 16, fontface = "bold"))
  # dev.off()
  
  # plot dotmatrix
  dm <- plot_gsea_dotmatrix(gsea_results[[visit]],
                            conditions = c("SMM", "MM", "NSCLC"),
                            padj_threshold = 0.01,
                            pathway2genes = PATHWAYS,
                            cluster = TRUE,
                            cluster_height = 0.96) +
    force_panelsizes(cols = unit(1.5, "cm")) + ggtitle(paste0("Diagnosis comparison at ", visit, " with healthy controls"))
  
  ggsave(paste0("figs/analysis/gsea/", visit, "_gsea_dotmatrix_by_DIAGNOSIS.png"), plot = dm, width = 8, height = 5, dpi = 300)
}

# generate legend colors with desaturated downregulated versions
legend_colors <- unlist(lapply(names(VISIT_COLKEY), function(visit) {
  base_color <- VISIT_COLKEY[[visit]]
  # desaturate base color
  down_color <- hex(colorspace::desaturate(hex2RGB(base_color), 0.5))
  
  # create named vector with up/down colors
  c(
    setNames(down_color, paste0(visit, "_Downregulated")),
    setNames(base_color, paste0(visit, "_Upregulated"))
  )
}))


deg_barplot <- plot_bar_by_group_and_condition(dgea_results,
                                groups = VISITS,
                                conditions = DIAGNOSES[-1],
                                ref_condition = "HC",
                                group_label = "Visit",
                                condition_label = "Diagnosis",
                                legend_colors = legend_colors)

ggsave("figs/analysis/dgea/COHORT_degs_bar_by_VISIT_and_DIAGNOSIS.png", plot = deg_barplot, width = 12, height = 8, dpi = 300)


# subject l2fc heatmap ---------------------------------------

# initialize a list to store clusters by visit
stored_clusters <- list()

plot_unique_degs_subject_heatmap <- function(l2fc_data,
                                             dgea_results,
                                             diagnoses,                   # Either a vector of diagnoses (e.g., c("HC", "NSCLC", "MM", "SMM")) or a single diagnosis (e.g., "NSCLC")
                                             visit,                       # e.g., "V2"
                                             variable_to_rank = "pvalue", # e.g., "padj", "pvalue", "log2FoldChange"
                                             top_n = 50,                  # Number of top genes to select per diagnosis
                                             scale_rows = FALSE,          # Whether genes will be Z-scored for better comparisons
                                             cluster_by = "rows",         # Accepts: "rows", "columns", or "both"
                                             n_row_clusters = 4,          # Number of row clusters
                                             n_column_clusters = 4,       # Number of column clusters
                                             additional_annotation_col = NULL  # Name of an extra metadata column to annotate (e.g., "BATCH")
) {
  # Required packages
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  if (!is.null(additional_annotation_col)) {
    library(RColorBrewer)
  }
  
  metadata <- l2fc_data$metadata
  count_matrix <- l2fc_data$count_matrix
  
  
  # Select Top Genes by variable for each diagnosis -------------------------
  all_top_genes <- character(0)
  
  for (diag in diagnoses) {
    if (!diag %in% names(dgea_results)) {
      warning(paste("Diagnosis", diag, "not found in dgea_results"))
      next
    }
    
    # Extract the full results for this diagnosis and visit
    diag_res <- dgea_results[[diag]][[visit]][["results"]]
    
    if (is.null(diag_res)) {
      warning(paste("No results found for diagnosis", diag, "at", visit))
      next
    }
    
    if (!variable_to_rank %in% colnames(diag_res)) {
      warning(paste("No", variable_to_rank, "column found in the results for diagnosis", diag, "at", visit))
      next
    }
    
    # Order the results by adjusted p-value (lowest first)
    ordered_res <- diag_res[order(diag_res[[variable_to_rank]], decreasing = FALSE), ]
    
    # Select the top_n genes and store the gene names
    top_genes <- head(ordered_res$row, top_n)
    all_top_genes <- c(all_top_genes, top_genes)
  }
  
  # Obtain the unique union of top genes
  unique_genes <- unique(all_top_genes)
  message("Number of unique top genes: ", length(unique_genes))
  
  # Subset the count matrix to include only the unique genes
  count_matrix <- count_matrix[rownames(count_matrix) %in% unique_genes, , drop = FALSE]
  
  
  
  # Filter Metadata and Subset Count Matrix ---------------------------------
  metadata_filtered <- metadata %>% 
    filter(VISIT == visit) %>%
    filter(DIAGNOSIS %in% diagnoses) %>%
    droplevels() %>%                      # This drops unused factor levels.
    mutate(CEGAT_ID = as.character(CEGAT_ID),
           STUDY_ID = as.character(STUDY_ID)) %>% # Converting factors storing IDs to characters here helps to keep them in sync
    arrange(DIAGNOSIS, match(CEGAT_ID, colnames(count_matrix)))
  
  count_matrix <- count_matrix[, metadata_filtered$CEGAT_ID, drop = FALSE]
  
  
  # Remove genes with NA or near-constant variation -------------------------
  row_sd <- apply(count_matrix, 1, function(x) sd(x, na.rm = TRUE))
  
  if (any(is.na(row_sd))) {
    message("Removing ", sum(is.na(row_sd)), " gene(s) with NA variance.")
    count_matrix <- count_matrix[!is.na(row_sd), , drop = FALSE]
    row_sd <- row_sd[!is.na(row_sd)]
  }
  
  message("Minimum SD before filtering: ", min(row_sd))
  epsilon <- 1e-8
  if (any(row_sd < epsilon, na.rm = TRUE)) {
    message("Removing ", sum(row_sd < epsilon, na.rm = TRUE), " gene(s) with zero or near-zero variance.")
    count_matrix <- count_matrix[row_sd > epsilon, , drop = FALSE]
  }
  
  # Create a named vector mapping from CEGAT_IDs to STUDY_IDs
  id_mapping <- setNames(as.character(metadata_filtered$STUDY_ID), metadata_filtered$CEGAT_ID)
  
  # Reorder the STUDY_IDs to match the order of the columns in count_matrix
  new_colnames <- id_mapping[colnames(count_matrix)]
  
  # Apply the new column names to count_matrix
  colnames(count_matrix) <- new_colnames
  
  
  # Row Scaling (Normalization) ---------------------------------------------
  if (scale_rows == TRUE) {
    count_matrix <- t(scale(t(count_matrix), center = TRUE, scale = TRUE))
  }
  
  
  
  # Build the Top Annotation Data Frame -------------------------------------
  ann_df <- data.frame(Visit = rep(visit, nrow(metadata_filtered)),
                       stringsAsFactors = FALSE)
  
  # Convert Diagnosis column to a factor using the filtered set levels
  if (length(diagnoses) == 1) {
    ann_df$Diagnosis <- factor(rep(diagnoses[1], nrow(metadata_filtered)), levels = diagnoses)
  } else {
    ann_df$Diagnosis <- factor(metadata_filtered$DIAGNOSIS, levels = unique(metadata_filtered$DIAGNOSIS))
  }
  
  if (!is.null(additional_annotation_col)) {
    if (!(additional_annotation_col %in% colnames(metadata_filtered))) {
      stop(paste("Column", additional_annotation_col, "not found in metadata."))
    }
    ann_df[[additional_annotation_col]] <- metadata_filtered[[additional_annotation_col]]
  }
  
  
  # Build Color Mapping for Annotations -------------------------------------
  
  ann_colors <- list()
  
  ann_colors[["Visit"]] <- VISIT_COLKEY
  
  ann_colors[["Diagnosis"]] <- DIAGNOSIS_COLKEY
  
  if (!is.null(additional_annotation_col)) {
    extra_vals <- ann_df[[additional_annotation_col]]
    if (is.numeric(extra_vals)) {
      ann_colors[[additional_annotation_col]] <- colorRamp2(c(min(extra_vals, na.rm = TRUE),
                                                              max(extra_vals, na.rm = TRUE)),
                                                            c("blue", "red"))
    } else {
      extra_levels <- unique(as.character(extra_vals))
      n_levels <- length(extra_levels)
      if(n_levels < 3){
        extra_pal <- brewer.pal(n = 3, name = "Set1")[seq_len(n_levels)]
      } else {
        extra_pal <- brewer.pal(n = n_levels, name = "Set1")
      }
      ann_colors[[additional_annotation_col]] <- setNames(extra_pal, extra_levels)
    }
  }
  
  top_annotation <- HeatmapAnnotation(
    df = ann_df,
    col = ann_colors,
    annotation_name_gp = gpar(fontsize = 10),
    annotation_legend_param = list(
      Diagnosis = list(title = "Diagnosis", at = levels(ann_df$Diagnosis), labels = levels(ann_df$Diagnosis))
    )
  )
  
  
  # Setup the Heatmap Color Function ----------------------------------------
  
  # maximum and minimum l2fc boundaries
  # col_fun <- colorRamp2(
  #   c(min(count_matrix, na.rm = TRUE), 0, max(count_matrix, na.rm = TRUE)),
  #   c("blue", "white", "red")
  # )
  
  # -3 and +3 l2fc boundaries
  col_fun <- colorRamp2(
    c(-3, 0, 3),
    c("blue", "white", "red")
  )
  
  
  # Clustering and Heatmap Creation with Split Clusters ---------------------
  
  ht <- NULL
  clusters <- list()
  
  if (cluster_by == "rows") {
    row_dist <- dist(count_matrix, method = "euclidean")
    row_clust <- hclust(row_dist, method = "ward.D2")
    row_order <- seriation::get_order(seriation::seriate(row_dist, method = "OLO"))
    clusters_row <- cutree(row_clust, k = n_row_clusters)
    
    ht <- Heatmap(
      count_matrix,
      name = "log2FoldChange",
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
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 6),
      top_annotation = top_annotation,
      row_dend_side = "right",
      row_dend_width = unit(3, "cm"),
      show_parent_dend_line = TRUE,
      use_raster = TRUE
    )
    clusters$row_clusters <- clusters_row
    
  } else if (cluster_by == "columns") {
    col_dist <- dist(t(count_matrix), method = "euclidean")
    col_clust <- hclust(col_dist, method = "ward.D2")
    col_order <- seriation::get_order(seriation::seriate(col_dist, method = "OLO"))
    clusters_col <- cutree(col_clust, k = n_column_clusters)
    
    ht <- Heatmap(
      count_matrix,
      name = "log2FoldChange",
      col = col_fun,
      cluster_rows = FALSE,
      
      cluster_columns = TRUE,
      clustering_distance_columns = function(m) dist(m, method = "euclidean"),
      clustering_method_columns = "ward.D2",
      column_order = col_order,
      column_split = clusters_col,     # split columns into n_column_clusters
      column_dend_reorder = FALSE,
      
      show_column_names = TRUE,
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 8),
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
    row_dist <- dist(count_matrix, method = "euclidean")
    row_clust <- hclust(row_dist, method = "ward.D2")
    row_order <- seriation::get_order(seriation::seriate(row_dist, method = "OLO"))
    clusters_row <- cutree(row_clust, k = n_row_clusters)
    
    # Column clustering
    col_dist <- dist(t(count_matrix), method = "euclidean")
    col_clust <- hclust(col_dist, method = "ward.D2")
    col_order <- seriation::get_order(seriation::seriate(col_dist, method = "OLO"))
    clusters_col <- cutree(col_clust, k = n_column_clusters)
    
    ht <- Heatmap(
      count_matrix,
      name = "log2FoldChange",
      col = col_fun,
      
      cluster_rows = TRUE, # perform clustering internally
      clustering_distance_rows = function(m) dist(m, method = "euclidean"),
      clustering_method_rows = "ward.D2",
      row_order = row_order,
      row_split = clusters_row,        # split rows into n_row_clusters
      row_dend_reorder = FALSE,
      
      cluster_columns = TRUE,
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

# produce single heatmap
subject_heatmap <- plot_unique_degs_subject_heatmap(l2fc_data_imputed, dgea_results,
                                                    diagnoses = "MM",
                                                    visit = "V3",
                                                    top_n = 50,
                                                    scale_rows = FALSE,
                                                    cluster_by = "both",
                                                    n_row_clusters = 3,
                                                    n_column_clusters = 3,
                                                    additional_annotation_col = "LENALIDOMIDE_AT_V1")

# inspect the heatmap and adjust number of clusters
subject_heatmap$ht

png(paste0("figs/analysis/subject_l2fc_heatmap/custom_subject_heatmap.png"), width = 4400, height = 3600, res = 300)
draw(subject_heatmap$ht)
dev.off()

# define number of row and column clusters
heatmap_clusters <- list(
  "V2" = list(n_row_clusters = 3,
              n_column_clusters = 3),
  "V3" = list(n_row_clusters = 4,
              n_column_clusters = 3),
  "V4" = list(n_row_clusters = 4,
              n_column_clusters = 4),
  "V5" = list(n_row_clusters = 4,
              n_column_clusters = 8),
  "V6" = list(n_row_clusters = 5,
              n_column_clusters = 5)
)


# loop over visits to create relevant plots
for (visit in VISITS[-1]) {
  
  n_row_clusters <- heatmap_clusters[[visit]]$n_row_clusters
  n_column_clusters <- heatmap_clusters[[visit]]$n_column_clusters
  
  # create heatmap
  subject_heatmap <- plot_unique_degs_subject_heatmap(l2fc_data_imputed, dgea_results,
                                                      diagnoses = DIAGNOSES,
                                                      visit = visit,
                                                      top_n = 50,
                                                      scale_rows = FALSE,
                                                      cluster_by = "both",
                                                      n_row_clusters = n_row_clusters,
                                                      n_column_clusters = n_column_clusters,
                                                      additional_annotation_col = NULL)
  
  # save the heatmap as a PNG file
  png(paste0("figs/analysis/subject_l2fc_heatmap/", visit, "_subject_heatmap.png"), width = 4400, height = 3600, res = 300)
  draw(subject_heatmap$ht)
  dev.off()
  
  # store row and column clusters for the current visit
  row_clusters <- subject_heatmap$clusters$row_clusters
  column_clusters <- subject_heatmap$clusters$column_clusters
  
  # get unique row cluster ids
  row_cluster_ids <- sort(unique(row_clusters))
  
  # loop through each gene cluster
  for (cluster_id in row_cluster_ids) {
    
    # get the gene names for the current cluster
    cluster_gene_set <- names(row_clusters)[
      row_clusters == cluster_id
    ]
    
    # generate the GO term dotplot
    plot <- dotplot_enrichment(cluster_gene_set, gmt_file = BTM_FILE)$plot
    
    if (!is.null(plot)) {
      
      plot <- plot + ggtitle(paste("Gene Ontology Terms enriched in gene-cluster", cluster_id, "at", visit))
      
      # save to file
      ggsave(
        filename = paste0("figs/analysis/subject_l2fc_heatmap/", visit, "_subject_heatmap_gene_cluster_", cluster_id, "_enrichment.png"),
        plot = plot,
        width = 8,
        height = 6
      )
    }
  }
  
  # get total samples at this visit (all clusters)
  metadata_visit <- l2fc_data_imputed$metadata %>% dplyr::filter(VISIT == visit)
  total_visit_samples <- nrow(metadata_visit)
  
  # assemble cluster assignments with metadata
  cluster_df <- data.frame(
    STUDY_ID = names(column_clusters),
    cluster = as.factor(column_clusters)
  ) %>%
    left_join(l2fc_data_imputed$metadata, by = c("STUDY_ID" = "STUDY_ID")) %>%
    dplyr::filter(VISIT == visit)
  
  # count per cluster and diagnosis, percent relative to all visit samples
  visit_diag_cluster <- cluster_df %>%
    dplyr::count(cluster, DIAGNOSIS) %>%
    mutate(percentage_of_visit = 100 * n / total_visit_samples)
  
  # plot: x = cluster, y = % of visit, fill = diagnosis
  p_visit <- ggplot(visit_diag_cluster, aes(x = cluster, y = percentage_of_visit, fill = DIAGNOSIS)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_manual(values = DIAGNOSIS_COLKEY) +
    labs(
      title = paste0("Diagnosis proportions across clusters at ", visit),
      x = "Cluster ID",
      y = "Percentage of all samples (%)",
      fill = "Diagnosis"
    ) +
    theme_minimal(base_size = 14)
  
  # save the visit-level cluster distribution plot
  ggsave(
    filename = paste0("figs/analysis/subject_l2fc_heatmap/", visit, "_diagnoses_percentage_by_cluster.png"),
    plot = p_visit,
    width = 8,
    height = 6
  )
}



# plot specific gene sets
# helper function, retrieves unique set of genes from pathways
extract_gene_set_from_pathways <- function(pathways, pathway_names) {
  
  # Filter pathways to only those with a Name in pathway_names
  selected <- pathways[pathway_names]
  
  # Unlist the gene names from the "Value" column and remove duplicates
  gene_set <- unique(unlist(selected))
  
  return(gene_set)
}

gene_set <- extract_gene_set_from_pathways(PATHWAYS, pathway_names = c("plasma cells & B cells, immunoglobulins (M156.0)",
                                                                       "plasma cells, immunoglobulins (M156.1)",
                                                                       "Plasma cell surface signature (S3)"))

# fixed set of genes heatmap
plot_genes_subject_heatmap <- function(l2fc_data,
                                       diagnoses,                  # Either a vector of diagnoses or a single diagnosis (e.g., "NSCLC")
                                       visit,                      # e.g., "V2"
                                       gene_set,                   # A fixed vector of genes to plot (e.g., c("CD27", "CYAT1", "DERL3"))
                                       scale_rows = FALSE,          # Whether to Z-score gene rows for better comparison
                                       cluster_by = "rows",        # Accepts: "rows", "columns", or "both"
                                       n_row_clusters = 4,         # Number of row clusters
                                       n_column_clusters = 4,      # Number of column clusters
                                       additional_annotation_col = NULL  # Name of an extra metadata column to annotate (e.g., "AGE")
) {
  # Required packages
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  if (!is.null(additional_annotation_col)) {
    library(RColorBrewer)
  }
  
  metadata <- l2fc_data$metadata
  count_matrix <- l2fc_data$count_matrix
  
  message("Number of genes provided: ", length(gene_set))
  
  # Subset the count matrix to include only those genes in gene_set
  count_matrix <- count_matrix[rownames(count_matrix) %in% gene_set, , drop = FALSE]
  
  # Filter metadata for the specified visit and diagnoses, and order by DIAGNOSIS then by CEGAT_ID.
  metadata_filtered <- metadata %>% 
    filter(VISIT == visit) %>%
    filter(DIAGNOSIS %in% diagnoses) %>%
    droplevels() %>%                      # This drops unused factor levels.
    mutate(CEGAT_ID = as.character(CEGAT_ID),
           STUDY_ID = as.character(STUDY_ID)) %>% # Converting factors storing IDs to characters here helps to keep them in sync
    arrange(DIAGNOSIS, match(CEGAT_ID, colnames(count_matrix)))
  
  # Subset count_matrix to include only columns corresponding to filtered metadata
  count_matrix <- count_matrix[, metadata_filtered$CEGAT_ID, drop = FALSE]
  
  ##### Remove genes with NA or near-constant variation #####
  row_sd <- apply(count_matrix, 1, function(x) sd(x, na.rm = TRUE))
  
  if (any(is.na(row_sd))) {
    message("Removing ", sum(is.na(row_sd)), " gene(s) with NA variance.")
    count_matrix <- count_matrix[!is.na(row_sd), , drop = FALSE]
    row_sd <- row_sd[!is.na(row_sd)]
  }
  
  message("Minimum SD before filtering: ", min(row_sd, na.rm = TRUE))
  epsilon <- 1e-8
  if (any(row_sd < epsilon, na.rm = TRUE)) {
    message("Removing ", sum(row_sd < epsilon, na.rm = TRUE), " gene(s) with zero or near-zero variance.")
    count_matrix <- count_matrix[row_sd > epsilon, , drop = FALSE]
  }
  
  # Create a named vector mapping from CEGAT_IDs to STUDY_IDs
  id_mapping <- setNames(as.character(metadata_filtered$STUDY_ID), metadata_filtered$CEGAT_ID)
  
  # Reorder the STUDY_IDs to match the order of the columns in count_matrix
  new_colnames <- id_mapping[colnames(count_matrix)]
  
  # Apply the new column names to count_matrix
  colnames(count_matrix) <- new_colnames
  
  # Scale each row (gene normalization)
  if (scale_rows == TRUE) {
    count_matrix <- t(scale(t(count_matrix), center = TRUE, scale = TRUE))
  }
  
  ##### Build the top annotation data frame #####
  ann_df <- data.frame(Visit = rep(visit, nrow(metadata_filtered)),
                       stringsAsFactors = FALSE)
  
  if (length(diagnoses) == 1) {
    ann_df$Diagnosis <- factor(rep(diagnoses[1], nrow(metadata_filtered)), levels = diagnoses)
  } else {
    ann_df$Diagnosis <- factor(metadata_filtered$DIAGNOSIS, levels = unique(metadata_filtered$DIAGNOSIS))
  }
  
  if (!is.null(additional_annotation_col)) {
    if (!(additional_annotation_col %in% colnames(metadata_filtered))) {
      stop(paste("Column", additional_annotation_col, "not found in metadata."))
    }
    ann_df[[additional_annotation_col]] <- metadata_filtered[[additional_annotation_col]]
  }
  
  ##### Build color mapping for annotations #####
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
      if(n_levels < 3){
        extra_pal <- brewer.pal(n = 3, name = "Set1")[seq_len(n_levels)]
      } else {
        extra_pal <- brewer.pal(n = n_levels, name = "Set1")
      }
      ann_colors[[additional_annotation_col]] <- setNames(extra_pal, extra_levels)
    }
  }
  
  top_annotation <- HeatmapAnnotation(
    df = ann_df,
    col = ann_colors,
    annotation_name_gp = gpar(fontsize = 10),
    annotation_legend_param = list(
      Diagnosis = list(title = "Diagnosis", at = levels(ann_df$Diagnosis), labels = levels(ann_df$Diagnosis))
    )
  )
  
  ##### Set up the heatmap color function #####
  col_fun <- colorRamp2(
    c(min(count_matrix, na.rm = TRUE), 0, max(count_matrix, na.rm = TRUE)),
    c("blue", "white", "red")
  )
  
  ##### Clustering and Heatmap Creation with Split Clusters #####
  ht <- NULL
  clusters <- list()
  
  if (cluster_by == "rows") {
    row_dist <- dist(count_matrix, method = "euclidean")
    row_clust <- hclust(row_dist, method = "ward.D2")
    row_order <- seriation::get_order(seriation::seriate(row_dist, method = "OLO"))
    clusters_val <- cutree(row_clust, k = n_row_clusters)
    
    ht <- Heatmap(
      count_matrix,
      name = "log2FoldChange",
      col = col_fun,
      cluster_rows = row_clust,     # supply hclust object
      row_order = row_order,        # custom order
      row_split = n_row_clusters,         # pass numeric value to split rows
      cluster_columns = FALSE,
      show_column_names = TRUE,
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 6),
      top_annotation = top_annotation,
      row_dend_side = "right",
      row_dend_width = unit(3, "cm"),
      use_raster = TRUE
    )
    clusters$row_clusters <- clusters_val
    
  } else if (cluster_by == "columns") {
    col_dist <- dist(t(count_matrix), method = "euclidean")
    col_clust <- hclust(col_dist, method = "ward.D2")
    col_order <- seriation::get_order(seriation::seriate(col_dist, method = "OLO"))
    clusters_val <- cutree(col_clust, k = n_column_clusters)
    
    ht <- Heatmap(
      count_matrix,
      name = "log2FoldChange",
      col = col_fun,
      cluster_rows = FALSE,
      cluster_columns = col_clust,  # supply hclust object
      column_order = col_order,     # custom order
      column_split = n_column_clusters,      # pass numeric value to split columns
      show_column_names = TRUE,
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 6),
      top_annotation = top_annotation,
      row_dend_side = "right",
      row_dend_width = unit(3, "cm"),
      use_raster = TRUE
    )
    clusters$column_clusters <- clusters_val
    
  } else if (cluster_by == "both") {
    # Row clustering
    row_dist <- dist(count_matrix, method = "euclidean")
    row_clust <- hclust(row_dist, method = "ward.D2")
    row_order <- seriation::get_order(seriation::seriate(row_dist, method = "OLO"))
    clusters_row <- cutree(row_clust, k = n_row_clusters)
    
    # Column clustering
    col_dist <- dist(t(count_matrix), method = "euclidean")
    col_clust <- hclust(col_dist, method = "ward.D2")
    col_order <- seriation::get_order(seriation::seriate(col_dist, method = "OLO"))
    clusters_col <- cutree(col_clust, k = n_column_clusters)
    
    ht <- Heatmap(
      count_matrix,
      name = "log2FoldChange",
      col = col_fun,
      cluster_rows = row_clust,
      row_order = row_order,
      row_split = n_row_clusters,        # split rows into n_row_clusters
      cluster_columns = col_clust,
      column_order = col_order,
      column_split = n_column_clusters,     # split columns into n_column_clusters
      show_column_names = TRUE,
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 6),
      top_annotation = top_annotation,
      row_dend_side = "right",
      row_dend_width = unit(3, "cm"),
      use_raster = TRUE
    )
    clusters$row_clusters <- clusters_row
    clusters$column_clusters <- clusters_col
    
  } else {
    stop("Invalid cluster_by parameter. Please choose 'rows', 'columns', or 'both'.")
  }
  
  return(list(ht = ht, clusters = clusters))
}


subject_heatmap <- plot_genes_subject_heatmap(l2fc_data,
                                              diagnoses = DIAGNOSES,
                                              visit = "V3",
                                              gene_set = gene_set,
                                              scale_rows = FALSE,
                                              cluster_by = "columns",
                                              n_column_clusters = 3,
                                              additional_annotation_col = NULL)
subject_heatmap$ht


# save the heatmap as a PNG file
png("figs/analysis/subject_l2fc_heatmap/custom_geneset_1.png", width = 4000, height = 3600, res = 300) # adjust width, height, and resolution
draw(subject_heatmap$ht)
dev.off()


# baseline pathway and gene count boxplots ----------------------------------------------------

# define baseline pathways of interest (poi)
baseline_poi = c(
  "cell cycle and transcription (M4.0)",
  "enriched in myeloid cells and monocytes (M81)",
  "immune activation - generic cluster (M37.0)",
  "NK cell surface signature (S1)",
  "T cell activation (I) (M7.1)",
  "enriched in B cells (I) (M47.0)",
  "plasma cells & B cells, immunoglobulins (M156.0)",
  "extracellular matrix, complement (M140)",
  "Resting dendritic cell surface signature (S10)"
)

# create GSVA plots comparing diagnoses
# save a pca biplot with loadings for each pathway and a boxplot by diagnosis
for (pathway in baseline_poi) {
  ### PCA BIPLOT ###
  
  # plot pc1 and pc2 in biplot with loadings to get an idea of the genes driving this pathway
  # get genes associated with that pathway
  pathway_genes <- PATHWAYS[[pathway]]
  
  # subset count matrix to genes of pathway
  # remove those genes defined in the pathway but not present in our count_matrix
  vsd_pathway_subset <- drop_na(V1_vsd[pathway_genes, ])
  
  # run PCA on the module
  # removeVar=0 means keep all genes (no variance filtering)
  p <- PCAtools::pca(
    mat       = as.matrix(vsd_pathway_subset),
    metadata  = V1_metadata,
    removeVar = 0
  )
  
  png(filename = paste0("figs/analysis/baseline_pca/COHORT_V1_",
                        sub(".*\\(([^()]*)\\)\\s*$", "\\1", pathway),
                        "_pca_pc1_pc2.png"), width = 2600, height = 2400, res = 300)
  print(biplot(p,
               colby = "DIAGNOSIS",
               colkey = DIAGNOSIS_COLKEY,
               x = "PC1",
               y = "PC2",
               lab = NULL,
               showLoadings = TRUE,
               legendPosition = "top",
               title = paste0("PCA - ", "‘", pathway, "’"))
  )
  dev.off()
  
  ### GSVA BOXPLOT ###
  
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
  boxp <- ggplot(pathway_df, aes(x = DIAGNOSIS, y = pathway_gsva_score, colour = DIAGNOSIS)) +
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
    stat_compare_means(method = "kruskal.test",
                       label.y = max(pathway_df$pathway_gsva_score) * 1.05,
                       show.legend = FALSE)
  
  ggsave(paste0("figs/analysis/baseline_boxplot/COHORT_V1_",
                sub(".*\\(([^()]*)\\)\\s*$", "\\1", pathway),
                "_gsva_score.png"), plot = boxp, width = 8, height = 6, dpi = 300)
  
}

# create boxplots of genes
genes_of_interest <- c(
  # M81
  "CD14",
  "SNCA",
  "CCR2",
  
  # M47.0
  "PAX5",
  "CD22",
  "CD79A",
  
  # M37.0
  "LTF",
  "DEFA1",
  "CD177",
  "OLFM4",
  "MMP8",
  
  # M4.0
  "E2F2",
  "CEACAM8",
  "PGLYRP1",
  "MKI67",
  "S100A12",
  
  # M7.1
  "CCR7",
  "LEF1",
  "MS4A1",
  "TCF7",
  "IL7R",
  
  "CD8A",
  "CD8B",
  "EOMES",
  
  # M11.1
  "F13A1",
  "FCER1G",
  "SLC7A7",
  
  # M140
  "CRISP3",
  "PRTN3",
  "LAMA2",
  
  # M145.1
  "TPM1",
  "LMOD1",
  "MYLK",
  "MYL9",
  "CALD1",
  
  # NK cells, S1
  "GZMH",
  "NKG7",
  "GNLY",
  "NCAM1", # CD56
  "KIR3DL2", # killer immunoglobulin-like receptors (KIR), for cytotoxic activity
  "FCGR3A",
  
  # S10
  "SLAMF8",
  "CLEC10A",
  "EMP1",
  "STAB1",
  "FZD5",
  
  # M49
  "ANK1",
  "KLF1",
  "SPTB",
  "EPB42",
  
  # something erythropoiesis related
  "SPTB",
  "KLF1",
  
  # personal additions to plasma cells, t-cells etc., likely won't be significant
  "IGJ",
  "MZB1",
  "IL7"
)

for (gene in genes_of_interest) {
  # boxplot gene by diagnosis
  boxp <- boxplot_gene_count_with_test(count_matrix = V1_count_matrix,
                                       metadata = V1_metadata,
                                       gene = gene,
                                       group_col = "DIAGNOSIS", 
                                       comparisons = list(c("HC", "SMM"), c("HC", "NSCLC"), c("HC", "MM")),
                                       test_method = "wilcox.test") +
    labs(
      title = gene,
      # title = paste0("COHORT V1 - ‘", gene, "’", " gene counts"),
      x     = "",
      y     = "d0 normalized gene counts"
    ) +
    scale_color_manual(values = DIAGNOSIS_COLKEY)
  ggsave(paste0("figs/analysis/baseline_boxplot/COHORT_V1_", gene, "_colby_DIAGNOSIS.png"), plot = boxp, width = 6, height = 6, dpi = 300)
}


# induced pathway and subject l2fc gene count boxplots -----------------------------------------------------

# define induced pathways of interest (poi)
induced_poi = c(
  "activated dendritic cells (M67)", # interferon related?
  "Activated (LPS) dendritic cell surface signature (S11)", # interferon related?
  "MHC-TLR7-TLR8 cluster (M146)", # interferon & innate related?
  "lysosomal/endosomal proteins (M139)", # interferon related?
  
  "platelet activation (II) (M32.1)",
  "immune activation - generic cluster (M37.0)", # innate related, at V3?
  
  "enriched in T cells (I) (M7.0)",
  
  "cell cycle (I) (M4.1)",
  "proteasome (M226)",
  "Plasma cell surface signature (S3)",
  "plasma cells, immunoglobulins (M156.1)"
)

# create boxplots of genes (5 genes each, for 5 interesting pathways)
genes_of_interest <- c(
  ### V2 innate and inflammation related ###
  # M67 interferon related genes
  "RSAD2",
  "SERPING1",
  # S11 interferon related genes
  "PDCD1LG2",
  "CD274",
  # more interferon related genes
  "IFIT1",
  "ATF3",
  "OAS1",
  "IRF7",
  "MX2",
  "STAT1",
  
  # M146 innate genes and antigen presentation
  "HLA-G",
  "HLA-C",
  "HLA-DRA",
  "HLA-DQA1",
  "HLA-DPA1",
  "TLR7",
  "TLR8",
  
  "TLR4",
  "CLEC4E",
  
  # M139 lysosome, endosome, maybe related to antigen presentation?
  "SORT1",
  "CTSS",
  "CTSB",
  
  ### V3 plasma cell and proliferation related ###
  # S3 plasma cell surface signature
  "KCNN3",
  "CAV1",
  "SDC1",
  # M156.1 plasma cells and immunoglobulins
  "TNFRSF17",
  "TXNDC5",
  "POU2AF1",
  "PNOC",
  "DERL3",
  "CD27",
  
  # M37.0 at V3, neutrophil related?
  "CD177",
  "ARG1",
  "CTSG",
  "PRTN3",
  "DEFA1",
  "MMP8",
  "AZU1",
  "LTF",
  "ELANE",
  "CHIT1",
  "OLR1",
  "BMX",
  "MMP9",
  
  # M4.1 at V3
  "MKI67",
  "KIF20A",
  "KIF20B",
  "E2F8",
  "CDC45",
  "CDC20",
  "RRM2",
  "CDCA3",
  "CENPA",
  
  "SKA3",
  "NUF2",
  "CENPE",
  "CENPK",
  "SGOL2",
  
  # M7.0 at V3 SMM
  "KLRB1",
  "XCL1",
  "XCL2",
  "GZMH",
  "ICOS",
  "CD3D",
  
  # M226 at V3 NSCLC
  "PSMA4",
  "POLR2K",
  "PSMA3"
)


# create GSVA plots comparing diagnoses
# save a pca biplot with loadings for each pathway and a boxplot by diagnosis
for (pathway in induced_poi) {
  ### PCA BIPLOT ###
  
  # plot pc1 and pc2 in biplot with loadings to get an idea of the genes driving this pathway
  # get genes associated with that pathway
  pathway_genes <- PATHWAYS[[pathway]]
  
  for (visit in VISITS[2:4]) {
    # subset count matrix to genes of pathway
    # remove those genes defined in the pathway but not present in our count_matrix
    l2fc_data_imputed_subset <- subset_l2fc_data(l2fc_data_imputed, condition = paste0("VISIT == ", "'", visit, "'"))
    
    common_genes <- intersect(pathway_genes, rownames(l2fc_data_imputed_subset$count_matrix))
    if (length(common_genes) < 2) {
      warning(sprintf("Fewer than 2 genes found for pathway '%s' at visit %s – skipping PCA.", pathway, visit))
      next
    }
    mat_subset <- as.matrix(l2fc_data_imputed_subset$count_matrix[common_genes, , drop = FALSE])
    
    # run PCA on the module
    # removeVar=0 means keep all genes (no variance filtering)
    p <- PCAtools::pca(
      mat       = mat_subset,
      metadata  = l2fc_data_imputed_subset$metadata,
      removeVar = 0
    )
    
    png(filename = paste0("figs/analysis/l2fc_pca/COHORT_", visit, "_",
                          sub(".*\\(([^()]*)\\)\\s*$", "\\1", pathway),
                          "_pca_pc1_pc2.png"), width = 2600, height = 2400, res = 300)
    print(biplot(p,
                 colby = "DIAGNOSIS",
                 colkey = DIAGNOSIS_COLKEY,
                 x = "PC1",
                 y = "PC2",
                 lab = NULL,
                 showLoadings = TRUE,
                 legendPosition = "top",
                 title = paste0("PCA V1 to ", visit, " l2fc - ", "‘", pathway, "’"))
    )
    dev.off()
    
    ### GSVA BOXPLOT ###
    
    # create visit subset
    l2fc_pathway_data_subset <- subset_pathway_data(l2fc_pathway_data, condition = paste0("VISIT == ", "'", visit, "'"))
    
    # extract the pathway's gsva score
    pathway_df <- data.frame(
      CEGAT_ID    = rownames(l2fc_pathway_data_subset$metadata),
      pathway_gsva_score = l2fc_pathway_data_subset$gsva_matrix[pathway, ],
      stringsAsFactors = FALSE
    )
    
    # merge with metadata
    pathway_df <- merge(
      pathway_df,
      data.frame(l2fc_pathway_data_subset$metadata, stringsAsFactors = FALSE),
      by = "CEGAT_ID"
    )
    
    # boxplot of module gsva score by diagnosis
    boxp <- ggplot(pathway_df, aes(x = DIAGNOSIS, y = pathway_gsva_score, colour = DIAGNOSIS)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
      labs(
        title = paste0("Module l2fc GSVA score for ‘", pathway, "’", " at ", visit),
        x     = "Diagnosis",
        y     = "GSVA score"
      ) +
      theme_bw() +
      theme(
        plot.title   = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 45, hjust = 1)
      ) +
      scale_colour_manual(values = DIAGNOSIS_COLKEY) +
      stat_compare_means(method = "kruskal.test",
                         label.y = max(pathway_df$pathway_gsva_score) * 1.05,
                         show.legend = FALSE)
    
    ggsave(paste0("figs/analysis/l2fc_boxplot/COHORT_", visit, "_",
                  sub(".*\\(([^()]*)\\)\\s*$", "\\1", pathway),
                  "_gsva_score.png"), plot = boxp, width = 8, height = 6, dpi = 300)
  }
}

for (visit in VISITS[2:4]) {
  
  # initialize day variable
  day_variable <- "d0"
  
  if (visit == "V2") {
    day_variable <- "d1"
  } else if (visit == "V3") {
    day_variable <- "d7"
  } else if (visit == "V4") {
    day_variable <- "d30"
  }
  
  # subset count matrix to genes of pathway
  # remove those genes defined in the pathway but not present in our count_matrix
  l2fc_data_imputed_subset <- subset_l2fc_data(l2fc_data_imputed, condition = paste0("VISIT == ", "'", visit, "'"))
  
  for (gene in genes_of_interest) {
    # boxplot gene by diagnosis
    boxp <- boxplot_gene_count_with_test(l2fc_data_imputed_subset$count_matrix,
                                         l2fc_data_imputed_subset$metadata,
                                         gene = gene,
                                         group_col = "DIAGNOSIS", 
                                         comparisons = list(c("HC", "SMM"), c("HC", "NSCLC"), c("HC", "MM")),
                                         test_method = "wilcox.test") +
      labs(
        title = gene,
        # title = paste0("COHORT V1 - ‘", gene, "’", " l2fc in gene counts"),
        x     = "",
        y     = paste0(day_variable, " vs d0 log2FoldChange in normalized gene counts")
      ) +
      scale_color_manual(values = DIAGNOSIS_COLKEY)
    ggsave(paste0("figs/analysis/l2fc_boxplot/COHORT_", visit, "_", gene, "_colby_DIAGNOSIS.png"), plot = boxp, width = 8, height = 6, dpi = 300)
  }
}

# manually check loadings for pathway pca to identify interesting genes
l2fc_data_imputed_subset <- subset_l2fc_data(l2fc_data_imputed, condition = paste0("VISIT == ", "'", "V3","'"))

pathway <- "enriched in T cells (I) (M7.0)"

pathway_genes <- PATHWAYS[[pathway]]

common_genes <- intersect(pathway_genes, rownames(l2fc_data_imputed_subset$count_matrix))
if (length(common_genes) < 2) {
  warning(sprintf("Fewer than 2 genes found for pathway '%s' at visit %s – skipping PCA.", pathway, visit))
  next
}
mat_subset <- as.matrix(l2fc_data_imputed_subset$count_matrix[common_genes, , drop = FALSE])

# run PCA on the module
# removeVar=0 means keep all genes (no variance filtering)
p <- PCAtools::pca(
  mat       = mat_subset,
  metadata  = l2fc_data_imputed_subset$metadata,
  removeVar = 0
)

biplot(p,
       colby = "DIAGNOSIS",
       colkey = DIAGNOSIS_COLKEY,
       lab = NULL,
       showLoadings = TRUE)





