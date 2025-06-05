# ==========================
# aux.R
# ==========================

# Load all required libraries
library(Matrix)
library(DropletUtils)
library(SingleCellExperiment)
library(rtracklayer)
library(scuttle)
library(scater)
library(scran)
library(scDblFinder)
library(edgeR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reticulate)
library(sceasy)
library(scales)
library(celda)
library(Seurat)


# 1. Load and annotate 10x data with GTF
load_data <- function(dataset_path, gtf_file) {
  gtf <- rtracklayer::readGFF(gtf_file)
  genes <- gtf[gtf$type == "gene", ]
  gene_names <- setNames(genes$gene_symbol, genes$gene_id)
  sce <- DropletUtils::read10xCounts(dataset_path, col.names = TRUE)
  valid <- intersect(rowData(sce)$ID, genes$gene_id)
  sce <- sce[match(valid, rowData(sce)$ID), ]
  rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID,
                                        gene_names[rowData(sce)$ID])
  rowData(sce)$Symbol <- gene_names[rowData(sce)$ID]
  rm(gtf, genes); gc()
  return(sce)
}

# 2. emptyDrops filtering + MA-based ambient gene removal
empty_droplet <- function(sce_raw, params, negative_features, output_dir, dataset_name) {
  sce <- sce_raw[, colSums(counts(sce_raw)) > 0]
  sce <- sce[rowSums(counts(sce)) > 0, ]

  # First pass: identify droplets
  ed1 <- emptyDrops(counts(sce),
                    lower = params$lower_bound,
                    retain = params$retain,
                    test.ambient = TRUE)
    
    
    plot_data <- data.frame(
    Total = ed1$Total,
    FDR = ed1$FDR
  )

    plot_data <- plot_data %>%
    mutate(
        CellType = case_when(
        !is.na(FDR) & FDR <= params$significance_threshold ~ "Non-empty",
        TRUE ~ "Empty"  
        )
    ) %>%
    filter(!is.na(Total) & Total > 0)  

    counts_table <- plot_data %>%
    group_by(CellType) %>%
    summarise(
        n = n(),
        max_y = max(Total, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(label_y = max_y * 1.2)  
    violin_ed <- ggplot(plot_data, aes(x = CellType, y = Total, fill = CellType)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.6, adjust = 1.2) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
    geom_text(data = counts_table, aes(x = CellType, y = label_y, label = paste0("n = ", n)), 
                inherit.aes = FALSE, vjust = 0) +
    scale_y_log10(labels = scales::comma) +
    scale_fill_manual(values = c("grey70", "dodgerblue3")) +
    labs(
        title = "Library Size Distribution After emptyDrops Cell Calling",
        x = "",
        y = "Library Size (log10 scale)"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
    ggsave(file.path(output_dir, paste0(dataset_name, "_ed_violin1.png")), violin_ed)


  if (!is.null(negative_features)) {
    cells <- which(ed1$FDR <= params$significance_threshold)
    empty <- which(is.na(ed1$FDR))
    ma <- edgeR::maPlot(rowSums(counts(sce)[, cells]),
                        rowSums(counts(sce)[, empty]),
                        normalize = TRUE, plot.it = FALSE)
    to_remove <- rownames(sce)[ma$M > 0 & rownames(sce) %in% negative_features]
    df_ma <- data.frame(A = ma$A, M = ma$M,
                        Removed = rownames(sce) %in% to_remove)
    p_ma <- ggplot(df_ma, aes(A, M, color = Removed)) +
            geom_point(alpha = 0.6) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            scale_color_manual(values = c("grey50","red")) +
            labs(title = paste(dataset_name, "emptyDrops MA plot"),
                 x = "Average log-counts (A)",
                 y = "Log-fold change (M)") +
            theme_minimal()
    ggsave(file.path(output_dir, paste0(dataset_name, "_MAplot.png")), p_ma)
    sce <- sce[!rownames(sce) %in% to_remove, ]
  }

  ed2 <- emptyDrops(counts(sce),
                    lower = params$lower_bound,
                    retain = Inf,
                    test.ambient = TRUE)

  plot_data <- data.frame(
    Total = ed2$Total,
    FDR = ed2$FDR
  )

plot_data <- plot_data %>%
    mutate(
        CellType = case_when(
        !is.na(FDR) & FDR <= params$significance_threshold ~ "Non-empty",
        TRUE ~ "Empty"  
        )
    ) %>%
    filter(!is.na(Total) & Total > 0)  

    counts_table <- plot_data %>%
    group_by(CellType) %>%
    summarise(
        n = n(),
        max_y = max(Total, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(label_y = max_y * 1.2)  
    violin_ed_post_ma <- ggplot(plot_data, aes(x = CellType, y = Total, fill = CellType)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.6, adjust = 1.2) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
    geom_text(data = counts_table, aes(x = CellType, y = label_y, label = paste0("n = ", n)), 
                inherit.aes = FALSE, vjust = 0) +
    scale_y_log10(labels = scales::comma) +
    scale_fill_manual(values = c("grey70", "dodgerblue3")) +
    labs(
        title = "Library Size Distribution After emptyDrops Cell Calling",
        x = "",
        y = "Library Size (log10 scale)"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
    ggsave(file.path(output_dir, paste0(dataset_name, "_ed_violin2.png")), violin_ed_post_ma)
  

  keep <- which(ed2$FDR <= params$significance_threshold)
  sce <- sce[, keep]
  return(sce)
}

QC <- function(sce,
                assay = "counts",
                params = list(nmads = 3),
                output_dir = ".",
                dataset_name = "sample",
                palette = c(`FALSE` = "#4E79A7", `TRUE` = "#E15759")) {
  
  
  sce <- sce[, colSums(assay(sce, assay)) > 0]
  
  is_mito <- grepl("^MT:", rowData(sce)$Symbol, ignore.case = TRUE)
  
  sce <- scuttle::logNormCounts(sce,
                                assay.type   = assay,
                                exprs_values = assay)
  
  
  stats <- scuttle::perCellQCMetrics(sce,
                                     assay.type = "logcounts",
                                     subsets    = list(Mito = is_mito))
  colData(sce) <- cbind(colData(sce), stats)
  
  thr <- list(
    Detected = isOutlier(stats$detected,              type = "lower",  nmads = params$nmads),
    UMIs     = isOutlier(stats$sum,                   type = "lower",  nmads = params$nmads),
    Mito     = isOutlier(stats$subsets_Mito_percent,  type = "higher", nmads = params$nmads)
  )
  sce$discard <- Reduce(`|`, thr)
  

  df <- tibble::tibble(
    Detected = stats$detected,
    UMIs     = stats$sum,
    Mito     = stats$subsets_Mito_percent,
    discard  = sce$discard
  )
  n_kept  <- sum(!df$discard)
  n_disc  <- sum(df$discard)
  legend_lbls <- c(
    `FALSE` = paste0("Kept (n=", n_kept, ")"),
    `TRUE`  = paste0("Discarded (n=", n_disc, ")")
  )
  
  
  
  df_plot <- df
  GEOM_POINT <- if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::geom_point_rast
  } else {
    message("Package 'ggrastr' not installed – using geom_point()")
    ggplot2::geom_point
  }
  n_kept      <- sum(!df$discard)
  n_discarded <- sum(df$discard)
  legend_lbls <- c(
    `FALSE` = paste0("Kept (n=", n_kept, ")"),
    `TRUE`  = paste0("Discarded (n=", n_discarded, ")")
  )
  
  viol_plot <- function(data, y, ylab, log = FALSE, thr_line = NULL, show_legend = FALSE) {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = "", y = .data[[y]], colour = discard)) +
      ggplot2::geom_violin(fill = "grey92", colour = "grey60", width = 0.9, alpha = 0.8) +
      GEOM_POINT(position = ggplot2::position_jitter(width = 0.15), size = 0.28, alpha = 0.6) +
      ggplot2::scale_colour_manual(
        values  = palette,
        labels  = legend_lbls,
        name    = "Cell status",
        guide   = if (show_legend) "legend" else "none"
      ) +
      ggplot2::labs(x = NULL, y = ylab)
    
    if (log) p <- p + ggplot2::scale_y_log10(labels = scales::comma)
    if (!is.null(thr_line)) p <- p + ggplot2::geom_hline(yintercept = thr_line,
                                                         linetype = "dashed", colour = "red")
    p + ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }
  
  p_detected <- viol_plot(df_plot, "Detected",
                          "Detected genes",
                          log = TRUE)
  
  p_umis <- viol_plot(df_plot, "UMIs",
                      "Total UMIs",
                      log = TRUE,
                      show_legend = TRUE,)
  
  p_mito <- viol_plot(df_plot, "Mito",
                      "Mitochondrial %",
                      show_legend = FALSE)
  
  
  p_final <- (p_detected | p_mito | p_umis) +
    patchwork::plot_annotation(title = paste(dataset_name, "QC metrics"),
                               theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5)))
  
  ggplot2::ggsave(file.path(output_dir, paste0(dataset_name, "_QC.png")),
                  plot   = p_final,
                  width  = 11, height = 4, dpi = 320)
  
  sce_filtered <- sce[, !sce$discard]
  saveRDS(sce_filtered, file.path(output_dir, paste0(dataset_name, "_filtered.rds")))
  invisible(sce_filtered)
}




# 4. Decontamination (optional)
decontamination <- function(sce_filtered, sce_raw, params, output_dir, dataset_name) {
  sce_bg <- sce_filtered[, colSums(counts(sce_filtered)) > 0]
  valid_genes <- intersect(rownames(sce_bg), rownames(sce_raw))
  sce_bg <- sce_bg[valid_genes, ]
  sce_decont <- decontX(sce_bg, background = sce_raw[valid_genes, ])
  assay(sce_decont, "decontXcounts") <- as(round(assay(sce_decont, "decontXcounts")), "dgCMatrix")
  saveRDS(sce_decont, file.path(output_dir, paste0(dataset_name, "_decontam.rds")))
  return(sce_decont)
}

# 5. Doublet detection & filtering
gene_n_doublet <- function(sce, params, output_dir, dataset_name) {
  assay_name <- if ("decontXcounts" %in% assayNames(sce)) "decontXcounts" else "counts"
  
  expr_mat   <- assay(sce, assay_name)
  keep_genes <- Matrix::rowSums(expr_mat > 0) > 0
  
  sce <- sce[keep_genes, , drop = FALSE]
  
  sce <- scDblFinder(sce, propRandom = params$propRandom, verbose = TRUE)

  df <- as.data.frame(colData(sce))
  p1 <- ggplot(df, aes(scDblFinder.score, fill=scDblFinder.class)) +
        geom_histogram(position="identity", alpha=0.6) +
        theme_minimal()
  ggsave(file.path(output_dir, paste0(dataset_name, "_dblScores.png")), p1)

  saveRDS(sce, file.path(output_dir, paste0(dataset_name, "_post_doublet.rds")))
  return(sce)
}

# 6. VST normalization & dimension reduction (older version included vst matrix in SCE)
vst_norm_dim_red <- function(sce, params, output_dir, dataset_name) {

  assay_name <- if ("decontXcounts" %in% assayNames(sce)) "decontXcounts" else "counts"

  meta <- colData(sce)
  cell_attr <- data.frame(  
    MitoPct  = meta$subsets_Mito_percent, 
    row.names = colnames(sce),
    check.names = FALSE
  )

 
  vst_out <- sctransform::vst(
    assay(sce, assay_name),
    cell_attr    = cell_attr,
    latent_var   = c("log_umi", "MitoPct"),
    return_cell_attr = TRUE,
    return_corrected_umi = FALSE,
    return_gene_attr     = TRUE,
    vst.flavor    = "v3",
    verbosity     = 2
  )

  # # Subset the SCE to match the vst output dimensions
  # sce <- sce[rownames(vst_out$umi_corrected), colnames(vst_out$umi_corrected), drop = FALSE]

  # # Insert the corrected UMI assay and gene attributes
  # assay(sce, "vst", withDimnames = FALSE) <- vst_out$umi_corrected
  # rowData(sce) <- cbind(rowData(sce), vst_out$gene_attr)

  sce <- runPCA(sce, exprs_values = "logcounts", name = "PCA")
  sce <- runUMAP(sce, dimred = "PCA", name = "UMAP")

  saveRDS(sce, file.path(output_dir, paste0(dataset_name, "_vstred.rds")))
  return(sce)
}

marker_plotting <- function(sce, markers, output_dir, dataset_name) {
  for (m in intersect(markers, rownames(sce))) {
    p <- plotReducedDim(sce, "UMAP", colour_by=m) +
         labs(title=paste(dataset_name, m)) +
         theme_minimal()
    ggsave(file.path(output_dir, paste0(dataset_name, "_marker_", m, ".png")), p)
  }
}

export_h5ad <- function(sce, output_dir, dataset_name) {
  use_virtualenv("r-sceasy-env", required=TRUE)
  outF <- file.path(output_dir, paste0(dataset_name, "_sce.h5ad"))
  sceasy::convertFormat(sce, "sce", "anndata", outF)
  message("Saved h5ad: ", outF)
}