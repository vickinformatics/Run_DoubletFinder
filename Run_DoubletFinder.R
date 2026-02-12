#' Run_DoubletFinder
#'
#' The `Run_DoubletFinder()` function runs the full DoubletFinder pipeline on one or more Seurat objects present in the R environment.
#' It handles preprocessing (normalization, variable feature selection, scaling, PCA), automatic PC selection, pK optimization via parameter sweep, homotypic doublet proportion estimation, and doublet classification using a two-pass approach.
#' Results are stored in a list called `DoubletFinder_processed` in the global environment.
#'
#' @param sample_names A character vector of one or more Seurat object names present in the R environment (e.g., "sample1" or c("sample1", "sample2", "sample3")).
#' @param doublet_rate Either a numeric value between 0 and 1 specifying the expected proportion of doublets (e.g., 0.075 for 7.5%), or "dynamic" to automatically estimate the rate per sample using the
#'  10x Genomics standard formula: 0.008 * (ncol(seurat) / 1000).
#'  Defaults to "dynamic". Override with a numeric value when the doublet rate is known.
#' @param PCs The number of principal components to use. If NULL, PCs will be calculated automatically using the `FindMinimumPCs()` function (default is NULL).
#' @param nfeatures The number of variable features to select during preprocessing (default is 2000).
#' @param cluster_resolution Resolution passed to `FindClusters()` for generating the cluster annotations used in homotypic proportion estimation (default is 0.1).
#' @param sct Logical. Whether the Seurat object was processed using SCTransform. Passed to `paramSweep()` and `doubletFinder()` (default is FALSE).
#' @param return_singlets_only Logical.
#'  If TRUE, each object in `DoubletFinder_processed` contains only singlet cells.
#'  If FALSE, all cells are retained with a `DoubletFinder` metadata column containing "Singlet" or "Doublet" classifications (default is TRUE).
#' @param variables_to_regress A character vector of metadata variables to regress out during scaling (default is NULL).
#' @param normalize.args A named list of additional arguments to pass to NormalizeData() (default is NULL).
#' @param variable.features.args A named list of additional arguments to pass to FindVariableFeatures() (default is NULL).
#' @param scale.args A named list of additional arguments to pass to ScaleData() (default is NULL).
#' @param pca.args A named list of additional arguments to pass to RunPCA() (default is NULL).
#' @param umap.args A named list of additional arguments to pass to RunUMAP() (default is NULL).
#' @param neighbors.args A named list of additional arguments to pass to FindNeighbors() (default is NULL).
#' @param clusters.args A named list of additional arguments to pass to FindClusters() (default is NULL).
#' @param doubletfinder.args A named list of additional arguments to pass to doubletFinder() (default is NULL).
#'
#' @return Invisibly returns the `DoubletFinder_processed` list, which is also assigned to the global environment.
#'  Each element is named after its corresponding input sample and contains the processed Seurat object.
#'  If `return_singlets_only = TRUE`, only singlet cells are retained.
#'  If `return_singlets_only = FALSE`, all cells are retained with a `DoubletFinder` metadata column added.
#'
#' @note If `PCs` is NULL, the `FindMinimumPCs()` function must be sourced before running `Run_DoubletFinder()`, as it is used to automatically determine the minimum number of informative PCs.
#'
#' @references McGinnis CS, Murrow LM, Gartner ZJ (2019). "DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors." *Cell Systems*, 8(4), 329-337. doi:10.1016/j.cels.2019.03.003.
#'
#' @author Vicki Do
#' @lastUpdated 2026-02-12
#'
#' @examples
#' # Run on a single sample with a known doublet rate
#' Run_DoubletFinder(sample_names = "sample1", doublet_rate = 0.075)
#'
#' # Run on multiple samples with a known doublet rate
#' Run_DoubletFinder(sample_names = c("sample1", "sample2", "sample3"),
#'                   doublet_rate = 0.075)
#'
#' # Run on multiple samples with an unknown doublet rate (dynamic estimation)
#' Run_DoubletFinder(sample_names = c("sample1", "sample2", "sample3"),
#'                   doublet_rate = "dynamic")
#'
#' # Return full annotated objects instead of just singlets (useful for QC/UMAP inspection)
#' Run_DoubletFinder(sample_names = c("sample1", "sample2"), doublet_rate = 0.075,
#'                   return_singlets_only = FALSE)
#'
#' # Access results after running
#' DoubletFinder_processed[["sample1"]]

Run_DoubletFinder <- function(sample_names,
                              doublet_rate = "dynamic",
                              PCs = NULL,
                              nfeatures = 2000,
                              cluster_resolution = 0.1,
                              sct = FALSE,
                              return_singlets_only = TRUE,
                              variables_to_regress = NULL,
                              normalize.args = NULL,
                              variable.features.args = NULL,
                              scale.args = NULL,
                              pca.args = NULL,
                              umap.args = NULL,
                              neighbors.args = NULL,
                              clusters.args = NULL,
                              doubletfinder.args = NULL) {
  
  # Validate doublet_rate input
  if (!identical(doublet_rate, "dynamic") && (!is.numeric(doublet_rate) || doublet_rate <= 0 || doublet_rate >= 1)) {
    stop("doublet_rate must be a numeric value between 0 and 1 (e.g., 0.075) or 'dynamic'.")
  }
  
  # Validate sample_names input
  if (!is.character(sample_names) || length(sample_names) == 0) {
    stop("sample_names must be a non-empty character vector of Seurat object names.")
  }
  
  # Initialize output list
  DoubletFinder_processed <- list()
  
  for (sample_name in sample_names) {
    cat("\n------ Processing:", sample_name, "------\n")
    
    # Retrieve object from global environment
    if (!exists(sample_name, envir = .GlobalEnv)) {
      warning("Object '", sample_name, "' not found in the R environment. Skipping.")
      next
    }
    seurat <- get(sample_name, envir = .GlobalEnv)
    
    # Preprocessing
    
    # Normalization
    if (is.null(normalize.args)) {
      seurat <- NormalizeData(seurat)
    } else {
      norm_args <- modifyList(list(object = seurat), normalize.args)
      seurat <- do.call(NormalizeData, norm_args)
    }
    
    # Find variable features
    if (is.null(variable.features.args)) {
      seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = nfeatures)
    } else {
      fvf_args <- modifyList(list(object = seurat, selection.method = "vst", nfeatures = nfeatures),
                             variable.features.args)
      seurat <- do.call(FindVariableFeatures, fvf_args)
    }
    
    # Scale data
    if (is.null(scale.args)) {
      seurat <- ScaleData(seurat, vars.to.regress = variables_to_regress)
    } else {
      scale_args <- modifyList(list(object = seurat, vars.to.regress = variables_to_regress), scale.args)
      seurat <- do.call(ScaleData, scale_args)
    }
    
    # Run PCA
    if (is.null(pca.args)) {
      seurat <- RunPCA(seurat)
    } else {
      pca_args <- modifyList(list(object = seurat), pca.args)
      seurat <- do.call(RunPCA, pca_args)
    }
    
    # Use provided PCs or calculate using FindMinimumPCs
    pcs_to_use <- PCs
    if (is.null(pcs_to_use)) {
      pcs_to_use <- FindMinimumPCs(seurat, reduction_type = "pca")
    }
    cat("Using", pcs_to_use, "PCs\n")
    
    # Clustering (required for homotypic proportion estimate)
    
    # UMAP
    if (is.null(umap.args)) {
      seurat <- RunUMAP(seurat, dims = 1:pcs_to_use)
    } else {
      umap_args <- modifyList(list(object = seurat, dims = 1:pcs_to_use), umap.args)
      seurat <- do.call(RunUMAP, umap_args)
    }
    
    # Find neighbors
    if (is.null(neighbors.args)) {
      seurat <- FindNeighbors(seurat, dims = 1:pcs_to_use)
    } else {
      fn_args <- modifyList(list(object = seurat, dims = 1:pcs_to_use), neighbors.args)
      seurat <- do.call(FindNeighbors, fn_args)
    }
    
    # Find clusters
    if (is.null(clusters.args)) {
      seurat <- FindClusters(seurat, resolution = cluster_resolution)
    } else {
      fc_args <- modifyList(list(object = seurat, resolution = cluster_resolution), clusters.args)
      seurat <- do.call(FindClusters, fc_args)
    }
    
    # pK Optimization
    
    sweep.list <- paramSweep(seurat, PCs = 1:pcs_to_use, sct = sct)
    sweep.stats <- summarizeSweep(sweep.list)
    bcmvn <- find.pK(sweep.stats)
    optimal.pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    cat("Optimal pK:", optimal.pk, "\n")
    
    # Expected doublet estimation
    
    # Resolve doublet rate (dynamic = 0.8% per 1,000 cells, per 10x Genomics standard)
    resolved_rate <- if (identical(doublet_rate, "dynamic")) {
      0.008 * (ncol(seurat) / 1000)
    } else {
      doublet_rate
    }
    cat("Doublet rate used:", resolved_rate * 100, "%\n")
    
    # Homotypic doublet proportion
    homotypic.prop <- modelHomotypic(seurat$seurat_clusters)
    
    # nExp based on resolved doublet_rate
    nExp.poi <- round(resolved_rate * ncol(seurat))
    nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
    cat("Expected doublets (adjusted):", nExp.poi.adj, "\n")
    
    # Two-pass DoubletFinder
    
    # Pass 1: generate pANN scores using unadjusted nExp
    if (is.null(doubletfinder.args)) {
      seurat <- doubletFinder(seurat, PCs = 1:pcs_to_use, pK = optimal.pk,
                              nExp = nExp.poi, sct = sct)
    } else {
      df_args <- modifyList(list(seu = seurat, PCs = 1:pcs_to_use, pK = optimal.pk,
                                 nExp = nExp.poi, sct = sct), doubletfinder.args)
      seurat <- do.call(doubletFinder, df_args)
    }
    
    # Retrieve pANN column from Pass 1
    pANN_col <- grep("^pANN", colnames(seurat@meta.data), value = TRUE)
    if (length(pANN_col) == 0) stop("pANN column not found after Pass 1 of DoubletFinder.")
    pANN_col <- pANN_col[length(pANN_col)]  # take most recent if multiple exist
    
    # Pass 2: reclassify using homotypic-adjusted nExp, reusing pANN scores
    if (is.null(doubletfinder.args)) {
      seurat <- doubletFinder(seurat, PCs = 1:pcs_to_use, pK = optimal.pk,
                              nExp = nExp.poi.adj, reuse.pANN = pANN_col, sct = sct)
    } else {
      df_args2 <- modifyList(list(seu = seurat, PCs = 1:pcs_to_use, pK = optimal.pk,
                                  nExp = nExp.poi.adj, reuse.pANN = pANN_col, sct = sct),
                             doubletfinder.args)
      seurat <- do.call(doubletFinder, df_args2)
    }
    
    # Rename DoubletFinder classification column
    
    doublet_col <- grep("^DF.classifications", colnames(seurat@meta.data), value = TRUE)
    if (length(doublet_col) == 0) stop("No DF.classifications column found. DoubletFinder may have failed.")
    doublet_col <- doublet_col[length(doublet_col)]  # take most recent if multiple exist
    seurat@meta.data$DoubletFinder <- seurat@meta.data[[doublet_col]]
    
    # Report results
    n_singlets <- sum(seurat$DoubletFinder == "Singlet", na.rm = TRUE)
    n_doublets <- sum(seurat$DoubletFinder == "Doublet", na.rm = TRUE)
    cat("Singlets:", n_singlets, "| Doublets:", n_doublets, "\n")
    
    # Store results
    
    if (return_singlets_only) {
      DoubletFinder_processed[[sample_name]] <- subset(seurat, DoubletFinder == "Singlet")
    } else {
      DoubletFinder_processed[[sample_name]] <- seurat
    }
  }
  
  # Assign to global environment and return invisibly
  assign("DoubletFinder_processed", DoubletFinder_processed, envir = .GlobalEnv)
  message("\nDone. Results stored in: DoubletFinder_processed")
  invisible(DoubletFinder_processed)
}