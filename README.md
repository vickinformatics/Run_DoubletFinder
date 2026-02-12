# Run_DoubletFinder
- If the user chooses not to define the number of PCs for PCA or Harmony (i.e., setting `PCs_pca` or `PCs_harmony` as NULL), they must source the custom `FindMinimumPCs()` function before running `Run_Seurat()`, as this will automatically determine the minimum number of PCs. Please refer to the documentation for [FindMinimumPCs](https://github.com/vickinformatics/FindMinimumPCs).

## Description
The `Run_DoubletFinder()` function runs the full DoubletFinder pipeline on one or more Seurat objects present in the R environment. It handles preprocessing (normalization, variable feature selection, scaling, PCA), automatic PC selection, pK optimization via parameter sweep, homotypic doublet proportion estimation, and doublet classification using a two-pass approach. Results are stored in a list called `DoubletFinder_processed` in the global environment.

## Arguments
- `sample_names` A character vector of one or more Seurat object names present in the R environment (e.g., "sample1" or c("sample1", "sample2", "sample3")).
- `doublet_rate` Either a numeric value between 0 and 1 specifying the expected proportion of doublets (e.g., 0.075 for 7.5%), or "dynamic" to automatically estimate the rate per sample using the 10x Genomics standard formula: 0.008 * (ncol(seurat) / 1000). Defaults to "dynamic". Override with a numeric value when the doublet rate is known.
- `PCs` The number of principal components to use. If NULL, PCs will be calculated automatically using the `FindMinimumPCs()` function (default is NULL).
- `nfeatures` The number of variable features to select during preprocessing (default is 2000).
- `cluster_resolution` Resolution passed to `FindClusters()` for generating the cluster annotations used in homotypic proportion estimation (default is 0.1).
- `sct` Logical. Whether the Seurat object was processed using SCTransform. Passed to `paramSweep()` and `doubletFinder()` (default is FALSE).
- `return_singlets_only` Logical. If TRUE, each object in `DoubletFinder_processed` contains only singlet cells. If FALSE, all cells are retained with a `DoubletFinder` metadata column containing "Singlet" or "Doublet" classifications (default is TRUE).
- `variables_to_regress` A character vector of metadata variables to regress out during scaling (default is NULL).
- `normalize.args` A named list of additional arguments to pass to NormalizeData() (default is NULL).
- `variable.features.args` A named list of additional arguments to pass to FindVariableFeatures() (default is NULL).
- `scale.args` A named list of additional arguments to pass to ScaleData() (default is NULL).
- `pca.args` A named list of additional arguments to pass to RunPCA() (default is NULL).
- `umap.args` A named list of additional arguments to pass to RunUMAP() (default is NULL).
- `neighbors.args` A named list of additional arguments to pass to FindNeighbors() (default is NULL).
- `clusters.args` A named list of additional arguments to pass to FindClusters() (default is NULL).
- `doubletfinder.args` A named list of additional arguments to pass to doubletFinder() (default is NULL).

## References
- McGinnis CS, Murrow LM, Gartner ZJ. DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors. Cell Systems (2019).
