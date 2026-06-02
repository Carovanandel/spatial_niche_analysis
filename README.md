# Spatial Transcriptomics Niche Analysis
Scripts used for my MSc thesis titled 'From Tissue Niches to Trajectories: Modeling Inclusion Body Myositis Pathogenesis using Spatial Transcriptomics'. \
Nov 2025 - June 2026

## Conda environments
For installation instructions of the used conda environments see `conda_environments/envs_installation.md`

## Description of notebooks

### Napari environment
- `480_5000_fiber_integration_part1.ipynb`: Per-lobe landmark-based alignment of fibers in the 5000 slide to the 480 slide.
### CellCharter environment
- `480_5000_fiber_integration_part2.ipynb`: Matching aligned fiber segmentation shapes between the two consecutive tissue slides based on min shape overlap and max size difference, followed by analysis of overlapping gene expression in matched fibers. Based on this analyis a third fiber matching criteria was added: pearson's r >0.6 for overlapping gene expression.
- `480_5000_fiber_integration_part3.ipynb`: Integration of matched fibers in the  two consecutive tissue slides by appending gene expression from the 480 genes dataset to the 5000 genes dataset and calculating mean counts for overlapping genes. Results in a new integrated fiber dataset anndata file.
- `480_5000_fiber_integration_QC.ipynb`: Calculating QC metrics for the integrated fiber dataset, size normalization of integrated gene counts and plotting PCA and trVAE UMAPs to prepare for pseudotime analysis with Palantir.
- `combine_spatial_shapes.ipynb`: Create combined segmentation shapes layer/geodataframe from the separate shape layers in the 480 and 5000 dataset. (Each sample + 'damaged' myofibers and all nuclei were in separate shape layers in the spatialdata object)
- `immunecells_neighborhood_analysis`: Analyis of fiber gene expression correlating with CD8/CD4/B-cell counts, using the integrated fiber dataset.
- `immunecells_neighborhood_expanded_shapes.ipynb`: Calculating number of immune cells in neighborhood of fibers using the segmentation masks +  initial analysis of genes correlating with CD8 counts.
- `scVI_trVAE_umaps.ipynb`: Comparing cellcharter embedding UMAPs with scVI vs trVAE.
- `stability_test_GMM_leiden_fibnuc.py`: Calculating clustering stability for fiber+nuclei spatial clustering for 480 and 5000 dataset, using both GMM clustering and Leiden clustering.
- `stability_test_leiden_fibnuc_fibers.py`: Calculating clustering stability for fiber only  spatial clustering for 480 and 5000 dataset, using Leiden clustering.
- `spatial_clustering_plot_stability.ipynb`: Plot clustering stability (ARI) of all tested spatial clusterings (using csv's from stability_test scripts).
- `spatial_clustering_final.ipynb`: TO DO
- `spatial_clustering_integrated_fibers_analysis.ipynb`: TO DO
- `spatial_clustering_plotting.ipynb`: Plotting spatial (fiber) niches and immune cell counts over the segmentation shapes - Plots used in final thesis.

### Palantir environment
- `palantir_pseudotime.ipynb`: TO DO

### scArches environment
- `trVAE_final.ipynb`: TO DO
- `trvae_integrated_fibers.ipynb`: TO DO

### SPArrOW environment
- `immunecells_neighborhood_masks_visualise.ipynb`: TO DO
