# RESOLVE_tools
Tentative set of tools and scripts for analysing spatial transcriptomic data with the resolve platform.

Assumes the input cell segmentation was generated with:
https://codebase.helmholtz.cloud/resolve_tools/resolve-processing

However, it should be adaptable to use segmentation output from other tools.


## Seurat 
+ [Functions for loading resolve data](https://codebase.helmholtz.cloud/resolve_tools/resolve-analysis/-/blob/master/load_function.R)
+ [Example script](https://codebase.helmholtz.cloud/resolve_tools/resolve-analysis/-/blob/master/resolveseurat.R)

The script is based on this [Seurat vignette](https://satijalab.org/seurat/articles/spatial_vignette_2.html#human-lymph-node-akoya-codex-system).

The script currently requires a development version of Seurat from the [feat/imaging branch](https://github.com/satijalab/seurat/tree/feat/imaging)


