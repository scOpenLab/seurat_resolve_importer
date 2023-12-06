# RESOLVE_analysis
Tentative set of tools and scripts for analysing spatial transcriptomic data with the resolve platform.

## Seurat Loading
+ [Functions for loading resolve data](https://codebase.helmholtz.cloud/resolve_tools/resolve-analysis/-/blob/master/load_function.R)
+ [Example script](https://codebase.helmholtz.cloud/resolve_tools/resolve-analysis/-/blob/master/resolveseurat.R)

The script is based on this [Seurat vignette](https://satijalab.org/seurat/articles/spatial_vignette_2.html#human-lymph-node-akoya-codex-system).

The script currently requires a development version of Seurat from the [feat/imaging branch](https://github.com/satijalab/seurat/tree/feat/imaging)

Assumes the input cell segmentation was generated with:
https://codebase.helmholtz.cloud/resolve_tools/resolve-processing

However, it should be adaptable to use segmentation output from other tools.

The `ReadResolve` function expects a folder with the following files:
- `*-cell_data.csv`:
 csv file with this header:
  ```
  cell,area,centroid.y,centroid.x,label,GENE1,GENE2,...
  ```
    - `cell` = table index, not used.
    - `centroid.y,centroid.x` = centroid coordinates (µm or pixel, see the `use.micron` argument)
    - `area` = area in pixel^2
    -  `GENE1,GENE2,...`= transcript counts per cell.
- `*-filtered_transcripts.txt` (either the raw output from resolve or the deduplicated output from [MindaGap](https://github.com/ViriatoII/MindaGap)):
    csv file with no header and these column:
    - x: pixels 
    - y: pixels
    - z: not used, but required
    - gene name
    - quality: not used, optional
 - `*-roi.zip` (optional, if not provided the centroids are used), zip format used by FiJi.
 - `*_mask.tiff`: segmentation mask from cellpose or similar tool (tiff file, used only to get the total width and height, so any image would work)
 - `*-gridfilled.tiff`: optional image to be added to the Seurat object for visualization. Generally we use the output from [MindaGap](https://github.com/ViriatoII/MindaGap), but any image would work. 
