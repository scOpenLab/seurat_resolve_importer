# RESOLVE_tools
Tentative set of tools and scripts for analysing spatial transcriptomic data with the resolve platform

## Contents:
+ (1) [Seurat](##Seurat)
+ (2) [Giotto](##Giotto)

## 1) Seurat <a name="#Seurat"></a>
[Seurat example script](https://github.com/MicheleBortol/RESOLVE_tools/blob/main/bin/resolveseurat.R)

The script is based on this [Seurat vignette](https://satijalab.org/seurat/articles/spatial_vignette_2.html#human-lymph-node-akoya-codex-system).

The script currently requires a development version of Seurat from the [feat/imaging branch](https://github.com/satijalab/seurat/tree/feat/imaging)

To do:
+ Support for cells coming for more than one image in the same Seurat object.
+ Add image to Seurat Object? (Probably not worth it if we use Giotto)

## 2) Giotto <a name="#Giotto"></a>

[Giotto from Seurat object example script](https://github.com/MicheleBortol/RESOLVE_tools/blob/main/bin/resolvegiottoseurat.R)

The script is based on the [seuratToGiotto function](https://github.com/RubD/Giotto/blob/suite/R/interoperability.R).

The script currently requires a development version of Giotto from the [suite branch](https://github.com/RubD/Giotto/tree/suite)

To do:
+ Support for cells coming for more than one image in the same objet?
+ Add direct creation of a Giotto object without passing from Seurat?

Note:
+ The coordinate systems for cells are different between Seurat and Giotto (in Seurat X and Y are inverted). 
+ The transcript coordinates are also imported by `seuratToGiotto` but never used by the package.
