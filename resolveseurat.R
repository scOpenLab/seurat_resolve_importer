rm(list = ls())
library(future)
plan("multisession", workers = 10)
library(data.table)
library(Seurat)
library(tiff)
library(RImageJROI)
library(docstring) # Just for the docstrings

# Source the Load functions
source("load_function.R")

# https://satijalab.org/seurat/articles/spatial_vignette_2.html#human-lymph-node-akoya-codex-system

### Path to files to read:
# OUTPUT_PATH: Path to the output of the Resolve Processing pipeline

OUTPUT_PATH <- ""

# SAMPLE NAME TO TEST
# "_" not accepted
SAMPLE_NAME <- ""

resolve.obj <- LoadResolve(OUTPUT_PATH, sample.name = SAMPLE_NAME)

#### Just as a quick test. See the Seurat docs for how to perform the analysis 
resolve.obj <- subset(resolve.obj, subset = nCount_Resolve > 10) # Skip cells with no transcripts 
resolve.obj <- SCTransform(resolve.obj, assay = "Resolve", verbose = FALSE)
resolve.obj <- FindVariableFeatures(resolve.obj, assay = "Resolve")
resolve.obj <- RunPCA(object = resolve.obj, npcs = 20, verbose = FALSE)
resolve.obj <- RunUMAP(object = resolve.obj, dims = 1:20, verbose = FALSE)
resolve.obj <- FindNeighbors(object = resolve.obj, dims = 1:20, verbose = FALSE)
resolve.obj <- FindClusters(object = resolve.obj, verbose = FALSE,
  resolution = 0.9)

resolve.obj@active.assay <- "Resolve" # To show the real counts in the plots


### Plotting
my_features = c("GENE2", "GENE3", "GENE4", "GENE5")
my_gene = "GENE1"
# No spatial data
nos_dp <- DimPlot(resolve.obj, label = TRUE, label.box = TRUE) + NoLegend() 
nos_fp <- FeaturePlot(resolve.obj, features = my_features, min.cutoff = "q10", max.cutoff = "q90") 

# Plus spatial data (cells as circles on their centroids)
cen_dp <- ImageDimPlot(resolve.obj, cols = "parade", fov = SAMPLE_NAME, boundaries = "centroids", border.color = NA)
cen_fp <- ImageFeaturePlot(resolve.obj, fov = SAMPLE_NAME, features = my_features, boundaries = "centroids", border.color = NA)

# Plus spatial data (cells as segmentation outlines)
seg_dp <- ImageDimPlot(resolve.obj, cols = "parade", fov = SAMPLE_NAME, boundaries = "segmentation", border.color = NA)
seg_fp <- ImageFeaturePlot(resolve.obj, features = my_features,
    min.cutoff = "q11", max.cutoff = "q90", fov = SAMPLE_NAME, boundaries = "segmentation", border.color = NA)

# Raw transcript visualization
# - size = 0 avoids displaying the cells
# - nmols = max number of molecules to display

cen_mp <- ImageFeaturePlot(resolve.obj, fov = SAMPLE_NAME, boundaries = "centroids",
    features = my_gene, size = 0, molecules = my_features, mols.cols = c("#FF0000",
 "#00FF00", "#0000FF", "#FFFF00"), nmols = 9999999999, border.color = NA)

seg_mp <- ImageFeaturePlot(resolve.obj, fov = SAMPLE_NAME, boundaries = "segmentation",
    features = my_gene, molecules = my_features[[4]], mols.cols = c("#0000FF"), nmols = 9999999999,
  min.cutoff = "q10", max.cutoff = "q90", border.color = NA)


####################
# Find Spatially variable features
# The coordinates for the spatial analysis can come from either the centroids or
# the polygons from the cell ROIs. The  GetTissueCoordinates.Centroids or
# GetTissueCoordinates.Segmentation need to be modified. See below.
#
# Both the "moransi" and the "markvariogram" work but the "moransi" method
# requires the installation of the "ape" or the "Rfast2" package

###############
# Change the definitions of these 2 functions from the SeuratObject package
# This is not good practice, but it works as a quick workaround for the
# FindSpatiallyVariableFeatures function which require the returned object
# to have rownames. 
#
# see the fortify.Centroids and the FindSpatiallyVariableFeatures functions
GetTissueCoordinates.Centroids <- function(object, full = TRUE, ...){
  a <- as.data.frame(object@coords)
  a$cell <- object@cells
  rownames(a) <- object@cells
  return(a)
}

assignInNamespace("GetTissueCoordinates.Centroids",
                  GetTissueCoordinates.Centroids, ns = "SeuratObject")

# see the fortify.Segmentation function and the FindSpatiallyVariableFeatures functions
GetTissueCoordinates.Segmentation <- function(object, full = TRUE, ...){
  a <- rbindlist(lapply(object@polygons , function(o){
    as.data.frame(t(o@labpt))
  }), id = "cell")
  b <- data.frame(x = a$V1, y = a$V2, cell = a$cell) # Columns need to be in this order
  row.names(b) <- b$cell
  return(b)
}

assignInNamespace("GetTissueCoordinates.Segmentation",
                  GetTissueCoordinates.Segmentation, ns = "SeuratObject")
################

##### Remove this line
# Just to reduce the size of the dataset for testing locally
#resolve.obj <- subset(resolve.obj, subset = nCount_Spatial > 1000)
########
resolve.obj <- ScaleData(resolve.obj) # 

sfs_method <- "markvariogram" # moransi also runs 
resolve.obj <- FindSpatiallyVariableFeatures(resolve.obj, assay = "Resolve", image = "seg",
  selection.method = sfs_method)
top_features <- SpatiallyVariableFeatures(resolve.obj, decreasing = T,
  selection.method = sfs_method)

top_fp <- ImageFeaturePlot(resolve.obj, fov = "centroids", features = "Vip", size = 0,
  molecules = top_features[1:3], nmols = nrow(transcripts)) 
