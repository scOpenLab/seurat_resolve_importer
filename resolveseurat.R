rm(list = ls())
library(Seurat)
library(future)
plan("multisession", workers = 10)
library(data.table)
library(RImageJROI)

# Source the Load functions
source("load_function.R")

# https://satijalab.org/seurat/articles/spatial_vignette_2.html#human-lymph-node-akoya-codex-system

### Path to files to read:
# PANORAMA_PATH: Path to the folder with the Panoramas from resolve
# OUTPUT_PATH: Path to the output of the Resolve Processing pipeline

PANORAMA_PATH <- ""
OUTPUT_PATH <- ""


resolve.objs <- LoadResolveFromFolders("/home/bq_mbortolomeazzi/NEUROBLASTOMA/output_mindagap",
    "/home/bq_mbortolomeazzi/NEUROBLASTOMA/panoramas_DAPI")


# SAMPLE NAME TO TEST
SAMPLE_NAME <- ""

resolve.obj <-resolve.objs[[SAMPLE_NAME]]


#### Quick test single cell analysis, not meaningful just to have something to plot
resolve.obj <- subset(resolve.obj, subset = nCount_Spatial > 10) # Skip cells with no transcripts 
resolve.obj <- SCTransform(resolve.obj, assay = "Resolve", verbose = FALSE)
resolve.obj <- FindVariableFeatures(resolve.obj, assay = "Resolve")
resolve.obj <- RunPCA(object = resolve.obj, npcs = 20, verbose = FALSE)
resolve.obj <- RunUMAP(object = resolve.obj, dims = 1:20, verbose = FALSE)
resolve.obj <- FindNeighbors(object = resolve.obj, dims = 1:20, verbose = FALSE)
resolve.obj <- FindClusters(object = resolve.obj, verbose = FALSE,
  resolution = 0.1)

resolve.obj@active.assay <- "Resolve" # To show the real counts in the plots




### Plotting
# FEATURES_TO_PLOT: vector of gene names
FEATURES_TO_PLOT <- ""

# No spatial data
nos_dp <- DimPlot(resolve.obj, label = TRUE, label.box = TRUE) + NoLegend() 
nos_fp <- FeaturePlot(resolve.obj, features = FEATURES_TO_PLOT, min.cutoff = "q10", max.cutoff = "q90") 

# Plus spatial data (cells as circles on their centroids)
cen_dp <- ImageDimPlot(resolve.obj, cols = "parade", fov = "centroids")
cen_fp <- ImageFeaturePlot(resolve.obj, fov = "centroids", features = FEATURES_TO_PLOT)

# Plus spatial data (cells as segmentation outlines)
seg_dp <- ImageDimPlot(resolve.obj, cols = "parade", fov = "segmentation")
seg_fp <- ImageFeaturePlot(resolve.obj, features = FEATURES_TO_PLOT,
    min.cutoff = "q10", max.cutoff = "q90", fov = "segmentation")

# Raw transcript visualization
# - size = 0 avoids displaying the cells
# - nmols = 9999999999 displays all the transcripts, instead of a sample 
cen_mp <- ImageFeaturePlot(resolve.obj, fov = "centroids", features = FEATURES_TO_PLOT[1],
  size = 0, molecules = FEATURES_TO_PLOT, mols.cols = c("#FF0000",
  "#00FF00", "#0000FF", "#FFFF00"), nmols = 9999999999)

seg_mp <- ImageFeaturePlot(resolve.obj, fov = "segmentation", features = FEATURES_TO_PLOT[1], 
  molecules = FEATURES_TO_PLOT, mols.cols = c("#0000FF"), nmols = 9999999999,
  min.cutoff = "q10", max.cutoff = "q90")




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
