rm(list = ls())
library(Seurat)
library(future)
plan("multisession", workers = 10)
library(data.table)
library(RImageJROI)

LoadResolve <- function(cells.file, mask.file, image.file = NULL, ROI.file = NULL, transcript.file = NULL, sample.name = NULL){
    if (is.null(sample.name)){
        sample.name =  sub(".+/([^-]+)-.+", "\\1", mask.file)
    }
    mask_metadata <- readTIFF(mask.file, payload = FALSE)
    width <- mask_metadata$width
    height <- mask_metadata$length
    cells <- fread(cells.file)
    cells[, x := DAPI_HEIGHT - x] 
    cells[, x_tile := DAPI_HEIGHT - x_tile]
    cells[, cell_id := paste0(sample.name, "-cell-", label)]
    if (!is.null(image.file)){
        image <- readTIFF(image.file, native = FALSE, all = FALSE, convert = FALSE,
         info = FALSE, indexed = FALSE, as.is = FALSE, payload = TRUE)
    }
    if (!is.null(ROI.file)){
        cell_roi <- read.ijzip(ROI.file, names = TRUE, list.files = FALSE, verbose = FALSE)
        cell_roi <- lapply(cell_roi, `[[`, "coords")
        cell_roi <- lapply(cell_roi, as.data.table)
        cell_roi <- rbindlist(cell_roi, idcol = "cell_id")
        cell_roi[, cell_id := paste0(sample.name, "-cell-", as.integer(cell_id) + 1)]
        cell_roi[, x2 := y] # reorienting the coordinates
        cell_roi[, y := x]
        cell_roi[, x := height - x2]
        cell_roi[, x2 := NULL]
        cell_roi <- as.data.frame(cell_roi)
    }
    if (!is.null(transcript.file)){
        transcripts <- fread(transcript.file)
        transcripts <- transcripts[, .(x = height - V2, y = V1, gene = V4)]
        transcripts <- as.data.frame(transcripts)
    }
    
    markers <- colnames(cells)[!(colnames(cells) %in% c("sample_name", "group",
          "cell_id", "area", "centroid.x", "centroid.y", "label"))]
    cells[, c(markers) := lapply(.SD, function(x){ifelse(is.na(x), 0, x)}), .SDcols = markers]
    
    cell_coordinaets[]
    
    
    measurements <- cells[, c(markers), with = F]
    measurements <- t(as.matrix(measurements))

    row.names(measurements) <- markers
    colnames(measurements) <- cells$cell_id
    
    resolve.obj <- CreateSeuratObject(counts = measurements, assay = "Spatial")

    resolve.obj@meta.data["cell_id"] <- rownames(resolve.obj@meta.data)
    resolve.obj@meta.data["sample_name"] <- cells$sample.name
    resolve.obj@meta.data <- as.data.table(resolve.obj@meta.data)

    resolve.obj@meta.data[, region := 1]
    resolve.obj@meta.data[, z := 0]
    resolve.obj@meta.data[, tile_num := 0]

    resolve.obj@meta.data <- as.data.frame(resolve.obj@meta.data)
    rownames(resolve.obj@meta.data) <- resolve.obj@meta.data$cell_id 

    return(resolve.obj)
}
