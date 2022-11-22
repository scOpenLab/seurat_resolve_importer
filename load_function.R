FormatResolve <- function(sample.name, mask.file, cells.file, transcript.file = NULL,
    image.file = NULL, roi.file = NULL){
    #' Format Resolve Data for Input
    #' 
    #' @description Helper function for ReadResolve, This function reads the provided file names
    #' and returns a list of formatted data for loading a SeuratObject.
    #'   
    #' @param sample.name character. Name of the sample (used for naming the cells, see @details).
    #' @param mask.file character. Path to the segmentation mask.
    #' @param cells.file character. Path to the single cell data.
    #' @param transcript.file character. Path to the transcript coordinates. Optional.
    #' @param image.file character. Path to an additional image to load. Optional.
    #' @param roi.file character. Path to the zip file with the segmentation ROIs. Optional.
    #' 
    #' @usage Meant to be called by ReadResolve. It assumes the input files have been created with
    #' ResolveTools pipeline: https://codebase.helmholtz.cloud/resolve_tools/resolve-processing
    #'
    #' @return A named list with the following objects:
    #' \itemize{
    #' \item expression: Expression matrix rows -> gene, columns -> cells.
    #' \item coordinates: Dataframe with the coordinates of the cell centroids.
    #' \item transcripts: Dataframe with the coordinates of the transcripts.
    #' \item image: Additional image
    #' \item rois: Dataframe with the coordinates of the segmentation ROIs.
    #' }
    #'
    #' @details The cells are named with the "$SAMPLE-cell-$N" pattern,
    #' with 1 < N < max(segmentation mask).
    #'
    outs <- list()
    
    mask_metadata <- readTIFF(mask.file, payload = FALSE)
    width <- mask_metadata$width
    height <- mask_metadata$length
    
    cells <- fread(cells.file)
    cells[, x := height - centroid.y]
    cells[, y := centroid.x]
    cells[, cell_id := paste0(sample.name, "-cell-", label)]
    
    markers <- colnames(cells)[!(colnames(cells) %in% c("sample_name", "group",
          "cell_id", "area", "centroid.x", "centroid.y", "label", "cell", "x", "y"))]
    cells[, c(markers) := lapply(.SD, function(x){ifelse(is.na(x), 0, x)}), .SDcols = markers]
    
    fov_coordinates <- as.data.frame(cells[, .(x, y, cell = cell_id)])
    rownames(fov_coordinates) <- fov_coordinates$cell

    measurements <- cells[, c(markers), with = F]
    measurements <- t(as.matrix(measurements))

    row.names(measurements) <- markers
    colnames(measurements) <- cells$cell_id
    
    outs[["expression"]] <- measurements
    outs[["coordinates"]] <- fov_coordinates
    
    if (!is.null(transcript.file)){
        transcripts <- fread(transcript.file)
        transcripts <- transcripts[, .(x = height - V2, y = V1, gene = V4)]
        transcripts <- as.data.frame(transcripts)
        outs[["transcripts"]] <- transcripts
    }

    if (!is.null(image.file)){
        image <- readTIFF(image.file, native = FALSE, all = FALSE, convert = FALSE,
         info = FALSE, indexed = FALSE, as.is = FALSE, payload = TRUE)
        outs[["image"]] <- image
    }
    
    if (!is.null(roi.file)){
        cell_roi <- read.ijzip(roi.file, names = TRUE, list.files = FALSE, verbose = FALSE)
        cell_roi <- lapply(cell_roi, `[[`, "coords")
        cell_roi <- lapply(cell_roi, as.data.table)
        cell_roi <- rbindlist(cell_roi, idcol = "cell_id")
        cell_roi[, cell_id := paste0(sample.name, "-cell-", as.integer(cell_id) + 1)]
        cell_roi[, y2 := y]
        cell_roi[, y := x]
        cell_roi[, x := height - y2]
        cell_roi[, y2 := NULL]
        cell_roi <- as.data.frame(cell_roi)
        outs[["rois"]] <- cell_roi
    }
    return(outs)
}    

ReadResolve <- function(data.dir, sample.name = NULL){
    #' Read Resolve Data from a given folder and return them formatted for import with Seurat
    #' 
    #' @description This function reads the provided file names from a given folder
    #' and returns a list of formatted data for loading a SeuratObject.
    #'
    #' @param data.dir character. Path to the directory containing the data to be loaded.
    #' @param sample.name character. Name of the sample (used for naming the cells, see @details).
    #' Optional, if not provided, it is generated from the `data.dir` parameter.
    #' 
    #' @usage It assumes the input files in `data.dir` have been created with:
    #' ResolveTools pipeline: https://codebase.helmholtz.cloud/resolve_tools/resolve-processing
    #'
    #' @return A named list with the following objects:
    #' \itemize{
    #' \item expression: Expression matrix rows -> genes, columns -> cells.
    #' \item coordinates: Dataframe with the coordinates of the cell centroids.
    #' \item transcripts: Dataframe with the coordinates of the transcripts.
    #' \item image: Additional image
    #' \item rois: Dataframe with the coordinates of the segmentation ROIs.
    #' }
    #'
    #' @details The cells are named with the "$SAMPLE-cell-$N" pattern,
    #' with 1 < N < max(segmentation mask).
    #'
    metadata <- dir(data.dir, recursive = FALSE, full.names = TRUE)
    metadata <- data.table(file_name = metadata)
    metadata[, sample_name := sub(".+/(.+)-.+", "\\1", file_name)]
    metadata[, file_type := sub(".+-(.+)\\..+", "\\1", file_name)]
    metadata <- dcast(metadata, sample_name ~ file_type, value.var = "file_name")
    setnames(metadata, grep("mask", colnames(metadata), value = TRUE), "mask")
    
    if (!is.null(sample.name)){
        metadata$sample_name <- sample.name
    }
    
    outs <- FormatResolve(metadata$sample_name, metadata$mask, metadata$cell_data,
        metadata$filtered_transcripts, metadata$gapfilled, metadata$roi)
    return(outs)
}


LoadResolve <- function(data.dir, assay = "Resolve", sample.name = NULL){
    #' Read Resolve Data from a given folder and create a SeuratObject.
    #' 
    #' @description This function reads the provided file names from a given folder
    #' and returns a list of formatted data for loading a SeuratObject.
    #'
    #' @param data.dir character. Path to the directory containing the data to be loaded.
    #' @param sample.name character. Name of the sample (used for naming the cells, see @details).
    #' Optional, if not provided, it is generated from the `data.dir` parameter.
    #' 
    #' @usage It assumes the input files in `data.dir` have been created with:
    #' ResolveTools pipeline: https://codebase.helmholtz.cloud/resolve_tools/resolve-processing
    #'
    #' @return A SeuratObject with:
    #' \itemize{
    #' \item assays = assay with @data containing the matrix with the transcript counts in each cell.
    #' \item images = FOVs for centroids and, if available, segmentation ROIs.  
    #'    If the trancritps are available they are loaded in the molecules slot.
    #' \item  meta.data
    #' }
    #'
    #' @details The cells are named with the "$sample.name-cell-$N" pattern.
    #'
    
    data <- ReadResolve(data.dir, sample.name)
    
    resolve.obj <- CreateSeuratObject(counts = data[["expression"]], assay = assay)

    resolve.obj@meta.data["cell_id"] <- rownames(resolve.obj@meta.data)
    resolve.obj@meta.data["sample_name"] <- sample.name
    resolve.obj@meta.data <- as.data.table(resolve.obj@meta.data)

    resolve.obj@meta.data[, region := 1]
    resolve.obj@meta.data[, z := 0]
    resolve.obj@meta.data[, tile_num := 0]

    resolve.obj@meta.data <- as.data.frame(resolve.obj@meta.data)
    rownames(resolve.obj@meta.data) <- resolve.obj@meta.data$cell_id 
    
    resolve.obj@images <- list()
    resolve.obj@images["centroids"] <- CreateFOV(
        coords = data[["coordinates"]],
        type = "centroids",
        nsides = 0L,
        radius = 1L,
        theta = 0L,
        molecules = data[["transcripts"]],
        assay = assay,
        key = NULL,
        name = NULL)
        
    if(!is.null(data[["rois"]])){
        resolve.obj@images["segmentation"] <- CreateFOV(
            coords = data[["rois"]],
            type = "segmentation",
            nsides = 0L,
            radius = 1L,
            theta = 0L,
            molecules = data[["transcripts"]], 
            assay = assay,
            key = NULL,
            name = NULL)
    }   
    return(resolve.obj)
}
