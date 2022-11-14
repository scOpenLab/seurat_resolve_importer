LoadResolve <- function(cells.file, mask.file, image.file = NULL, ROI.file = NULL, transcript.file = NULL, sample.name = NULL){
    if (is.null(sample.name)){
        sample.name =  sub(".+/([^-]+)-.+", "\\1", mask.file)
    }
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
    
    resolve.obj <- CreateSeuratObject(counts = measurements, assay = "Spatial")

    resolve.obj@meta.data["cell_id"] <- rownames(resolve.obj@meta.data)
    resolve.obj@meta.data["sample_name"] <- cells$sample.name
    resolve.obj@meta.data <- as.data.table(resolve.obj@meta.data)

    resolve.obj@meta.data[, region := 1]
    resolve.obj@meta.data[, z := 0]
    resolve.obj@meta.data[, tile_num := 0]

    resolve.obj@meta.data <- as.data.frame(resolve.obj@meta.data)
    rownames(resolve.obj@meta.data) <- resolve.obj@meta.data$cell_id 
    
    
    if (!is.null(transcript.file)){
        transcripts <- fread(transcript.file)
        transcripts <- transcripts[, .(x = height - V2, y = V1, gene = V4)]
        transcripts <- as.data.frame(transcripts)
    }
    else{
        transcripts <- NULL
    }
    
    resolve.obj@images <- list()
    resolve.obj@images["cen"] <- CreateFOV(
        coords = fov_coordinates,
        type = "centroids",
        nsides = 0L,
        radius = 1L,
        theta = 0L,
        molecules = transcripts, # Only for visualising the raw transcripts, can be skipped
        assay = "Spatial",
        key = NULL,
        name = NULL)

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
        cell_roi[, y2 := y]
        cell_roi[, y := x]
        cell_roi[, x := height - y2]
        cell_roi[, y2 := NULL]
        cell_roi <- as.data.frame(cell_roi)
        
        resolve.obj@images["seg"] <- CreateFOV(
        coords = cell_roi,
        type = "segmentation",
        nsides = 0L,
        radius = 1L,
        theta = 0L,
        molecules = transcripts, # Only for visualising the raw transcripts, can be skipped
        assay = "Spatial",
        key = NULL,
        name = NULL)
    }
    return(resolve.obj)
}

LoadResolveFromFolders <- function(cell.folder, panorama.folder = NULL){
    metadata <- dir(cell.folder, recursive = T, full.names = T, pattern = "EXP")
    metadata <- data.table(file_name = metadata)
    metadata[, sample_name := sub(".+/(.+)-.+", "\\1", file_name)]
    metadata[, file_type := sub(".+/.+-(.+)\\..+", "\\1", file_name)]
    if (!is.null(panorama.folder)){
        panorama_meta <- dir(panorama.folder, recursive = TRUE, full.names = TRUE,
            pattern = "(_Channel3_R8_.tiff|results_withFP.txt)")
        panorama_meta <- data.table(file_name = panorama_meta)
        panorama_meta[, sample_name := sub(".+Panorama_(.+)_(Channel3_R8_.tiff|results_withFP.txt)", "\\1", file_name)]
        panorama_meta[, file_type := ifelse(grepl("txt", file_name), "transcripts", "DAPI")]
        metadata <- rbind(metadata, panorama_meta)
    }
    metadata <- dcast(metadata, sample_name ~ file_type, value.var = "file_name")
    sample.names <- sort(unique(metadata$sample_name))
    resolve.objs <-lapply(sample.names, function(sname){
        mdata <- metadata[sample_name == sname]
        if (!is.null(panorama.folder)){
            obj <- LoadResolve(mdata$cell_data, mdata$mask, mdata$gridfilled, mdata$roi, mdata$transcripts, sname)
        } else {
            obj <- LoadResolve(mdata$cell_data, mdata$mask, mdata$gridfilled, mdata$roi, sample.name = sname)
        }
        print(sname)
        return(obj)
    })
    names(resolve.objs) <- sample.names
    return(resolve.objs)
}
