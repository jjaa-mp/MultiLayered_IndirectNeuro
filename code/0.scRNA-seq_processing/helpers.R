# biomaRt conversion

## Input 1) rownames from count matrix and 2) biomart conversion dataframe

genenames_matrix <- function(matrix_genesID, biomart_conversion){
  genes_mtx = data.frame(ENS_id = matrix_genesID)
  colnames(conversion) <- c("ENS_id", "hgnc_symbol")
  merged_dataframes <- left_join(x = genes_mtx, y = conversion, by = "ENS_id")
  
  #Cleaning
  merged_dataframes <- merged_dataframes %>% 
    group_by(ENS_id) %>% 
    mutate(ENS_id_2 = paste0(hgnc_symbol, collapse = ";"))
  #If empty or NA, keep original Ensembl IDs
  merged_dataframes[which(merged_dataframes$ENS_id_2=="", arr.ind = TRUE),]$ENS_id_2 <- merged_dataframes[which(merged_dataframes$ENS_id_2=="", arr.ind = TRUE),]$ENS_id
  merged_dataframes[which(merged_dataframes$ENS_id_2=="NA", arr.ind = TRUE),]$ENS_id_2 <- merged_dataframes[which(merged_dataframes$ENS_id_2=="NA", arr.ind = TRUE),]$ENS_id
  #Remove duplicates from ENS_id
  merged_dataframes <- merged_dataframes[!duplicated(merged_dataframes$ENS_id),]
  #Keep gene names belonging to different Ensembl IDs
  merged_dataframes$ENS_id_2 <- make.unique(merged_dataframes$ENS_id_2)
  
  return(merged_dataframes$ENS_id_2)
}

# QC Filtering based on total counts, gene counts, percentage of mitochondrial gene and zero counts
## Inspired by Moreau et al 2021 study (https://doi.org/10.1242/dev.197962;https://matthieuxmoreau.github.io/EarlyPallialNeurogenesis/)

QC.Filtering <- function(seurat_object, var.split, nMAD = 3,
                         verbose = T, 
                         by.mito = F, 
                         by.count = F, 
                         by.gene = F,
                         by.zeros = F){
  list_seurats0 <- SplitObject(seurat_object, split.by = var.split)
  list_seurats <- lapply(list_seurats0, FUN = function(x) {
    x <- as.data.frame(x@meta.data)
  })
  cells.to.keep <- c()
  if (by.mito) {
    message(paste0(Sys.time(), ": Filtering genes by percent.mito."))
    for (x in list_seurats) {

      message(paste0("Selected variable: ", unique(x[[var.split]])))
      Cell.QC.Stat <- x
      max.mito.thr <- median(Cell.QC.Stat$nCount_RNA) + nMAD*mad(Cell.QC.Stat$percent.mito)
      min.mito.thr <- median(Cell.QC.Stat$nCount_RNA) - nMAD*mad(Cell.QC.Stat$percent.mito)
      Cell.QC.Stat <- Cell.QC.Stat %>% filter(nCount_RNA < max.mito.thr)
      Cell.QC.Stat <- Cell.QC.Stat %>% filter(nCount_RNA > min.mito.thr)
      cells.to.keep <- c(cells.to.keep, rownames(Cell.QC.Stat))
    }
  }
  if (by.count) {
    message(paste0(Sys.time(), ": Filtering cells by counts."))
    for (x in list_seurats) {

      message(paste0("Selected variable: ", unique(x[[var.split]])))
      Cell.QC.Stat <- x
      max.count.thr <- median(Cell.QC.Stat$nCount_RNA) + nMAD*mad(Cell.QC.Stat$nCount_RNA)
      min.count.thr <- median(Cell.QC.Stat$nCount_RNA) - nMAD*mad(Cell.QC.Stat$nCount_RNA)
      Cell.QC.Stat <- Cell.QC.Stat %>% filter(nCount_RNA < max.count.thr)
      Cell.QC.Stat <- Cell.QC.Stat %>% filter(nCount_RNA > min.count.thr)
      cells.to.keep <- c(cells.to.keep, rownames(Cell.QC.Stat))
    }
  }
  if (by.gene) {
    message(paste0(Sys.time(), ": Filtering cells by nFeatures"))
    for (x in list_seurats) {

      message(paste0("Selected variable: ", unique(x[[var.split]])))
      Cell.QC.Stat <- x
      max.feature.thr <- median(Cell.QC.Stat$nFeature_RNA) + nMAD*mad(Cell.QC.Stat$nFeature_RNA)
      min.feature.thr <- median(Cell.QC.Stat$nFeature_RNA) - nMAD*mad(Cell.QC.Stat$nFeature_RNA)
      Cell.QC.Stat <- Cell.QC.Stat %>% filter(nFeature_RNA < max.feature.thr)
      Cell.QC.Stat <- Cell.QC.Stat %>% filter(nFeature_RNA > min.feature.thr)
      cells.to.keep <- c(cells.to.keep, rownames(Cell.QC.Stat))
    }
  }
  if (by.zeros){
    message(paste0(Sys.time(), ": Filtering cells by percent.zeros"))
    for (x in list_seurats) {

      message(paste0("Selected variable: ", unique(x[[var.split]])))
      Cell.QC.Stat <- x
      max.zero.thr <- median(Cell.QC.Stat$percent.zeros) + nMAD*mad(Cell.QC.Stat$percent.zeros)
      min.zero.thr <- median(Cell.QC.Stat$percent.zeros) - nMAD*mad(Cell.QC.Stat$percent.zeros)
      Cell.QC.Stat <- Cell.QC.Stat %>% filter(percent.zeros < max.zero.thr)
      Cell.QC.Stat <- Cell.QC.Stat %>% filter(percent.zeros > min.zero.thr)
      cells.to.keep <- c(cells.to.keep, rownames(Cell.QC.Stat))
    }
  }
  
  if (length(list_seurats0) == 1){
    seurat.combined <- list_seurats0[[1]]
  }
  
  else{
    if (verbose) message(paste0(Sys.time(), " Merging ", length(list_seurats), " subdatasets."))
    seurat.combined <- merge(x= list_seurats0[[1]], y=list_seurats0[2:length(list_seurats)])
  }
  if (verbose) message(paste0(Sys.time(), " Nº Cells kept ", length(cells.to.keep[!duplicated(cells.to.keep)])))
  seurat.combined <- subset(seurat.combined, 
                            cells = cells.to.keep[!duplicated(cells.to.keep)])
  return(seurat.combined)
}

# Seurat SCT and Integration

sct_int_Seurat <- function(seurat_object, var.split){
  
  list_seurats <- SplitObject(seurat_object, split.by = var.split)
  list_seurats <- lapply(X = list_seurats, FUN = SCTransform, 
                         variable.features.n = 3000,
                         return.only.var.genes = FALSE)

  features <- SelectIntegrationFeatures(object.list = list_seurats, nfeatures = 3000)
  list_seurats <- PrepSCTIntegration(object.list = list_seurats, anchor.features = features)
  seurat.anchors <- FindIntegrationAnchors(object.list = list_seurats, 
                                           normalization.method = "SCT",
                                           anchor.features = features)
  seurat.combined <- IntegrateData(anchorset = seurat.anchors, 
                                   normalization.method = "SCT")
  return(seurat.combined)
}
