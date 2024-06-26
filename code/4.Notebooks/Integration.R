# V1

library(Seurat)
library(dplyr)
set.seed(1) 

setwd(file.path('/home/rstudio/jm_rstudio/jm_2023_R/'))

d2.exp <- ReadMtx(
  mtx = "./GW18_V1/matrix.mtx", features = "./GW18_V1/genes.tsv",
  cells = "./GW18_V1/barcodes.tsv"
)
bh21 <- CreateSeuratObject(counts = d2.exp, project = 'BhaduriV1', min.cells = 3, assay='RNA')
bh21[["RNA"]] <- as(object = bh21[["RNA"]], Class = "Assay")

bh21[["percent.mt"]] <- PercentageFeatureSet(bh21, pattern = "^MT-")
bh21@meta.data[["percent.zeros"]] <- (colSums(GetAssayData(bh21) == 0)*100)/
  (dim(GetAssayData(bh21))[1])

source("/home/rstudio/jm_rstudio/data_CBL/indNeuro_tmp/0.scRNA-seq_processing/helpers.R")
bh21 <- QC.Filtering(bh21, var.split = "orig.ident", 
                     by.mito = TRUE, 
                     by.count = TRUE, 
                     by.gene = TRUE, 
                     by.zeros = TRUE)

bh21 <- SCTransform(bh21, vars.to.regress = "percent.mt", verbose = FALSE)


bh21 <- RunPCA(bh21, verbose = FALSE)
DimPlot(bh21, reduction = "pca")
ElbowPlot(bh21, ndims = 50)

bh21 <- RunUMAP(bh21, dims = 1:20, verbose = FALSE)
DimPlot(bh21, reduction = "umap", label = TRUE)


bh21 <- FindNeighbors(bh21, dims = 1:20, verbose = FALSE)
bh21 <- FindClusters(bh21,verbose = FALSE)
DimPlot(bh21, label = TRUE)

# find markers
markers <- FindAllMarkers(bh21, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

DoHeatmap(bh21, features = top5$gene) + NoLegend()
FeaturePlot(bh21, features = c('HOPX', 'VIM', 'EOMES', 'NEUROD2'))

bh21 <- SetIdent(bh21, value = bh21@meta.data$seurat_clusters)

Idents(bh21)
DimPlot(bh21, group.by = 'seurat_clusters', label = TRUE)
bh21.sub.v1 <- subset(x = bh21, idents = c("10", "6", "8"))

DefaultAssay(bh21.sub.v1) <- 'RNA'

bh21.sub.v1 <- DietSeurat(
  bh21.sub.v1,
  layers = 'counts',
  features = NULL,
  assays = 'RNA',
  dimreducs = NULL,
  graphs = NULL,
  misc = TRUE)
gc()

bh21.sub.v1

# PFC
exp <- ReadMtx(
  mtx = "./GW18_PFC/matrix.mtx", features = "./GW18_PFC/enes.tsv",
  cells = "./GW18_PFC/barcodes.tsv"
)
bh21 <- CreateSeuratObject(counts = exp, project = 'BhaduriPFC', min.cells = 3, assay='RNA')

bh21[["RNA"]] <- as(object = bh21[["RNA"]], Class = "Assay")
bh21[["percent.mt"]] <- PercentageFeatureSet(bh21, pattern = "^MT-")
bh21@meta.data[["percent.zeros"]] <- (colSums(GetAssayData(bh21) == 0)*100)/(dim(GetAssayData(bh21))[1])

source("/home/rstudio/jm_rstudio/data_CBL/indNeuro_tmp/0.scRNA-seq_processing/helpers.R")
bh21 <- QC.Filtering(bh21, var.split = "orig.ident", 
                     by.mito = TRUE, 
                     by.count = TRUE, 
                     by.gene = TRUE, 
                     by.zeros = TRUE)

bh21 <- SCTransform(bh21, vars.to.regress = "percent.mt", verbose = FALSE)


bh21 <- RunPCA(bh21, verbose = FALSE)
DimPlot(bh21, reduction = "pca")
ElbowPlot(bh21, ndims = 50)

bh21 <- RunUMAP(bh21, dims = 1:20, verbose = FALSE)
DimPlot(bh21, reduction = "umap")


bh21 <- FindNeighbors(bh21, dims = 1:20, verbose = FALSE)
bh21 <- FindClusters(bh21,verbose = FALSE)
DimPlot(bh21, label = TRUE)

# find markers
PFC.markers <- FindAllMarkers(bh21, only.pos = TRUE)
PFC.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

DoHeatmap(bh21, features = top5$gene) + NoLegend()
FeaturePlot(bh21, features = c('HOPX', 'VIM', 'EOMES', 'NEUROD2'))

bh21 <- FindSubCluster(
  bh21,
  cluster = '5',
  graph.name = 'SCT_snn',
  subcluster.name = "sub.cluster",
  resolution = 0.5)

bh21 <- SetIdent(bh21, value = bh21@meta.data$sub.cluster)
PFC.markers <- FindAllMarkers(bh21, only.pos = TRUE)
PFC.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

DoHeatmap(bh21, features = top5$gene) + NoLegend()
DimPlot(bh21, label = TRUE)
FeaturePlot(bh21, features = c('PAX6', 'SOX2', 'EOMES', 
                               'SSTR2', 'PPP1R17'),label = FALSE)

FeaturePlot(bh21, features = c('EOMES', 'SSTR2', 'PPP1R17'))

Idents(bh21)
DimPlot(bh21, group.by = 'sub.cluster', label = TRUE)

bh21.sub.pfc <- subset(x = bh21, idents = c("8", "5_1"))

DefaultAssay(bh21.sub.pfc) <- 'RNA'

bh21.sub.pfc <- DietSeurat(
  bh21.sub.pfc,
  layers = 'counts',
  features = NULL,
  assays = 'RNA',
  dimreducs = NULL,
  graphs = NULL,
  misc = TRUE)
gc()

bh21.sub.pfc

bh21.sub.v1

combined <- merge(bh21.sub.pfc, y = bh21.sub.v1, add.cell.ids = c("pfc", "v1"), project = "Bhaduri21")


# Trevino 2021
tr21 <- readRDS('/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/tr21_raw_with_QC_rna_counts.rds')
tr21
combined <- subset(x = combined, subset = TOP2A < 1)

s.list <- c(tr21, combined)

s.list <- lapply(X = s.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = s.list, nfeatures = 3000)
s.list <- PrepSCTIntegration(object.list = s.list, anchor.features = features)
s.list <- lapply(X = s.list, FUN = RunPCA, features = features)
anchors <- FindIntegrationAnchors(object.list = s.list, normalization.method = "SCT",
                                  anchor.features = features, 
                                  dims = 1:30, reduction = "rpca", k.anchor = 20)
combined.sct <- IntegrateData(anchorset = anchors, 
                              normalization.method = "SCT", dims = 1:30)

combined.sct <- RunPCA(combined.sct, verbose = FALSE)
DimPlot(combined.sct)

combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:30)
combined.sct <- FindClusters(combined.sct, resolution = 0.1) #Broad clustering
DimPlot(combined.sct)
DimPlot(combined.sct, group.by = 'integrated_snn_res.0.1')

markers <- FindAllMarkers(combined.sct, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
DoHeatmap(combined.sct, features = top5$gene) + NoLegend()

# Unknown cluster 3
Idents(combined.sct)
combined.sct <- subset(x = combined.sct, idents = c("3"), invert=TRUE)

# Name clusters
combined.sct@meta.data <- mutate(combined.sct@meta.data,
                                 Cluster.Joint = case_when(
                                   integrated_snn_res.0.1 == "1" ~ "oRG",
                                   integrated_snn_res.0.1 == "2" ~ "vRG",
                                   integrated_snn_res.0.1  == "0" ~ "IP"))


DimPlot(combined.sct, group.by = 'Cluster.Joint')
saveRDS(combined.sct, file = "/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/combined_integrated-all.rds")
combined.sct <- readRDS(file = "/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/combined_integrated-all.rds")

# Aesthetics
copy_seurat <- combined.sct
copy_seurat@reductions$pca@cell.embeddings[,1] = copy_seurat@reductions$pca@cell.embeddings[,1]*-1
copy_seurat@reductions$pca@cell.embeddings[,2] = copy_seurat@reductions$pca@cell.embeddings[,2]*-1
DimPlot(copy_seurat, reduction = "pca", group.by = "Cluster.Joint")

write.table(combined.sct@reductions$pca@cell.embeddings[,0:2], 
            file = "/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/combined_all_cellembeddings.tsv",
            row.names = TRUE, col.names = TRUE,  sep = '\t')

library(SeuratDisk)

DefaultAssay(combined.sct) <- 'RNA'
combined.sct <- DietSeurat(combined.sct, assays  = 'RNA')

SaveH5Seurat(combined.sct, filename = "/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/combined_count_all.h5Seurat")

Convert("/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/combined_count_all.h5Seurat", dest = "h5ad")
file.remove("/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/combined_count_all.h5Seurat") 
