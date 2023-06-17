# Libraries
library(Seurat)
library(dplyr)
library(SeuratDisk)

# Helpers
source("/home/rstudio/jm_rstudio/data_CBL/indNeuro_tmp/0.scRNA-seq_processing/helpers.R")

# Read raw counts - Data downloaded from http://solo.bmap.ucla.edu/shiny/webapp/
load("~/jm_rstudio/pol2019/raw_counts_mat.rdata") 


# Read meta data
meta = read.table("~/jm_rstudio/pol2019/cell_metadata.csv", 
                  sep = ",",
                  header = TRUE)
meta <- meta[,c(1:8)] #I'd take Donor as batch and Index as Sample.ID
rownames(meta) <- meta$Cell
meta$Cell <- NULL
table(meta$Cluster)

# Subset for neural progenitor cells
meta <- meta[meta$Cluster %in% c("oRG", "vRG", "IP"),]
colremove <- setdiff(colnames(raw_counts_mat),rownames(meta))
raw_counts_mat1 <- raw_counts_mat[, !colnames(raw_counts_mat) %in% colremove]

# Create Seurat object with min filtering
pol19 <- CreateSeuratObject(counts = raw_counts_mat1, project = "pol19",
                             assay = "RNA",
                             meta.data = meta)


pol19[["percent.mt"]] <- PercentageFeatureSet(pol19, pattern = "^MT-")
pol19@meta.data[["percent.zeros"]] <- (colSums(GetAssayData(pol19) == 0)*100)/(dim(GetAssayData(pol19))[1])
pol19 <- subset(x = pol19, subset = TOP2A < 1)

columns.keep <- c("orig.ident", "Cluster", "Subcluster", "Donor", 
                  "Layer", "Gestation_week", "Library", "percent.mt", "percent.zeros")
pol19@meta.data <- pol19@meta.data[,columns.keep]
pol19

## QC
source("~/jm_rstudio/data_CBL/sc_codes/project_layeringNCX/src/modules/QCFiltering.R")
pol19 <- QC.Filtering(pol19, var.split = "Library", 
                       by.mito = TRUE, 
                       by.count = TRUE, 
                       by.gene = TRUE, 
                       by.zeros = TRUE)

# Simple normalization
DefaultAssay(pol19) <- "RNA"

### --- ###

#dir.create("/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files")
#SaveH5Seurat(pol19, filename = "/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/pol19_counts.h5Seurat")
#Convert("/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/pol19_counts.h5Seurat", dest = "h5ad")
#file.remove("/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/pol19_counts.h5Seurat")

### --- ###


# Align with Trevino et al 21 dataset
tr21.rna <- readRDS(file = "/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/tr21_raw_with_QC_rna_counts.rds")
tr21.rna@meta.data$orig.ident = "Trevino"
pol19@meta.data$orig.ident = "Polioudakis"
head(pol19@meta.data,2)
head(tr21.rna@meta.data,2)
keep.columns1 <- c("orig.ident", "nCount_RNA",
"nFeature_RNA", "percent.mt",
"Cluster", "Gestation_week",  "Library")
keep.columns2 <- c("orig.ident", "nCount_RNA",
"nFeature_RNA","percent.mt",
"Cluster.Name", "Age", "Batch")
pol19@meta.data <- pol19@meta.data[,keep.columns1]
tr21.rna@meta.data <- tr21.rna@meta.data[,keep.columns2]
colnames(tr21.rna@meta.data)
colnames(pol19@meta.data)[5] <- "Cluster.Name"
colnames(pol19@meta.data)[6] <- "Age"
colnames(pol19@meta.data)[7] <- "Batch"
colnames(pol19@meta.data) == colnames(tr21.rna@meta.data)
genes.keep <- intersect(rownames(pol19), rownames(tr21.rna))
length(genes.keep)
tr21.rna <- subset(tr21.rna, features = genes.keep)
pol19 <- subset(pol19, features = genes.keep)
tmp <- merge(tr21.rna, y = pol19, add.cell.ids = c("tr21", "pol19"), project = "merged")
table(tmp@meta.data$Age)
tmp <- DietSeurat(tmp)
gc()
merged <- CreateSeuratObject(counts = tmp@assays$RNA@counts, project = "merged",
assay = "RNA",
meta.data = tmp@meta.data)
table(merged@meta.data$Age)
Idents(merged) <- "orig.ident"
#Quality Checks
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), split.by = "Age",
ncol = 3)
plot1 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
integrated <- sct_int_Seurat(seurat_object = merged, var.split = "Batch")
# PCA
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
DimPlot(integrated, reduction = "pca")
DimPlot(integrated, reduction = "pca", split.by = "orig.ident", group.by = "Cluster.Name")

VizDimLoadings(integrated, dims = 1:2, reduction = "pca")


# For better geometry
copy_seurat <- integrated
copy_seurat@reductions$pca@cell.embeddings[,1] = copy_seurat@reductions$pca@cell.embeddings[,1]*-1
copy_seurat@reductions$pca@cell.embeddings[,2] = copy_seurat@reductions$pca@cell.embeddings[,2]*-1

DimPlot(copy_seurat, reduction = "pca", split.by = "orig.ident", group.by = "Cluster.Name")
copy_seurat@reductions$pca@feature.loadings[,1] = copy_seurat@reductions$pca@feature.loadings[,1]*-1
copy_seurat@reductions$pca@feature.loadings[,2] = copy_seurat@reductions$pca@feature.loadings[,2]*-1

VizDimLoadings(copy_seurat, dims = 1:2, reduction = "pca")

sub <- subset(x = copy_seurat, idents = "Polioudakis")
plot1.pol19 <- DimPlot(sub, reduction = "pca", group.by = "Cluster.Name")+labs(title = "")

### --- ###
#dir.create("/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files")
#write.table(integrated@reductions$pca@cell.embeddings[,0:2],
#file = "/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/integrated_polioudakis19_cellembeddings.tsv",
#row.names = TRUE, col.names = TRUE,  sep = '\t')
### --- ###



############## Analysis Polioudakis et al 2019 independently ###################
load("~/jm_rstudio/pol2019/raw_counts_mat.rdata") 


# Read meta data
meta = read.table("~/jm_rstudio/pol2019/cell_metadata.csv", 
                  sep = ",",
                  header = TRUE)
meta <- meta[,c(1:8)] #I'd take Donor as batch and Index as Sample.ID
rownames(meta) <- meta$Cell
meta$Cell <- NULL
table(meta$Cluster)

# Subset for neural progenitor cells
meta <- meta[meta$Cluster %in% c("oRG", "vRG", "IP"),]
colremove <- setdiff(colnames(raw_counts_mat),rownames(meta))
raw_counts_mat1 <- raw_counts_mat[, !colnames(raw_counts_mat) %in% colremove]

# Create Seurat object with min filtering
pol19 <- CreateSeuratObject(counts = raw_counts_mat1, project = "pol19",
                            assay = "RNA",
                            meta.data = meta,
                            min.cells = 200,
                            min.features = 500)


pol19[["percent.mt"]] <- PercentageFeatureSet(pol19, pattern = "^MT-")
pol19@meta.data[["percent.zeros"]] <- (colSums(GetAssayData(pol19) == 0)*100)/(dim(GetAssayData(pol19))[1])
pol19 <- subset(x = pol19, subset = TOP2A < 1)

columns.keep <- c("orig.ident", "Cluster", "Subcluster", "Donor", 
                  "Layer", "Gestation_week", "Library", "percent.mt", "percent.zeros")
pol19@meta.data <- pol19@meta.data[,columns.keep]
pol19

## QC
source("~/jm_rstudio/data_CBL/sc_codes/project_layeringNCX/src/modules/QCFiltering.R")
pol19 <- QC.Filtering(pol19, var.split = "Library", 
                      by.mito = TRUE, 
                      by.count = TRUE, 
                      by.gene = TRUE, 
                      by.zeros = TRUE)

# Simple normalization
DefaultAssay(pol19) <- "RNA"

pol19 <- NormalizeData(
  pol19,
  normalization.method = "LogNormalize",
  scale.factor = 10000)

pol19 <- FindVariableFeatures(pol19)
LabelPoints(
  plot = VariableFeaturePlot(pol19),
  points = head(VariableFeatures(pol19), 10),
  repel = TRUE)

pol19 <- ScaleData(pol19)
pol19 <- RunPCA(pol19)

DimPlot(pol19, reduction = "pca")
DimPlot(pol19, group.by = "Cluster", reduction = "pca")
#For better geometry
copy_2 <- pol19
#copy_2@reductions$pca@cell.embeddings[,1] = copy_2@reductions$pca@cell.embeddings[,1]*-1
copy_2@reductions$pca@cell.embeddings[,2] = copy_2@reductions$pca@cell.embeddings[,2]*-1

DimPlot(copy_2, reduction = "pca", group.by = "Cluster")
library(ggplot2)
plot2.pol19 <- DimPlot(copy_2, group.by = "Cluster", reduction = "pca")+NoLegend()+labs(title = "")

plot2.pol19+plot1.pol19

VizDimLoadings(pol19, dims = 1:2, reduction = "pca")

FeaturePlot(pol19, features=c("EGR1", "F0S", "PTN", "HOPX", "EOMES", "VIM"))
FeaturePlot(pol19, features=c("PPFIA2", "RBFOX2", "AKAP9", "TCF4", "GRIA2", "CHD3", "ROBO2"))

ElbowPlot(pol19, ndims = 50)

pol19 <- RunUMAP(pol19, reduction = "pca", dims = 1:10)
umap_age <- DimPlot(pol19, group.by = "Gestation_week")
umap_batch <- DimPlot(pol19, group.by = "Library")
umap_age + umap_batch
DimPlot(pol19, group.by = "Subcluster")

################################################################################ 


