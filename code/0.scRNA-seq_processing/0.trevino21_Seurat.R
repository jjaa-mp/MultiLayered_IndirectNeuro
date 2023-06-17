# Libraries
library(Seurat)
library(readxl)
library(tibble)
library(dplyr)
library(biomaRt)
library(SeuratDisk)
library(factoextra)

# Helpers
source("/home/rstudio/jm_rstudio/data_CBL/indNeuro_tmp/0.scRNA-seq_processing/helpers.R")

# Read raw data from tsv.gz file - obligatory if first time running this script

#mtx <- read.table("/home/rstudio/jm_rstudio/data_indNeuro/rna_counts.tsv.gz", sep='\t')
#saveRDS(mtx, file = "/home/rstudio/jm_rstudio/data_indNeuro/rna_counts.rds")

# Read raw data from rds file
mtx <- readRDS("/home/rstudio/jm_rstudio/data_indNeuro/rna_counts.rds")

# Read meta data
cell.meta <- read.table("/home/rstudio/jm_rstudio/data_indNeuro/rna_cell_metadata.txt", header = TRUE, sep='\t')

# Remove doublets using metadata
cell.meta <- cell.meta[cell.meta$DF_classification != "Doublet",]
cell.meta = cell.meta[cell.meta$Cell.ID %in% colnames(mtx),] #Cells metadata in matrix  

# Read supplementary data to include cluster annotation
cluster_annotation <- read_excel("/home/rstudio/jm_rstudio/data_indNeuro/ST_trevino21_rna.xlsx", sheet = "F")
colnames(cluster_annotation) <- cluster_annotation[1,]
cluster_annotation <- as_tibble(cluster_annotation)
cluster_annotation <- cluster_annotation[-1,]

# Select clusters of interested based on original annotation
cluster_annotation <- cluster_annotation[cluster_annotation$Name %in% c("Early RG", "Late RG", "nIPC"),]

cell.meta <- cell.meta[cell.meta$seurat_clusters %in% cluster_annotation$`Cluster ID`,]
cell.meta <-  cell.meta[cell.meta$Cell.ID %in% colnames(mtx),]
rownames(cell.meta) <- cell.meta$Cell.ID
cell.meta$Cell.ID <- NULL

mtx <- mtx[,colnames(mtx) %in% rownames(cell.meta)]
#nrow(cell.meta) == ncol(mtx)

cell.meta <- cell.meta[,0:7]
names(cell.meta)[7]  <- "tr21_clusters"

cluster_annotation <- cluster_annotation[,0:2]
names(cluster_annotation) <- c("tr21_clusters", "Cluster.Name")
cell.meta$Cell.ID <- rownames(cell.meta)
cell.meta <- inner_join(x=cell.meta,y=cluster_annotation,by="tr21_clusters")
rownames(cell.meta) <- cell.meta$Cell.ID
cell.meta$Cell.ID <- NULL


#Changing gene names from ENS in count matrix 

ensembl <- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                            dataset = 'hsapiens_gene_ensembl',
                            #mirror = "www"
                            host = 'https://grch37.ensembl.org')
conversion <- biomaRt::getBM(filters= "ensembl_gene_id", 
                             attributes= c("ensembl_gene_id","hgnc_symbol"),
                             values = rownames(mtx), 
                             mart=ensembl)

##Change Ensembl ID to HGNC symbol in rownames count matrix
rownames(mtx) <- genenames_matrix(rownames(mtx), conversion)

# Create Seurat object with min filtering
tr21 <- CreateSeuratObject(counts = mtx, project = "tr21", 
                           assay = "RNA", 
                           meta.data = cell.meta)

tr21[["percent.mt"]] <- PercentageFeatureSet(tr21, pattern = "^MT-")
tr21@meta.data[["percent.zeros"]] <- (colSums(GetAssayData(tr21) == 0)*100)/(dim(GetAssayData(tr21))[1])

Idents(tr21) <- "Age"
tr21 <- subset(x = tr21, idents = c("pcw16", "pcw20", "pcw21"))
tr21.rna <- subset(x = tr21, subset = TOP2A < 1)
columns.keep <- c("orig.ident", "Sample.ID", "Age", "Batch", "Cluster.Name", "percent.mt", "percent.zeros")
tr21.rna@meta.data <- tr21.rna@meta.data[,columns.keep]
tr21.rna

# QC
tr21.rna <- QC.Filtering(tr21.rna, var.split = "Batch", 
                         by.mito = TRUE, 
                         by.count = TRUE, 
                         by.gene = TRUE, 
                         by.zeros = TRUE)


DefaultAssay(tr21.rna) <- "RNA"

### --- ###
#dir.create("/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files")
#saveRDS(tr21.rna, file = "/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/tr21_raw_with_QC_rna_counts.rds")

#SaveH5Seurat(tr21.rna, filename = "/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/rna_counts.h5Seurat")
#Convert("/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/rna_counts.h5Seurat", dest = "h5ad")
#file.remove("/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/rna_counts.h5Seurat")

### --- ###

# Simple normalization

tr21.rna <- NormalizeData(
  tr21.rna,
  normalization.method = "LogNormalize",
  scale.factor = 10000)

tr21.rna <- FindVariableFeatures(tr21.rna)
LabelPoints(
  plot = VariableFeaturePlot(tr21.rna),
  points = head(VariableFeatures(tr21.rna), 10),
  repel = TRUE)

tr21.rna <- ScaleData(tr21.rna)
tr21.rna <- RunPCA(tr21.rna)
DimPlot(tr21.rna, reduction = "pca")

FeaturePlot(tr21.rna, features=c("CLU", "ID4", "PTN", "HOPX", "EOMES", "VIM"))
ElbowPlot(tr21.rna, ndims = 50)
tr21.rna <- RunUMAP(tr21.rna, reduction = "pca", dims = 1:20)
umap_age <- DimPlot(tr21.rna, group.by = "Age")
umap_batch <- DimPlot(tr21.rna, group.by = "Batch")
umap_age + umap_batch


# SCT per Batch and Integration by anchors
tr21.rna <- sct_int_Seurat(seurat_object = tr21.rna, var.split = "Batch")

# PCA
tr21.rna <- RunPCA(tr21.rna, npcs = 50, verbose = FALSE)
DimPlot(tr21.rna, reduction = "pca")
DimPlot(tr21.rna, reduction = "pca", split.by = "Age")

Loadings(tr21.rna[["pca"]])
VizDimLoadings(tr21.rna, dims = 1:2, reduction = "pca")

### --- ### - Cell embeddings

#write.table(tr21.rna@reductions$pca@cell.embeddings[,0:2], 
#          file = "/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/cellembeddings.tsv",
#          row.names = TRUE, col.names = TRUE,  sep = '\t')


### --- ###

### --- ### - Gene loadings

#write.table(tr21.rna@reductions$pca@feature.loadings[,0:2], 
#            file = "/home/rstudio/jm_rstudio/data_indNeuro/intermediate_files/featureloadings.tsv",
#            row.names = TRUE, col.names = TRUE,  sep = '\t')


### --- ###



##### prcomp PLOT ##### 
py.mtx <- t(x = tr21.rna@assays$SCT@scale.data)[,ncol(t(x = tr21.rna@assays$SCT@scale.data)):1] #for better geometry

pca.res <-  prcomp(py.mtx, center = F, scale. = F)
fviz_eig(pca.res)
#ggplot(data.frame(pca.res$x), aes(x=PC1, y=PC2)) + geom_point()

df1 <- data.frame(pca.res$x)
df1$Cell.ID <- rownames(df1)
tr21.rna@meta.data$Cell.ID <- rownames(tr21.rna@meta.data)
df_meta <- tr21.rna@meta.data
df_meta <- df_meta[rownames(pca.res$x),]

df2 <- merge(df1, df_meta[,c("Age", 
                             "Cluster.Name", 
                             "Cell.ID", "Sample.ID")], by="Cell.ID")
rownames(df2) <- df2$Cell.ID

p1 <- ggplot(df2, aes(x=PC1, y=PC2)) +
  geom_point(aes(colour = Cluster.Name))
p2 <- ggplot(df2, aes(x=PC1, y=PC2)) +
  geom_point(aes(colour = Age))
p1+p2

fviz_pca_ind(pca.res,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = df2$Age, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups")

fviz_pca_var(pca.res, col.var = "steelblue")


fviz_pca_var(pca.res, select.var = list(contrib = 14), #n=14 at least three per main axis
             repel=TRUE, 
             title="")+theme_classic(base_size = 15)