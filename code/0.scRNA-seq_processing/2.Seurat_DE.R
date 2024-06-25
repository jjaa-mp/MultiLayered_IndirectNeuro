# Libraries
library(Seurat)
library(readxl)
library(tibble)
library(dplyr)
library(biomaRt)
library(SeuratDisk)
library(ggplot2)

# Helpers
source("/home/rstudio/jm_rstudio/data_CBL/indNeuro_tmp/0.scRNA-seq_processing/helpers.R")

# Read raw data from rds file
mtx <- readRDS("/home/rstudio/jm_rstudio/data_indNeuro/rna_counts.rds")
dim(mtx)
# Read meta data
cell.meta <- read.table("/home/rstudio/jm_rstudio/data_indNeuro/rna_cell_metadata.txt", header = TRUE, sep='\t')

# Remove doublets using metadata
cell.meta <- cell.meta[cell.meta$DF_classification != "Doublet",]
cell.meta = cell.meta[cell.meta$Cell.ID %in% colnames(mtx),] 

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
Idents(tr21)

table(tr21$Cluster.Name)
table(tr21$Age)

columns.keep <- c("orig.ident", "Sample.ID", "Age", "Batch", "Cluster.Name", "percent.mt", "percent.zeros")
tr21@meta.data <- tr21@meta.data[,columns.keep]
tr21.rna

# QC
tr21 <- QC.Filtering(tr21, var.split = "Batch", 
                         by.mito = TRUE, 
                         by.count = TRUE, 
                         by.gene = TRUE, 
                         by.zeros = TRUE)


DefaultAssay(tr21) <- "RNA"

tr21 <- NormalizeData(
  tr21,
  normalization.method = "LogNormalize",
  scale.factor = 10000)

# DEG - Are there significant differences Late vs Early - stringent criteria
Idents(tr21) <- "Cluster.Name"
tr21.de.markers <- FindMarkers(tr21, ident.1 = "Late RG", ident.2 = "Early RG",
                               logfc.threshold = 1.25, min.pct = 0.25, min.diff.pct = 0.5)
head(tr21.de.markers)

VlnPlot(tr21, features = c('NR2F1', 'CLU'), idents = c("Early RG", "Late RG"), ncol = 2) & theme(axis.title.x = element_blank())

VlnPlot(tr21, features = c('HMGCR', 'KLF6', 'MVD', 'HMGCS1'), idents = c("Early RG", "Late RG"), ncol = 2) & theme(axis.title.x = element_blank())

VlnPlot(tr21, features = c('NR2F1', 'CLU', 'HMGCS1', 'KLF6', 'MVD', 'HMGCR'), idents = c("Early RG", "Late RG"), ncol = 6) & theme(axis.title.x = element_blank())
