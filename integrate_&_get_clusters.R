rm(list = ls())
gc()

setwd("D:/AWS_Instance/vizgen/Merscope_visualizer/6wkMbd1_combined/pipeline_1")

source("functions.R")

library(dplyr)
library(readr)
library(Matrix)
library(Seurat)
#library(Azimuth)
library(ggplot2)
library(patchwork)
library(SeuratWrappers)
options(future.globals.maxSize = 600 * 1024^2)
options(Seurat.object.assay.version = "v5")

###################################################################################################
# Load data
file_1 <- read_rds("KO_WT_slice_meta_added_1.rds")
file_2 <- read_rds("KO_WT_slice_meta_added_2.rds")
file_3 <- read_rds("KO_WT_slice_meta_added_3.rds")

###################################################################################################
# Custom normalization (defined in functions.R)
file_1 <- combined_normalization(file_1)
file_2 <- combined_normalization(file_2)
file_3 <- combined_normalization(file_3)

###########################################################################################
#Merge all 3 objects
seurat_combined <- merge(file_3, y = c(file_1,file_2), add.cell.ids = c("S3","S1","S2"))
DefaultAssay(seurat_combined) <- "RNA"

###########################################################################################
# scale, PCs and integration
seurat_combined <- FindVariableFeatures(seurat_combined, selection.method = "vst")
seurat_combined <- ScaleData(seurat_combined,vars.to.regress = c("volume"))
seurat_combined <- RunPCA(seurat_combined, features = VariableFeatures(object = seurat_combined))
ElbowPlot(seurat_combined,ndims = 50)

determinePCA(seurat_combined)

library(harmony)
seurat_combined <- RunHarmony(
  object = seurat_combined,
  group.by.vars = "dataset",
  dims = 1:34,
  plot_convergence = TRUE,
  theta=1,
  lambda=0.01,
  sigma=0.01,
  tau=0.9
)

###########################################################################################
# get clusters
seurat_combined <- FindNeighbors(seurat_combined, reduction = "harmony", dims = 1:34)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5,cluster.name = "harmony_clusters_v2")

###########################################################################################
# get umap embeddings and visualize clusters
seurat_combined <- RunUMAP(seurat_combined, 
                           dims = 1:10, 
                           reduction = "harmony",
                           n.threads=4)
DimPlot(seurat_combined, 
        reduction = "umap", 
        group.by = "harmony_clusters_v2",
        label = T,pt.size = 1.5)

saveRDS(seurat_combined,"6wk_processed_ALL_genes.rds")