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

seurat_combined <- read_rds("6wk_processed_wth_celltypes.rds")
#----------------------------------------------------------
#DEG analysis by Excitatory neuron layers
celltype <- c("ExN_L2/3","ExN_L4","ExN_L5/6")

cell_subset <- subset(seurat_combined, 
                      subset = celltype_v2 %in% celltype
                      #& dataset == "slice1" 
)

cell_subset <- JoinLayers(cell_subset)
DefaultAssay(cell_subset)

de_results_condition <- FindMarkers(
  cell_subset,
  slot = "data",
  test.use = "wilcox",
  ident.1 = "KO", 
  ident.2 = "WT", 
  group.by = "condition",
  min.pct = 0,            
  logfc.threshold = 0 
)

head(de_results_condition)
write.csv(de_results_condition, "DEG/ExN_WT_vs_KO_wilcox.csv")

#----------------------------------------------------------
#DEG analysis by Inhibitory neuron layers
celltype <- c("InN_Pvalb+","InN_Sst+","InN_Vip+","InN_Lamp5+")

cell_subset <- subset(seurat_combined, 
                      subset = celltype_v2 %in% celltype
                      #& dataset == "slice1" 
)

cell_subset <- JoinLayers(cell_subset)
DefaultAssay(cell_subset)

de_results_condition <- FindMarkers(
  cell_subset,
  slot = "data",
  test.use = "wilcox",
  ident.1 = "KO", 
  ident.2 = "WT", 
  group.by = "condition",
  min.pct = 0,            
  logfc.threshold = 0 
)

head(de_results_condition)
write.csv(de_results_condition, "DEG/InN_WT_vs_KO_wilcox.csv")