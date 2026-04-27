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

###########################################################################################
# transfer clusters from DEG-removed object to full gene object
seurat_combined1 <- read_rds("6wk_processed_ALL_genes.rds")
seurat_combined2 <- read_rds("6wk_processed_DEG_removed.rds")

seurat_combined1$harmony_clusters_v2 <- seurat_combined2$harmony_clusters
saveRDS(seurat_combined1,"../6wk_processed.rds")

###########################################################################################
# get marker gene plots

seurat_combined <- read_rds("../6wk_processed.rds")

ExN_genes <- c("Slc17a7","Satb2","Rorb","Cux2")
InN_genes <- c("Gad1","Gad2","Sst","Vip","Pvalb","Lamp5")
Ast_genes <- c("Gfap","Aqp4","S100b","Aldh1l1")
MG_genes <- c("Aif1","Apbb1ip")
Oligo_genes <- c("Olig1","Olig2","Cspg4","Mbp","Sox10")
set1_genes <- c("Fezf2","Sncg","Neurod1","Cldn5","Scube1","Tbr1","Bcl11b")
set2_genes <- c("Chat","Map2","Slc1a2","Rbfox3","Adora2a","Ptpru","Lhx6","Drd1")
all_gene_list <- c("Slc17a7","Cux2","Neurod1","Satb2","Rorb",
                   "Fezf2","Bcl11b","Tbr1","Gad1","Gad2","Pvalb",
                   "Sst","Vip","Lamp5","S100b","Aldh1l1","Gfap",
                   "Olig1","Olig2","Sox10","Mbp","Cspg4","Aif1","Apbb1ip","Cldn5")

VlnPlot(seurat_combined, features = all_gene_list, group.by = "harmony_clusters_v2", 
        layer = "data", stack = TRUE, flip = TRUE)

###########################################################################################
# manual cluster → cell type annotation

seurat_combined$celltype_v2 <- case_when(
  seurat_combined$harmony_clusters_v2 %in% c("0","22","6","19","23") ~ "ExN_L2_3",
  seurat_combined$harmony_clusters_v2 %in% c("1") ~ "ExN_L4",
  seurat_combined$harmony_clusters_v2 %in% c("3","7","12","13","17","18","20") ~ "ExN_L5_6",
  
  seurat_combined$harmony_clusters_v2 %in% c("11") ~ "InN_Pvalb+",
  seurat_combined$harmony_clusters_v2 %in% c("14") ~ "InN_Lamp5+",
  seurat_combined$harmony_clusters_v2 %in% c("15") ~ "InN_Sst+",
  seurat_combined$harmony_clusters_v2 %in% c("21") ~ "InN_Vip+",
  
  seurat_combined$harmony_clusters_v2 %in% c("5","10") ~ "Oligo",
  seurat_combined$harmony_clusters_v2 %in% c("2","9") ~ "Astro",
  seurat_combined$harmony_clusters_v2 %in% c("4","16") ~ "Endo",
  seurat_combined$harmony_clusters_v2 %in% "8" ~ "MG",
  TRUE ~ "unknown"
)

umap_plot <- DimPlot(seurat_combined, 
                     reduction = "umap", 
                     group.by = "celltype_v2", 
                     label = T, 
                     pt.size = 1.0,
                     raster = F) +
  scale_color_manual(values = cell_colors) +  
  theme_minimal()

umap_plot
ggsave(
  filename = "umap_with_labels.pdf", 
  plot = umap_plot, 
  width = 10,    
  height = 8,     
  dpi = 600,     
  device = "pdf"  
)

all_gene_list <- c("Slc17a7","Cux2","Neurod1","Satb2","Rorb",
                   "Fezf2","Bcl11b","Tbr1","Gad1","Gad2","Pvalb",
                   "Sst","Vip","Lamp5","S100b","Aldh1l1","Gfap",
                   "Olig1","Olig2","Sox10","Mbp","Cspg4","Aif1","Apbb1ip","Cldn5","Adora2a")

cluster_list <- c("ExN_L2/3","ExN_L4","ExN_L5/6",
                  "InN_Pvalb+","InN_Sst+","InN_Vip+","InN_Lamp5+",
                  "Oligo","Astro","Endo","MG")

seurat_combined$celltype_v2 <- factor(seurat_combined$celltype_v2, levels = cluster_list)

vln_plot <- VlnPlot(seurat_combined, features = all_gene_list, group.by = "celltype_v2", 
                    layer = "data", stack = TRUE, flip = TRUE)

ggsave("vln_plot.pdf", 
       plot = vln_plot, 
       width = 8,  
       height = 10,
       dpi = 600)

saveRDS("6wk_processed_wth_celltypes.rds")