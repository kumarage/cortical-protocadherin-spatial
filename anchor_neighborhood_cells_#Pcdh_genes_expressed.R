#-------------------------------------------------------------------------------------------------
# Get #Total_Pcdhg_genes data from the neighbouring cells in L2_3, L_4 and L5_6 cortical layers
#-------------------------------------------------------------------------------------------------

rm(list = ls())
gc()

setwd("D:/AWS_Instance/vizgen/Merscope_visualizer/6wkMbd1_combined/pipeline_1/Protocadherin_proximity/similar_expression_pattern/Pcdhg_along_layer_plots")

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
colnames(seurat_combined@meta.data)

# selected layers/slcies and conditions
main_cell_types <- c("ExN_L2_3","ExN_L4","ExN_L5_6")
datasets <- c("slice1","slice2","slice3")
conditions <- c("WT","KO")

# loop through layer / slice / condition
for (focus_cell_type in main_cell_types){
  for (dataset1 in datasets){
    for (condition1 in conditions){
      print('-------------------------------------------------------------------') 
      print('-------------------------------------------------------------------')
      print('-------------------------------------------------------------------')
      print('-------------------------------------------------------------------') 
      print(focus_cell_type)
      print('-------------------------------------------------------------------')
      print(dataset1)
      print('-------------------------------------------------------------------')
      print(condition1)
      print('-------------------------------------------------------------------')
      print('-------------------------------------------------------------------') 
      print('-------------------------------------------------------------------') 
      
      all_celltypes <- unique(seurat_combined@meta.data$celltype_v2)
      # keep target layer cells
      avoid_celltypes <- setdiff(main_cell_types,focus_cell_type)
      keep_celltypes <- setdiff(all_celltypes,avoid_celltypes)
      
      seurat_subset <- subset(seurat_combined, 
                              subset = celltype_v2 %in% keep_celltypes
                              & dataset == dataset1 
                              & condition == condition1
      )
      head(seurat_subset@meta.data)
      head(seurat_subset)
      
      #----------------------------------
      #Number of Pcdhg genes
      #----------------------------------
      cell_coords <- seurat_subset@meta.data[, c("center_x", "center_y")]
      cell_coords$center_y <- -cell_coords$center_y
      
      pcdh_vector <- c("Pcdhga1","Pcdhga2","Pcdhga3","Pcdhga4","Pcdhga5","Pcdhga6",
                        "Pcdhga7","Pcdhga8","Pcdhga9","Pcdhga10","Pcdhga11","Pcdhga12",
                        "Pcdhgb1","Pcdhgb2","Pcdhgb4","Pcdhgb5","Pcdhgb6","Pcdhgb7",
                        "Pcdhgb8","Pcdhgc3","Pcdhgc4","Pcdhgc5")
      
      # pcdh_vector <- c("Pcdha1","Pcdha2","Pcdha3","Pcdha4","Pcdha5","Pcdha6",
      #                   "Pcdha7","Pcdha8","Pcdha9","Pcdha11","Pcdha12","Pcdhac1",
      #                   "Pcdhac2")
      
      # pcdh_vector <- c("Pcdhb1","Pcdhb2","Pcdhb3","Pcdhb4","Pcdhb5","Pcdhb6",
      #                   "Pcdhb7","Pcdhb8","Pcdhb9","Pcdhb10","Pcdhb11","Pcdhb12",
      #                   "Pcdhb13","Pcdhb14","Pcdhb15","Pcdhb16","Pcdhb17","Pcdhb18",
      #                   "Pcdhb19","Pcdhb20","Pcdhb21","Pcdhb22")
      
      # binary expression (gene present/absent) and per-cell sum
      pcdh_expression <- seurat_subset@assays$RNA$data[pcdh_vector, ]
      # Convert the expression values to 1 (if > 0) or 0 (if = 0)
      pcdh_identity_matrix <- apply(pcdh_expression, 2, function(cell_expression) {
        as.integer(cell_expression > 0)  
      })
      pcdh_identity_matrix <- t(pcdh_identity_matrix)
      
      pcdh_identity_df <- data.frame(pcdh_identity_matrix)
      colnames(pcdh_identity_df) <- pcdh_vector
      head(pcdh_identity_df)
      
      pcdh_sums <- rowSums(pcdh_identity_df)
      seurat_subset <- AddMetaData(seurat_subset, metadata = pcdh_sums, col.name = "Total_pcdhg_exp")
      head(seurat_subset@meta.data)
      mean(seurat_subset@meta.data$Total_pcdhg_exp)
      
      #--------------------------------------
      #Go to a single celltype cortical layer
      cell_layer <- c(focus_cell_type)
      seurat_subset1 <- subset(seurat_subset, 
                               subset = celltype_v2 %in% focus_cell_type
      )
      print(unique(seurat_subset1$celltype_v2))
      cell_coords_cluster <- seurat_subset1@meta.data[, c("center_x", "center_y")]
      cell_coords_cluster$center_y <- -cell_coords_cluster$center_y
      #--------------------------------------
      #--------------------------------------
      # Load relevant anchor point coordinates
      fname1 <- paste0("break_coords_",dataset1,"_",condition1,"_",focus_cell_type,".csv")
      break_coords <- read.csv(fname1)
      #--------------------------------------
      #--------------------------------------
      #Calculate cumulative distance
      distances <- sqrt(diff(break_coords$x)^2+diff(break_coords$y)^2)
      cumulative_distance <- c(0,cumsum(distances))
      points_array <- cbind(break_coords$x, break_coords$y)
      colnames(points_array) <- c("x", "y")
      #--------------------------------------
      
      #--------------------------------------
      #--------------------------------------
      neighbor_types <- c("N","A")
      
      for (neighbor_type in neighbor_types) {
        print(neighbor_type)
        
        #--------------------------------------
        total_pcdh_exp <- NULL
        cell_coords_matrix <- NULL
        if (neighbor_type=="N"){
          #For neuronal cells RUN
          total_pcdh_exp <- seurat_subset1@meta.data$Total_pcdhg_exp
          cell_coords_matrix <- as.matrix(cell_coords_cluster)
        }else if (neighbor_type=="A"){
          #For all neighbour cells RUN
          total_pcdh_exp <- seurat_subset@meta.data$Total_pcdhg_exp
          cell_coords_matrix <- as.matrix(cell_coords)
        }
        
        break_coords_matrix <- as.matrix(points_array)
        break_coords$cumulative_distance <- cumulative_distance
        
        #--------------------------------------
        #--------------------------------------
        #calculate nearest neighbours
        #Here use object seurat_subset to check neighbourhood with all the cells and use
        #object seurat_subset1 to check neighbourhood with only the neuronal cells
        #--------------------------------------
        #Select a neighborhood of 30 cells
        #--------------------------------------
        library(FNN)
        k <- 20
        nn_result <- get.knnx(cell_coords_matrix, break_coords_matrix, k)
        
        neighbor_pcdh_exp <- NULL
        point_neighbors  <- NULL
        if (neighbor_type=="N"){
          #For neuronal cells RUN
          neighbor_pcdh_exp <- apply(nn_result$nn.index, 2, function(indices) {
            seurat_subset1@meta.data[rownames(seurat_subset1@meta.data)[indices], "Total_pcdhg_exp"]
          })
          point_neighbors <- apply(nn_result$nn.index, 2, function(indices) {
            rownames(seurat_subset1@meta.data)[indices]
          })
        }else if (neighbor_type=="A"){
          #For all neighbour cells RUN
          neighbor_pcdh_exp <- apply(nn_result$nn.index, 2, function(indices) {
            seurat_subset@meta.data[rownames(seurat_subset@meta.data)[indices], "Total_pcdhg_exp"]
          })
          point_neighbors <- apply(nn_result$nn.index, 2, function(indices) {
            rownames(seurat_subset@meta.data)[indices]
          })
        }
        
        head(point_neighbors)
        nrow <- dim(point_neighbors)[1]
        ncol <- dim(point_neighbors)[2]
        
        N <- nrow * ncol
        df <- data.frame(cell_Id = rep(NA_character_, N), Total_pcdh = rep(NA_real_, N),cell_Type = rep(NA_character_, N))
        
        for (i in seq_len(nrow)) {
          for (j in seq_len(ncol)) {
            cell_id <- point_neighbors[i, j]
            if (!is.na(cell_id)) {
              df[j + (i - 1) * ncol, 1] <- as.character(cell_id)
              df[j + (i - 1) * ncol, 2] <- seurat_subset@meta.data[cell_id, "Total_pcdhg_exp"]
              df[j + (i - 1) * ncol, 3] <- seurat_subset@meta.data[cell_id, "celltype_v2"]
            }
          }
        }
        str(df)
        
        if (neighbor_type=="N"){
          #For neuronal cells RUN, use Pcdhg for gamma, Pcdha for alpha and Pcdhb for beta
          fname2 <- paste0("Pcdhg_in_neighbors_",focus_cell_type,"_",dataset1,"_",condition1,"_ExN.csv")
          write.csv(df, fname2)
          
        }else if (neighbor_type=="A"){
          #For all neighbour cells RUN
          fname2 <- paste0("Pcdhg_in_neighbors_",focus_cell_type,"_",dataset1,"_",condition1,"_All.csv")
          write.csv(df, fname2)
        }
      }
      
      print('-------------------------------------------------------------------')
      print('--------------------------------DONE-------------------------------')
      print('-------------------------------------------------------------------') 
    }
  }
}
