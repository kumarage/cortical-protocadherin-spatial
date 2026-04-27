#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#Calculate cell-wise Jaccard Index in anchor(break coord) neighborhoods
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#For each layer calculate Jaccard Similarity Score 
#from the neighbouring cells in a given celltype and given neighbourhood per each cell
#-------------------------------------------------------------------------------------------------
rm(list = ls())
gc()

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

#Cell proximity study selection
main_cell_types <- c("ExN_L2_3","ExN_L4","ExN_L5_6")
datasets <- c("slice1","slice2","slice3")
conditions <- c("WT","KO")

#-------------------------------------------------------------------------------
# cellwise Jaccard Index function for binary gene vectors
#-------------------------------------------------------------------------------
calculate_JI_Index <- function(root, query, pcdhg_df) {
  
  row1 <- as.numeric(pcdhg_df[root, ])
  row2 <- as.numeric(pcdhg_df[query, ])
  
  dot_product <- sum(row1 * row2)
  pcd_union <- sum((row1 != 0)| (row2 != 0))
  
  if(pcd_union != 0){
    JI <- dot_product/pcd_union
  }else{
    JI <- 0
  }
  
  return(JI)
}

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
      #Number of Pcdh genes
      #----------------------------------
      cell_coords <- seurat_subset@meta.data[, c("center_x", "center_y")]
      cell_coords$center_y <- -cell_coords$center_y
      
      pcdh_vector <- c("Pcdhga1","Pcdhga2","Pcdhga3","Pcdhga4","Pcdhga5","Pcdhga6",
                        "Pcdhga7","Pcdhga8","Pcdhga9","Pcdhga10","Pcdhga11","Pcdhga12",
                        "Pcdhgb1","Pcdhgb2","Pcdhgb4","Pcdhgb5","Pcdhgb6","Pcdhgb7",
                        "Pcdhgb8","Pcdhgc3","Pcdhgc4","Pcdhgc5")
      
      #pcdh_vector <- c("Pcdha1","Pcdha2","Pcdha3","Pcdha4","Pcdha5","Pcdha6",
      #                   "Pcdha7","Pcdha8","Pcdha9","Pcdha11","Pcdha12","Pcdhac1",
      #                   "Pcdhac2")
      
      #pcdh_vector <- c("Pcdhb1","Pcdhb2","Pcdhb3","Pcdhb4","Pcdhb5","Pcdhb6",
      #                  "Pcdhb7","Pcdhb8","Pcdhb9","Pcdhb10","Pcdhb11","Pcdhb12",
      #                  "Pcdhb13","Pcdhb14","Pcdhb15","Pcdhb16","Pcdhb17","Pcdhb18",
      #                  "Pcdhb19","Pcdhb20","Pcdhb21","Pcdhb22")
      
      pcdh_expression <- seurat_subset@assays$RNA$data[pcdh_vector, ]
      # Convert the expression values to 1 (if > 0) or 0 (if = 0)
      pcdh_identity_matrix <- apply(pcdh_expression, 2, function(cell_expression) {
        as.integer(cell_expression > 0)  
      })
      pcdh_identity_matrix <- t(pcdh_identity_matrix)
      
      pcdh_identity_df <- data.frame(pcdh_identity_matrix)
      colnames(pcdh_identity_df) <- pcdh_vector
      head(pcdh_identity_df)
      #--------------------------------------
      #--------------------------------------
      # restrict to target single layer
      cell_layer <- c(focus_cell_type)
      seurat_subset1 <- subset(seurat_subset, 
                               subset = celltype_v2 %in% focus_cell_type
      )
      print(unique(seurat_subset1$celltype_v2))
      cell_coords_cluster <- seurat_subset1@meta.data[, c("center_x", "center_y")]
      cell_coords_cluster$center_y <- -cell_coords_cluster$center_y
      pcdh_identity_df_cluster <- pcdh_identity_df[rownames(seurat_subset1@meta.data),]
      #--------------------------------------
      #--------------------------------------
      # Load data dots coordinates
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
        pcdh_identity <- NULL
        cell_coords_matrix <- NULL
        if (neighbor_type=="N"){
          #For neuronal cells RUN
          pcdh_identity <- pcdh_identity_df_cluster
          cell_coords_matrix <- as.matrix(cell_coords_cluster)
        }else if (neighbor_type=="A"){
          #For all neighbour cells RUN
          pcdh_identity <- pcdh_identity_df
          cell_coords_matrix <- as.matrix(cell_coords)
        }
        
        break_coords_matrix <- as.matrix(points_array)
        break_coords$cumulative_distance <- cumulative_distance
        
        library(FNN)
        # global kNN graph for all cells
        kk <- 30
        knn_results_all <- get.knn(cell_coords_matrix, k = kk)
        
        # anchor-based neighborhoods
        k <- 20
        knn_results <- get.knnx(cell_coords_matrix, break_coords_matrix, k)
        neighborhood_count <- nrow(knn_results$nn.index)
        
        output_list <- list()
        
        for (n_idx in 1:neighborhood_count) {
          anchor_indices <- knn_results$nn.index[n_idx, ]
          anchor_cell_ids <- rownames(pcdh_identity)[anchor_indices]
          
          for (cell_idx in anchor_indices) {
            this_cell_id <- rownames(pcdh_identity)[cell_idx]
            neighbor_indices <- knn_results_all$nn.index[cell_idx, ]
            neighbor_cell_ids <- rownames(pcdh_identity)[neighbor_indices]
            
            # Get similarity of this cell to each of its 30 neighbors
            this_vector <- as.numeric(pcdh_identity[this_cell_id, ])
            
            similarity_vector <- sapply(neighbor_cell_ids, function(n_id) {
              calculate_JI_Index(this_cell_id, n_id, pcdh_identity)
            })
            
            similarity_named <- setNames(as.list(similarity_vector), paste0("JI_", seq_along(similarity_vector)))
            
            output_row <- data.frame(
              neighborhood = n_idx,
              cell_id = this_cell_id,
              similarity_named,
              stringsAsFactors = FALSE
            )
            
            output_list[[length(output_list) + 1]] <- output_row
          }
        }
        
        output_df <- do.call(rbind, output_list)
        colnames(output_df)[3:ncol(output_df)] <- paste0("JI_", 1:(ncol(output_df) - 2))
        
        if (neighbor_type=="N"){
          #For neuronal cells RUN
          fname2 <- paste0("gamma_Pcdh_Jaccard_score_of_neighbors_",focus_cell_type,"_",dataset1,"_",condition1,"_ExN.csv")
          write.csv(output_df, fname2)
          
        }else if (neighbor_type=="A"){
          #For all neighbor cells RUN
          fname2 <- paste0("gamma_Pcdh_Jaccard_score_of_neighbors_",focus_cell_type,"_",dataset1,"_",condition1,"_All.csv")
          write.csv(output_df, fname2)
        }
      }
      print('-------------------------------------------------------------------')
      print('--------------------------------DONE-------------------------------')
      print('-------------------------------------------------------------------') 
    }
  }
}
