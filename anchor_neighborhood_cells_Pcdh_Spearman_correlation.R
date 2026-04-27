#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#Calculate Spearman Correlation
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#Calculate Spearman Correlation from the neighbouring cells in L2_3, L_4 and L5_6 cortical layers
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
options(future.globals.maxSize = 10 * 1024^3)
options(Seurat.object.assay.version = "v5")

seurat_combined <- read_rds("6wk_processed_wth_celltypes.rds")
colnames(seurat_combined@meta.data)

#Cell proximity study selection
main_cell_types <- c("ExN_L2_3","ExN_L4","ExN_L5_6")
datasets <- c("slice1","slice2","slice3")
conditions <- c("WT","KO")

#-------------------------------------------------------------------------------
#cellwise Spearman Correlation function on normalized expression
#-------------------------------------------------------------------------------
calculate_Spearman_corr_alt <- function(root, query, pcdhg_df) {
  
  #Gets Pcdhg expression values of center cell
  row1 <- as.numeric(pcdhg_df[root, ])
  row2 <- as.numeric(pcdhg_df[query, ])
  
  SP_corr <- cor(row1, row2, method = "spearman", use = "complete.obs")
  
  return(SP_corr)
}
#-------------------------------------------------------------------------------
# average pairwise Spearman correlation within a neighborhood
#-------------------------------------------------------------------------------

calculate_SC_for_neighborhood <- function(indices, pcdhg_df) {
  
  cluster_index <- indices
  cluster_cells <- rownames(pcdhg_df)[cluster_index]
  total_cells <- length(cluster_cells)
  
  cum_SC <- 0
  counter <- 0
  for (cell1 in cluster_cells){
    for (cell2 in cluster_cells){
      if(cell1 != cell2){
        SC <- calculate_Spearman_corr_alt(cell1,cell2,pcdhg_df)
        if(!is.na(SC)){
          counter <- counter + 1
          cum_SC <- cum_SC + SC
        }
      }
    }
  }
  
  # Compute Average Spearman correlation for the neighborhood
  similarity_score <- cum_SC/counter
  return(similarity_score)
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
      #Avoid other layer neurons
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
      
      pcdha_vector <- c("Pcdha1","Pcdha2","Pcdha3","Pcdha4","Pcdha5","Pcdha6",
                        "Pcdha7","Pcdha8","Pcdha9","Pcdha11","Pcdha12","Pcdhac1",
                        "Pcdhac2")
      
      pcdhb_vector <- c("Pcdhb1","Pcdhb2","Pcdhb3","Pcdhb4","Pcdhb5","Pcdhb6",
                        "Pcdhb7","Pcdhb8","Pcdhb9","Pcdhb10","Pcdhb11","Pcdhb12",
                        "Pcdhb13","Pcdhb14","Pcdhb15","Pcdhb16","Pcdhb17","Pcdhb18",
                        "Pcdhb19","Pcdhb20","Pcdhb21","Pcdhb22")
      
      pcdhg_vector <- c("Pcdhga1","Pcdhga2","Pcdhga3","Pcdhga4","Pcdhga5","Pcdhga6",
                        "Pcdhga7","Pcdhga8","Pcdhga9","Pcdhga10","Pcdhga11","Pcdhga12",
                        "Pcdhgb1","Pcdhgb2","Pcdhgb4","Pcdhgb5","Pcdhgb6","Pcdhgb7",
                        "Pcdhgb8","Pcdhgc3","Pcdhgc4","Pcdhgc5")
      
      pcdh_all_vector <- c(pcdha_vector,pcdhb_vector,pcdhg_vector)
      pcdh_expression <- seurat_subset@assays$RNA$data[pcdh_all_vector, ]
      pcdh_matrix <- t(pcdh_expression)
      pcdh_df <- data.frame(pcdh_matrix)
      colnames(pcdh_df) <- pcdh_all_vector
      head(pcdh_df)
      #--------------------------------------
      #--------------------------------------
      #Go to a single celltype cortical layer
      cell_layer <- c(focus_cell_type)
      seurat_subset1 <- subset(seurat_subset, 
                               subset = celltype_v2 %in% focus_cell_type
      )
      print(unique(seurat_subset1$celltype_v2))
      cell_coords_cluster <- seurat_subset1@meta.data[, c("center_x", "center_y")]
      cell_coords_cluster$center_y <- -cell_coords_cluster$center_y
      pcdh_df_cluster <- pcdh_df[rownames(seurat_subset1@meta.data),]
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
        pcdhg_identity <- NULL
        cell_coords_matrix <- NULL
        if (neighbor_type=="N"){
          #For neuronal cells RUN
          pcdhg_identity <- pcdh_df_cluster
          cell_coords_matrix <- as.matrix(cell_coords_cluster)
        }else if (neighbor_type=="A"){
          #For all neighbour cells RUN
          pcdhg_identity <- pcdh_df
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
        #Select a neighborhood of 20 cells
        #--------------------------------------
        
        
        library(FNN)
        k <- 20
        knn_results <- get.knnx(cell_coords_matrix, break_coords_matrix, k)
        Spearman_scores <- apply(knn_results$nn.index, 1, function(indices) {
          #Call the function
          calculate_SC_for_neighborhood(indices, pcdhg_identity)
        })
        
        output_df <- NULL
        output_df <- data.frame(
          Spearman_scores = Spearman_scores,
          cumulative_distance = break_coords$cumulative_distance
        )
        str(output_df)
        
        if (neighbor_type=="N"){
          #For neuronal cells RUN
          fname2 <- paste0("Pcdh_Spearman_score_",focus_cell_type,"_",dataset1,"_",condition1,"_ExN.csv")
          write.csv(output_df, fname2)
          
        }else if (neighbor_type=="A"){
          #For all neighbor cells RUN
          fname2 <- paste0("Pcdh_Spearman_score_",focus_cell_type,"_",dataset1,"_",condition1,"_All.csv")
          write.csv(output_df, fname2)
        }
      }
      print('-------------------------------------------------------------------')
      print('--------------------------------DONE-------------------------------')
      print('-------------------------------------------------------------------') 
    }
  }
}