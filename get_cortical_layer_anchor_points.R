#-------------------------------------------------------------------------------------------------
# manually draw cortical layer path and generate evenly spaced anchor points
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

main_cell_types <- c("ExN_L2_3","ExN_L4","ExN_L5_6")
#datasets <- c("slice1","slice2","slice3")
#conditions <- c("WT","KO")

focus_cell_type <- "ExN_L2_3"
dataset1 <- "slice1"
condition1 <- "WT"

# keep target layer + other cell types
all_celltypes <- unique(seurat_combined@meta.data$celltype_v2)
avoid_celltypes <- setdiff(main_cell_types,focus_cell_type)
keep_celltypes <- setdiff(all_celltypes,avoid_celltypes)
seurat_subset <- subset(seurat_combined, 
                        subset = celltype_v2 %in% keep_celltypes
                        & dataset == dataset1 
                        & condition == condition1
)
head(seurat_subset@meta.data)
head(seurat_subset)

#--------------------------------------
#--------------------------------------
# isolate coordinates for selected layer
cell_layer <- c(focus_cell_type)
seurat_subset1 <- subset(seurat_subset, 
                         subset = celltype_v2 %in% cell_layer
)
cell_coords_cluster <- seurat_subset1@meta.data[, c("center_x", "center_y")]
#--------------------------------------
#--------------------------------------
# plot spatial coordinates
plot(cell_coords_cluster$center_x, -cell_coords_cluster$center_y, 
     main = paste0(dataset1," ",condition1," ",focus_cell_type," points"), 
     xlab = "center_x", ylab = "center_y", 
     pch = 20, col = "blue",
     asp = 1)

# manually click points along cortical layer
clicked_points <- locator(type = "p", col = "red", pch = 4)
break_coords <- clicked_points
points(break_coords$x, break_coords$y, col = "green", pch = 18, cex = 1.0)

#--------------------------------------
#--------------------------------------
# cumulative distance along the path
distances <- sqrt(diff(break_coords$x)^2+diff(break_coords$y)^2)
cumulative_distance <- c(0,cumsum(distances))
points_array <- cbind(break_coords$x, break_coords$y)
colnames(points_array) <- c("x", "y")
points_array

#seg_len <- 100
#n <- ceiling(cumulative_distance[length(cumulative_distance)]/seg_len)
n <- 50
#--------------------------------------
# Calculate total cumulative distance and create equally spaced points along curve
total_distance <- cumulative_distance[length(cumulative_distance)]
target_distances <- seq(0, total_distance, length.out = n)

# Find interpolation indices using cumulative distances
indices <- findInterval(target_distances, cumulative_distance, 
                        rightmost.closed = TRUE)
indices <- pmin(indices, length(cumulative_distance) - 1)  # Handle edge case

# Linear interpolation between original points
x_new <- numeric(n)
y_new <- numeric(n)

for (j in seq_along(target_distances)) {
  i <- indices[j]
  seg_start <- points_array[i, ]
  seg_end <- points_array[i + 1, ]
  
  segment_length <- cumulative_distance[i + 1] - cumulative_distance[i]
  if (segment_length == 0) {
    x_new[j] <- seg_start["x"]
    y_new[j] <- seg_start["y"]
  } else {
    fraction <- (target_distances[j] - cumulative_distance[i]) / segment_length
    x_new[j] <- seg_start["x"] + fraction * (seg_end["x"] - seg_start["x"])
    y_new[j] <- seg_start["y"] + fraction * (seg_end["y"] - seg_start["y"])
  }
}

equal_points <- cbind(x = x_new, y = y_new)
# save interpolated points
mod_points <- list(
  x = equal_points[, "x"],
  y = equal_points[, "y"] 
)
points(mod_points$x, mod_points$y, col = "yellow", pch = 18, cex = 1.0)
# save anchor points
fname1 <- paste0("break_coords_",dataset1,"_",condition1,"_",focus_cell_type,".csv")
write.csv(mod_points,fname1)
#-------------------------------------------------------------------------------------------------
