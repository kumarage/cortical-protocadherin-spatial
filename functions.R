#determine PCA components

determinePCA = function(obj){
  std.PCA = c()
  scaled_data <- GetAssayData(obj, assay = "RNA", slot = "scale.data")
  for(i in 10){
    scaled_data_random = scaled_data
    scaled_data_random = matrix(sample(scaled_data_random), nrow = nrow(scaled_data_random))
    colnames(scaled_data_random) = colnames(scaled_data)
    rownames(scaled_data_random) = rownames(scaled_data)
    
    obj[["RNA"]]$scale.data = scaled_data_random
    obj_random = RunPCA(obj, npcs = 50, verbose = TRUE)
    std.PCA = c(std.PCA, obj_random@reductions$pca@stdev[1])
  }
  std.PCA = median(std.PCA)
  
  ElbowPlot(obj, ndims=50) +
    geom_hline(yintercept = std.PCA, color="red")
  
  num.PCA = sum(obj@reductions$pca@stdev >= std.PCA)
  message(paste0("The principle components to choose:", num.PCA))
  return(std.PCA)
}


remove_genes <- function(seurat_obj, genes_to_remove) {
  genes_to_keep <- setdiff(rownames(seurat_obj), genes_to_remove)
  seurat_obj <- subset(seurat_obj, features = genes_to_keep)
  return(seurat_obj)
}

#Do both normalizations at the same time
combined_normalization<- function(object) {
  counts <- GetAssayData(object, assay = "RNA", layer = "counts")
  norm_data = t(apply(counts, 1, function(x) x / object@meta.data$volume)) * mean(object@meta.data$volume)
  norm_data <- t(apply(norm_data, 1, function(x) x / colSums(norm_data) * 500))
  norm_data <- as(log(norm_data + 1), "dgCMatrix")
  object <- SetAssayData(object, assay = "RNA", layer = "data", new.data = norm_data)
  return(object)
}