library(harmony)
library(scran)
library(scater)

#Load data
mouse_thymus_datasets_sce <- 
  readRDS(file = "thymus_mouse_datasets_sce_corrected.rds")

thymus_mouse_datasets_seurat <- 
  readRDS(file = "thymus_mouse_datasets_seurat_corrected.rds")

#Calculate pca

