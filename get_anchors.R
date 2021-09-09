library(tidyverse)
library(dplyr)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(scater)
library(limma)
library(DropletUtils)
library(data.table)
library(singleCellTK)
library(org.Mm.eg.db)
library(scMerge)
library(scran)

#Load seurat object
thymus_mouse_datasets_seurat <- 
  readRDS(file = "thymus_mouse_datasets_seurat_corrected.rds")

#Generate Seurat anchors to use in MAT2 integration
##Split into a list by study
mouse_thymus_list <- 
  SplitObject(thymus_mouse_datasets_seurat , split.by = "batch")

##Find variable features
mouse_thymus_list_with_features <- lapply(mouse_thymus_list, 
                                          FindVariableFeatures,
                                          selection.method = "vst",
                                          nfeatures = 2000,
                                          verbose = FALSE)

features <- 
  SelectIntegrationFeatures(object.list = mouse_thymus_list_with_features)

mouse_thymus_anchors <- 
  FindIntegrationAnchors(object.list = mouse_thymus_list_with_features, 
                         anchor.features = features)

saveRDS(mouse_thymus_anchors, "mouse_thymus_anchors_corrected.rds")

mouse_thymus_anchors_data <- Seurat::AnnotateAnchors(mouse_thymus_anchors)
#write.csv(mouse_thymus_anchors_data, "mouse_thymus_anchors_corrected_colnames.csv")

mouse_anchors <- mouse_thymus_anchors_data %>%
  mutate_all(~gsub("cell", "", .)) %>%
  subset(select = c("cell1", "cell2", "anchor.score")) %>%
  dplyr::rename(score = anchor.score) %>%
  mutate_if(is.character,as.numeric)

write.csv(mouse_anchors, "mouse_thymus_anchors_corrected_colnames.csv")
