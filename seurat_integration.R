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
library(patchwork)

#Load data
thymus_mouse_datasets_seurat <- 
  readRDS(file = "thymus_mouse_datasets_seurat_corrected.rds")

#Load anchors
mouse_thymus_anchors <- readRDS(file = "mouse_thymus_anchors_corrected.rds")

#Integrate data
# this command creates an 'integrated' data assay
mouse.thymus.combined <- IntegrateData(anchorset = mouse_thymus_anchors)

#Analyze integrated data
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(mouse.thymus.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
mouse.thymus.combined <- ScaleData(mouse.thymus.combined, verbose = FALSE)
mouse.thymus.combined <- RunPCA(mouse.thymus.combined, npcs = 30, verbose = FALSE)
mouse.thymus.combined <- RunUMAP(mouse.thymus.combined, reduction = "pca", dims = 1:30)
mouse.thymus.combined <- FindNeighbors(mouse.thymus.combined, reduction = "pca", dims = 1:30)
mouse.thymus.combined <- FindClusters(mouse.thymus.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(mouse.thymus.combined, reduction = "umap", group.by = "batch")
p2 <- DimPlot(mouse.thymus.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

#Convert to Single Cell Experiment Object
mouse.thymus.combined.sce <- as.SingleCellExperiment(mouse.thymus.combined)

#Make cluster proportions plot by batch
#Get cell proportions
# Function that returns the proportion of each cell cluster per batch
getCellProportions <- function(sce, celltype=FALSE){
  
  # Get the nb of control, AA, and PGE cells
  total.cells <- table(colData(sce)$batch)
  
  # Gives number of control, AA, and PGE cells per cluster
  cell_per_cluster <- lapply(levels(colData(sce)$seurat_clusters), function(l) {
    tmp <- sce[,colData(sce)$seurat_clusters == l]
    control <- tmp[,colData(tmp)$Condition == "control"]
    AA <- tmp[,colData(tmp)$Condition == "AA"]
    PGE <- tmp[,colData(tmp)$Condition == "PGE"]
    data.frame(control=ncol(control),
               AA=ncol(AA),
               PGE=ncol(PGE),
               label=l)
  }) %>%
    bind_rows
  
  # Get proportion of each cell cluster, out of the total number of WT or KO cells
  cell_per_cluster <- cell_per_cluster %>%
    mutate(control_proportion = control/total.cells['control'],
           AA_proportion = AA/total.cells['AA'],
           PGE_proportion = PGE/total.cells['PGE'])
  
  # Get long tibble
  cell_per_cluster <- cell_per_cluster %>%
    dplyr::select(-c(control, AA, PGE)) %>%
    pivot_longer(ends_with('proportion'), names_to ="condition", values_to="proportion") %>%
    mutate(label=fct_inorder(label)) # keep same order for plotting
  
  
  # If we have cell annotation, we can add the cell type annotation
  if(celltype==TRUE){
    cell_per_cluster <- colData(sce) %>%
      as.data.frame %>%
      dplyr::select(label, cell_type) %>% 
      distinct %>%
      right_join(cell_per_cluster, by="label") %>%
      arrange(as.numeric(label)) %>%
      mutate(label_cell=paste(label, cell_type, sep="_"),
             label_cell=fct_inorder(label_cell)) # have a column with cluster number + cell type
  }
  
  return(cell_per_cluster)
}

cellPerClust<- getCellProportions(cluster.walktrap.full.40$sce) %>%
  as.data.frame()

p <- ggplot(cellPerClust) + 
  geom_bar(aes(x=label, y=proportion, fill=condition), stat="identity", position="dodge") +
  scale_fill_manual(values=c("#abd9e9", "#fdae61", "#91cf60")) +
  theme_cowplot() +
  theme(legend.position = "right",
        legend.title=element_blank(),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        plot.title = element_text(hjust = 0.5, size=20)) +
  ggtitle("Cell Proportions")
p

