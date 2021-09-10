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
library(cowplot)

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

pdf("seurat_integration_umap.pdf",
    width=20, height=10)
p1 + p2
dev.off()


#Convert to Single Cell Experiment Object
mouse.thymus.combined.sce <- as.SingleCellExperiment(mouse.thymus.combined)

#Make cluster proportions plot by batch
#Get cell proportions
# Function that returns the proportion of each cell cluster per batch
getCellProportions <- function(sce, celltype=FALSE){
  
  # Get the nb of cells
  total.cells <- table(colData(sce)$batch)
  
  # Gives number of cells in each batch per cluster
  cell_per_cluster <- lapply(levels(colData(sce)$seurat_clusters), function(l) {
    tmp <- sce[,colData(sce)$seurat_clusters == l]
    Dhalla_2020 <- tmp[,colData(tmp)$batch == "Dhalla_2020"]
    Bornstein_2018 <- tmp[,colData(tmp)$batch == "Bornstein_2018"]
    Meredith_2015 <- tmp[,colData(tmp)$batch == "Meredith_2015"]
    Park_2020 <- tmp[,colData(tmp)$batch == "Park_2020"]
    Tabula_Muris_2020 <- tmp[,colData(tmp)$batch == "Tabula_Muris_2020"]
    Wells_2020 <- tmp[,colData(tmp)$batch == "Wells_2020"]
    data.frame(Dhalla_2020=ncol(Dhalla_2020),
               Bornstein_2018=ncol(Bornstein_2018),
               Meredith_2015=ncol(Meredith_2015),
               Park_2020=ncol(Park_2020),
               Tabula_Muris_2020=ncol(Tabula_Muris_2020),
               Wells_2020=ncol(Wells_2020),
               seurat_clusters=l)
  }) %>%
    bind_rows
  
  # Get proportion of each cell cluster, out of the total number of WT or KO cells
  cell_per_cluster <- cell_per_cluster %>%
    mutate(Dhalla_2020_proportion = Dhalla_2020/total.cells['Dhalla_2020'],
           Bornstein_2018_proportion = Bornstein_2018/total.cells['Bornstein_2018'],
           Meredith_2015_proportion = Meredith_2015/total.cells['Meredith_2015'],
           Park_2020_proportion = Park_2020/total.cells['Park_2020'],
           Tabula_Muris_2020_proportion = Tabula_Muris_2020/total.cells['Tabula_Muris_2020'],
           Wells_2020_proportion = Wells_2020/total.cells['Wells_2020'])
  
  # Get long tibble
  cell_per_cluster <- cell_per_cluster %>%
    dplyr::select(-c(Bornstein_2018,
                     Dhalla_2020,
                     Meredith_2015,
                     Park_2020,
                     Tabula_Muris_2020,
                     Wells_2020)) %>%
    pivot_longer(ends_with('proportion'), names_to ="batch", values_to="proportion") %>%
    mutate(seurat_clusters=fct_inorder(seurat_clusters)) # keep same order for plotting
  
  
  # If we have cell annotation, we can add the cell type annotation
  if(celltype==TRUE){
    cell_per_cluster <- colData(sce) %>%
      as.data.frame %>%
      dplyr::select(seurat_clusters, type) %>% 
      distinct %>%
      right_join(cell_per_cluster, by="seurat_clusters") %>%
      arrange(as.numeric(seurat_clusters)) %>%
      mutate(seurat_clusters_cell=paste(seurat_clusters, type, sep="_"),
             seurat_clusters_cell=fct_inorder(seurat_clusters_cell)) # have a column with cluster number + cell type
  }
  
  return(cell_per_cluster)
}

cellPerClust<- getCellProportions(mouse.thymus.combined.sce) %>%
  as.data.frame()

p <- ggplot(cellPerClust) + 
  geom_bar(aes(x=seurat_clusters, y=proportion, fill=batch), stat="identity", position="dodge") +
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33')) +
  theme_cowplot() +
  theme(legend.position = "right",
        legend.title=element_blank(),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        plot.title = element_text(hjust = 0.5, size=20)) +
  ggtitle("Cell Proportions")


pdf("seurat_integration_cell_proportions.pdf",
    width=20, height=10)
p
dev.off()

