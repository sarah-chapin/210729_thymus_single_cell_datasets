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


#Dhalla, Fatima, et al. 
#"Biologically indeterminate yet ordered promiscuous gene expression in single 
#medullary thymic epithelial cells." The EMBO journal 39.1 (2020): e101828.

##Sample info
sample_info <- read.delim(file = "Dhalla_2020/E-MTAB-8105.sdrf.txt", 
                          sep = "\t") %>%
  as.data.frame()

sample_source_array <- sample_info %>%
  select(c("Derived.Array.Data.File", "Source.Name")) %>%
  mutate(Array.Name = removeExt(Derived.Array.Data.File)) %>%
  select(-c("Derived.Array.Data.File"))
  

  
dhalla_2020_h5_files <-  paste("Dhalla_2020/", 
                               list.files(path = "Dhalla_2020/", 
                                          pattern = "h5"),
                               sep = "")

dhalla_2020_h5_names <- removeExt(dhalla_2020_h5_files)
dhalla_2020_h5_names <- gsub("Dhalla_2020/", "", dhalla_2020_h5_names)

dhalla_2020_sces <- lapply(dhalla_2020_h5_files, read10xCounts)
names(dhalla_2020_sces) <- dhalla_2020_h5_names

add_source_name <- function(sce, array) {
  source <- sample_source_array[sample_source_array$Array.Name == array,]$Source.Name
  sce$Source <- source
  return(sce)
}

dhalla_2020_sces_with_sources <- mapply(add_source_name, 
                            sce = dhalla_2020_sces,
                            array = dhalla_2020_h5_names)

dhalla_2020_full_sce <- do.call(cbind, dhalla_2020_sces_with_sources)
dhalla_2020_full_sce$Study <- "Dhalla_2020"

dhalla_2020_full_symbols_sce <- convertGeneIDs(inSCE = dhalla_2020_full_sce,
                                               inSymbol = "ENSEMBL",
                                               outSymbol = "SYMBOL",
                                               database = "org.Mm.eg.db")

##Normalize the counts data
#Function to get logcounts-only sce

get_logcounts_sce <- function(sce) {
  #Convert counts to df
  counts_df <- as.data.frame(as.matrix(assay(sce)))
  sce_corrected <- SingleCellExperiment(list(counts = counts_df),
                                        colData = colData(sce))
  #Quick cluster
  set.seed(100)
  clust.sce_corrected <- 
    quickCluster(sce_corrected)
  
  #Obtain Size Factors
  sce_corrected_sf <- computeSumFactors(
    sce_corrected, 
    cluster=clust.sce_corrected, 
    min.mean=0.1)
  
  #Perform normalization
  sce_corrected_normalized <- 
    logNormCounts(sce_corrected_sf)
  
  #Extract logcounts
  logcounts_data <- 
    assay(sce_corrected_normalized, "logcounts")
  
  #Get SCE with only logcounts
  sce_logcounts_only <-
    SingleCellExperiment(list(logcounts = logcounts_data),
                         colData = colData(sce))
  
  return(sce_logcounts_only)
}

get_quick_logcounts_sce <- function(sce) {
  #Convert counts to df
  counts_df <- as.data.frame(as.matrix(assay(sce)))
  sce_corrected <- SingleCellExperiment(list(counts = counts_df),
                                        colData = colData(sce))
  sce_corrected_normalized <- 
    logNormCounts(sce_corrected)
  
  #Extract logcounts
  logcounts_data <- 
    assay(sce_corrected_normalized, "logcounts")
  
  #Get SCE with only logcounts
  sce_logcounts_only <-
    SingleCellExperiment(list(logcounts = logcounts_data),
                         colData = colData(sce))
  
  return(sce_logcounts_only)
}


dhalla_2020_full_symbols_sce_logcounts_only <- 
  get_logcounts_sce(dhalla_2020_full_symbols_sce)

dhalla_2020_full_symbols_sce_quick_logcounts_only <- 
  get_quick_logcounts_sce(dhalla_2020_full_symbols_sce)

saveRDS(dhalla_2020_full_symbols_sce_quick_logcounts_only, 
        file = "dhalla_2020_full_symbols_sce_quick_logcounts_only.rds")

#Brennecke, Philip, et al. 
#"Single-cell transcriptome analysis reveals coordinated ectopic 
#gene-expression patterns in medullary thymic epithelial cells." 
#Nature immunology 16.9 (2015): 933-941.

muc1Pos.counts <- read.delim("Brennecke_2015/E-MTAB-3346.processed.1/muc1Pos.counts",
                             sep = "\t")

PB10.counts <- read.delim("Brennecke_2015/E-MTAB-3624.processed.1/PB10.counts",
                             sep = "\t",
                          header = FALSE)

brennecke_sample_info <- read.delim(file = "Brennecke_2015/E-MTAB-3624.sdrf.txt", 
                          sep = "\t") %>%
  as.data.frame()

brennecke_sample_info_mouse <- 
  brennecke_sample_info[brennecke_sample_info$Characteristics.organism. == "Mus musculus", ]

##LOOKS LIKE THESE ASSAYS ARE ONLY FOR HUMAN DATA. THE MOUSE DATA HAS FASTQ 
#FILES THAT WOULD NEED TO BE ANALYZED BEFORE BEING ADDED AS AN SCE

#Wells, Kristen L., et al. 
#"Combined transient ablation and single-cell RNA-sequencing reveals the 
#development of medullary thymic epithelial cells." Elife 9 (2020): e60188.

#We are only interested in the control data for this paper
wells_2020_control <- read_csv("Wells_2020/GSE137699_combinedControl.csv") %>%
  as.data.frame() 

#Remove duplicated genes
wells_2020_dup_genes <- wells_2020_control$Gene[duplicated(wells_2020_control$Gene)]

wells_2020_control <- wells_2020_control %>%
  filter(!(wells_2020_control$Gene %in% wells_2020_dup_genes)) %>%
  column_to_rownames(var = "Gene")

wells_2020_control_barcodes<- sub("^[^_]*_", "", colnames(wells_2020_control))
wells_2020_control_source<- sub("(\\_.*)", "", colnames(wells_2020_control))

wells_2020_control_sce <- 
  SingleCellExperiment(list(logcounts =wells_2020_control))

wells_2020_control_sce$Sample <- "GSE137699_combinedControl.csv"
wells_2020_control_sce$Barcode <- wells_2020_control_barcodes
wells_2020_control_sce$Study <- "Wells_2020"
wells_2020_control_sce$Source <- wells_2020_control_source

saveRDS(wells_2020_control_sce, file = "wells_2020_control_sce.rds")
#Park, Jong-Eun, et al. 
#"A cell atlas of human thymic development defines T cell repertoire formation." 
#Science 367.6480 (2020).

##Convert mouse h5ad to h5seurat

Convert("Park_2020/HTA08.v02.A04.Science_mouse_total.new.h5ad", 
        dest = "h5seurat", 
        overwrite = FALSE)

park_2020_mouse <- 
  LoadH5Seurat("Park_2020/HTA08.v02.A04.Science_mouse_total.new.h5seurat")

park_2020_mouse_sce <- as.SingleCellExperiment(park_2020_mouse)

##Extract barcodes
park_2020_mouse_sce_FCAI<- 
  park_2020_mouse_sce[,grep("FCAI", colnames(park_2020_mouse_sce))]
park_2020_mouse_sce_wholeThy <- 
  park_2020_mouse_sce[, grep("wholeThy", colnames(park_2020_mouse_sce))]

park_2020_mouse_FCAI_barcodes <- 
  sub("^[^_]*-", "", colnames(park_2020_mouse_sce_FCAI))
park_2020_mouse_sce_FCAI$Barcode <- park_2020_mouse_FCAI_barcodes

park_2020_mouse_wholeThy_barcodes <- 
  sub("\\-.*", "", colnames(park_2020_mouse_sce_wholeThy))
park_2020_mouse_sce_wholeThy$Barcode <- park_2020_mouse_wholeThy_barcodes

park_2020_mouse_sce_with_barcodes <- cbind(park_2020_mouse_sce_FCAI,
                                           park_2020_mouse_sce_wholeThy)
park_2020_mouse_sce_with_barcodes$Study <- "Park_2020"


park_2020_mouse_logcounts_data <- 
  as.data.frame(as.matrix(assay(park_2020_mouse_sce_with_barcodes, "logcounts")))

park_2020_mouse_sce_with_barcodes_logcounts_only <-
  SingleCellExperiment(list(logcounts = park_2020_mouse_logcounts_data),
                       colData = colData(park_2020_mouse_sce_with_barcodes))

saveRDS(park_2020_mouse_sce_with_barcodes_logcounts_only, 
        file = "park_2020_mouse_sce_with_barcodes_logcounts_only.rds")

#Collect human data into sce
park_2020_human <- 
  LoadH5Seurat("Park_2020/HTA08.v01.A05.Science_human_epi.new.h5seurat")

park_2020_human_sce <- as.SingleCellExperiment(park_2020_human)

#Extract human barcodes
park_2020_human_barcodes <- gsub(".*-","", colnames(park_2020_human_sce))
park_2020_human_sce$Barcode <- park_2020_human_barcodes
park_2020_human_sce$Study <- "Park_2020"

#Tabula Muris Consortium. 
#"A single cell transcriptomic atlas characterizes aging tissues in the mouse." 
#Nature 583.7817 (2020): 590.

tabula_muris_2020_sce <- 
  read10xCounts("Tabula_Muris_2020/GSE132042_RAW/Thymus-10X_P7_11")

tabula_muris_2020_sce$Study <- "Tabula_Muris_2020"

tabula_muris_2020_sce_logcounts_only <- 
  get_logcounts_sce(tabula_muris_2020_sce)

tabula_muris_2020_sce_quick_logcounts_only <- 
  get_quick_logcounts_sce(tabula_muris_2020_sce)

saveRDS(tabula_muris_2020_sce_quick_logcounts_only,
        file = "tabula_muris_2020_sce_quick_logcounts_only.rds")

#Meredith, Matthew, et al. 
#"Aire controls gene expression in the thymic epithelium with ordered 
#stochasticity." Nature immunology 16.9 (2015): 942-949.

meredith_2015 = fread("Meredith_2015/GSE70798_SCS_MEC.csv.gz") %>%
  as.data.frame() %>%
  column_to_rownames(var = "V1")

meredith_2015_metadata <- read_delim(
  file = "Meredith_2015/GSE70798_SCS_MEC_metadata_table.txt",
  delim = "\t") %>%
  as.data.frame()

meredith_2015_sce <- 
  SingleCellExperiment(list(counts = meredith_2015),
                       colData = meredith_2015_metadata)

meredith_2015_sce$Study <- "Meredith_2015"

meredith_2015_sce_logcounts_only <- get_logcounts_sce(meredith_2015_sce)

meredith_2015_sce_quick_logcounts_only <- 
  get_quick_logcounts_sce(meredith_2015_sce)

saveRDS(meredith_2015_sce_quick_logcounts_only,
        file = "meredith_2015_sce_quick_logcounts_only.rds")

#Bornstein, Chamutal, et al. 
#"Single-cell mapping of the thymic stroma identifies IL-25-producing tuft 
#epithelial cells." Nature 559.7715 (2018): 622-626.

#Get single cell data with tuft cell annotations
bornstein_2018_metadata <- read.table(
  gzfile("Bornstein_2018/GSE103967_metadata.txt.gz"),
  header = TRUE,
  sep="\t") %>%
  column_to_rownames(var = "well")



bornstein_2018_data_files <- list.files(path = "Bornstein_2018/GSE103967_RAW/",
                              pattern = ".txt.gz", full.names = TRUE)

open_txt_gz <- function(file) {
  df <- read.table(
    gzfile(file),
    header = TRUE,
    sep="\t")
  return(df)
}

bornstein_2018_data_frames <- lapply(bornstein_2018_data_files, open_txt_gz)

bornstein_2018_data_matrix <- bind_cols(bornstein_2018_data_frames) 

#Remove columns from data matrix that don't have annotations
bornstein_2018_data_cols_to_remove <- unique(
  colnames(bornstein_2018_data_matrix)[! colnames(bornstein_2018_data_matrix) %in% 
                               rownames(bornstein_2018_metadata)])

anno.rows.to.remove <- unique(
  rownames(bornstein_2018_metadata)[! rownames(bornstein_2018_metadata) %in% 
                                 colnames(bornstein_2018_data_matrix)])

bornstein_2018_data_matrix <- bornstein_2018_data_matrix %>%
  dplyr::select(-all_of(bornstein_2018_data_cols_to_remove)) %>%
  as.matrix()

bornstein_2018_metadata <- bornstein_2018_metadata %>%
  rownames_to_column(var = "well") %>%
  filter(!well %in% anno.rows.to.remove) %>%
  column_to_rownames(var = "well") %>%
  plyr::rename(c("Cell_barcode" = "Barcode"))

bornstein_2018_metadata$Study <- "Bornstein_2018"

bornstein_cell_types <- 
  read_delim("Bornstein_2018/thymus_epithel_clusts.txt") %>%
  as.data.frame() %>%
  column_to_rownames(var = "...1") %>%
  rownames_to_column(var = "cell")

bornstein_2018_metadata <- bornstein_2018_metadata %>%
  rownames_to_column(var = "cell")

bornstein_2018_metadata_full <- merge(bornstein_2018_metadata,
                                      bornstein_cell_types,
                                      by = "cell",
                                      all = T) %>%
  column_to_rownames(var = "cell") %>%
  rename(cell.types = group)
  
write_csv(bornstein_2018_metadata_full, 
          file = "Bornstein_2018/GSE103967_metadata.csv")


bornstein_2018_sce <- 
  SingleCellExperiment(list(counts = bornstein_2018_data_matrix),
                       colData = bornstein_2018_metadata_full)

bornstein_2018_sce$Batch <- bornstein_2018_sce$Amp_batch_ID

bornstein_2018_sce_logcounts_only <- get_logcounts_sce(bornstein_2018_sce)

bornstein_2018_sce_quick_logcounts_only <- 
  get_quick_logcounts_sce(bornstein_2018_sce)

saveRDS(bornstein_2018_sce_quick_logcounts_only,
        file = "bornstein_2018_sce_quick_logcounts_only.rds")

#Select relevant coldata
select_colData_combine_sces <- function(sce) {
  if (!"cell.types" %in% names(colData(sce))) {
    sce$cell.types <- ""
  }
  #if (!"Batch" %in% names(colData(sce))) {
  #  sce$Batch <- ""
  #}
  colData(sce) <- subset(colData(sce), , c("Study", "cell.types")) 
  return(sce)
}

sces <- list(dhalla_2020_full_symbols_sce_quick_logcounts_only, 
             wells_2020_control_sce, 
             park_2020_mouse_sce_with_barcodes_logcounts_only,
             tabula_muris_2020_sce_quick_logcounts_only,
             meredith_2015_sce_quick_logcounts_only,
             bornstein_2018_sce_quick_logcounts_only)

sces_subset <- lapply(sces, select_colData_combine_sces) 

#Remove genes not common to all datasets
genes <- lapply(sces, rownames)
common_genes <- Reduce(intersect, genes)

remove_uncommon_genes <- function(sce, gene_list) {
  sce_common_gene_only <- subset(sce, rownames(sce) %in% gene_list, )
  return(sce_common_gene_only)
}

sces_subset_common_genes_only <- lapply(sces_subset, 
                                        remove_uncommon_genes,
                                        gene_list = common_genes)

mouse_datasets_sce <- do.call(cbind, sces_subset_common_genes_only)

#Convert to Seurat and then back to get rid of duplicate cell names
mouse_datasets_sce_logcounts <- as.matrix(assay(mouse_datasets_sce))

mouse_datasets_sce_as_matrix <- 
  SingleCellExperiment(list(logcounts = mouse_datasets_sce_logcounts),
                       colData = colData(mouse_datasets_sce))

mouse_datasets_seurat <- as.Seurat(mouse_datasets_sce_as_matrix,
                                   counts = NULL,
                                   data = "logcounts")

mouse_datasets_sce_cleaned <- as.SingleCellExperiment(mouse_datasets_seurat)

mouse_datasets_logcounts_df <- as.data.frame(assay(mouse_datasets_sce_cleaned))
mouse_datasets_metadata_df <- 
  as.data.frame(colData(mouse_datasets_sce_cleaned)) %>%
  subset(select = -c(orig.ident, ident))

saveRDS(mouse_datasets_sce_cleaned, file = "thymus_mouse_datasets_sce.rds")
saveRDS(mouse_datasets_seurat, file = "thymus_mouse_datasets_seurat.rds")
write.csv(mouse_datasets_logcounts_df, 
          file = "thymus_mouse_logcounts.csv", 
          row.names = TRUE)
write.csv(mouse_datasets_metadata_df, 
          file = "thymus_mouse_metadata.csv",
          row.names = TRUE)

#Generate Seurat anchors to use in MAT2 integration
##Split into a list by study
mouse_thymus_list <- SplitObject(mouse_datasets_seurat , split.by = "Study")

##Find variable features
mouse_thymus_list_with_features <- lapply(mouse_thymus_list, 
                                          FindVariableFeatures,
                                          selection.method = "vst",
                                          nfeatures = 2000,
                                          verbose = FALSE)

for (i in 1:length(mouse_thymus_list)) {
  mouse_thymus_list_with_features[[i]] <- 
    FindVariableFeatures(mouse_thymus_list[[i]], 
                         selection.method = "vst",
                         nfeatures = 2000, 
                         verbose = FALSE)
}

#Clean metadata
metadata_mouse <- read.csv(file = "thymus_mouse_metadata.csv") %>%
  column_to_rownames(var = "X") %>%
  subset(select = c("Study", "cell.types")) %>%
  dplyr::rename(batch = Study, type = cell.types)

write.csv(metadata_mouse, 
          file = "thymus_mouse_metadata_corrected.csv",
          row.names = TRUE)
  

write.csv(mouse_datasets_metadata_df, 
          file = "thymus_mouse_metadata.csv",
          row.names = TRUE)

#Clean SCE
mouse_datasets <- readRDS(file = "thymus_mouse_datasets_sce.rds")

cell_numbers <- 1:length(colnames(mouse_datasets))
cell_number_names <- paste("cell", cell_numbers)
colnames(mouse_datasets) <- cell_number_names

colData(mouse_datasets) <- colData(mouse_datasets) %>%
  subset(select = c("Study", "cell.types")) 

names(colData(mouse_datasets)) <- c("batch", "type")

mouse_datasets_seurat_corrected <- as.Seurat(mouse_datasets,
                                   counts = NULL,
                                   data = "logcounts")

mouse_datasets_logcounts_df_corrected <- 
  as.data.frame(assay(mouse_datasets))

mouse_datasets_metadata_df_corrected <- 
  as.data.frame(colData(mouse_datasets)) 


%>%
  subset(select = c("Study", "cell.types")) %>%
  dplyr::rename(batch = Study,
         type = cell.types)
  

saveRDS(mouse_datasets, file = "thymus_mouse_datasets_sce_corrected.rds")
saveRDS(mouse_datasets_seurat_corrected, 
        file = "thymus_mouse_datasets_seurat_corrected.rds")

write.csv(mouse_datasets_logcounts_df_corrected, 
          file = "thymus_mouse_logcounts_corrected.csv", 
          row.names = TRUE)
write.csv(mouse_datasets_metadata_df_corrected, 
          file = "thymus_mouse_metadata_corrected_colnames.csv",
          row.names = TRUE)
