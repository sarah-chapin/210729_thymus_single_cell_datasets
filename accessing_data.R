library(tidyverse)
library(dplyr)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(scater)
library(limma)
library(DropletUtils)
library(data.table)


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

write_csv(bornstein_2018_metadata, 
          file = "Bornstein_2018/GSE103967_metadata.csv")


bornstein_2018_sce <- 
  SingleCellExperiment(list(counts = bornstein_2018_data_matrix),
                       colData = bornstein_2018_metadata)
