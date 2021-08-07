library(biomaRt)
library(AnnotationHub)
library(tidyverse)
library(dplyr)

#Get mouse protein coding genes
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

genemap_all_coding <- getBM(attributes=c('ensembl_transcript_id',
                                         "ensembl_gene_id",
                                         "external_gene_name"),
                        filters = c('biotype'),
                        values = list('protein_coding'),
                        mart = ensembl) %>%
  dplyr::rename(transcript_id = ensembl_transcript_id,
                gene_id = ensembl_gene_id,
                symbol = external_gene_name) %>%
  distinct()

write_csv(genemap_all_coding, file = "mouse_protein_coding_genes.csv")

#Get human protein coding genes
ensembl_human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

genemap_human_all_coding <- getBM(attributes=c('ensembl_transcript_id',
                                         "ensembl_gene_id",
                                         "external_gene_name"),
                            filters = c('biotype'),
                            values = list('protein_coding'),
                            mart = ensembl_human) %>%
  dplyr::rename(transcript_id = ensembl_transcript_id,
                gene_id = ensembl_gene_id,
                symbol = external_gene_name) %>%
  distinct()

write_csv(genemap_human_all_coding, file = "human_protein_coding_genes.csv")
