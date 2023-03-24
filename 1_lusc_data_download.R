# load library ------------------------------------------------------------

library(MultiAssayExperiment)
library(ELMER)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(EDASeq)
library(edgeR)
library(tidyverse)

working_dir = "/home/yusri/Documents/project/lusc_keratinization/"

setwd(working_dir)
genome <- "hg38"

# download transcription profiling data-----------------------------------------

query.rna <- GDCquery(project = "TCGA-LUSC",
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          workflow.type = "STAR - Counts",
                          sample.type = c("Primary Tumor","Solid Tissue Normal"),
                      # sample.type = c("Primary Tumor"),
                           legacy=FALSE)

# GDCdownload(query.rna)


rna <-GDCprepare(query.rna, 
                save = TRUE,
                save.filename=LUSC_rna_hg38.rda)


# download methylation data ------------------------------------------------------------
query.met <- GDCquery(project = "TCGA-LUSC",
                      data.category = "DNA Methylation",
                      data.type="Methylation Beta Value",
                      platform = "Illumina Human Methylation 450",
                      sample.type = c("Primary Tumor","Solid Tissue Normal"),
                      legacy=FALSE)
# GDCdownload(query.met)
met <-GDCprepare(query.met, 
                 save = TRUE,
                 save.filename = "LUSC_met_hg38.rda" 
)


# download mutation data --------------------------------------------------

query_mut <- GDCquery(
  project = "TCGA-LUSC", 
  data.category = "Simple Nucleotide Variation", 
  legacy = FALSE, 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
# GDCdownload(query_mut)
maf_df <- GDCprepare(query_mut,save = TRUE,
                  save.filename = "LUSC_mut_hg38.rda" )

# download mirna data

query_miRNA <- GDCquery(
  project = "TCGA-LUSC", 
  experimental.strategy = "miRNA-Seq",
  data.category = "Transcriptome Profiling", 
  data.type = "miRNA Expression Quantification",
  legacy = FALSE,
  sample.type = c("Primary Tumor")
)

GDCdownload(query = query_miRNA)

mirna <- GDCprepare(
  query = query_miRNA,save = TRUE,
  save.filename = "LUSC_mirna_hg38.rda"
)

# download mirna isoform data

query_isoform <- GDCquery(
  project = "TCGA-LUSC", 
  data.category = "Transcriptome Profiling", 
  data.type = "Isoform Expression Quantification",
  legacy = FALSE,
  sample.type = c("Primary Tumor")
)

GDCdownload(query = query_isoform)

isoform <- GDCprepare(
  query = query_isoform,save = TRUE,
  save.filename = "LUSC_isoform_hg38.rda"
)
