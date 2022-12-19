# load library ------------------------------------------------------------
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)
library(ELMER)
library(ggpubr)

working_dir = "/home/yusri/Documents/project/lusc_keratinization/"
#working_dir = "/home/yusri/pCloudDrive/project/lusc_immune"

setwd(working_dir)
genome <- "hg38"

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# load data-----------------------------------------------------------------
met <- loadRData(("LUSC_met_hg38.rda"))
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))
met <- subset(met,subset = !as.character(seqnames(met)) %in% c("chrNA","chrX","chrY"))
rna <- loadRData("rna_with_cluster.rda")
# mean_methylation-----------------------------------------------------------
df1= as_tibble_col(met$sample,column_name =c("sample"))
df2=as_tibble(colData(rna)) %>% 
  dplyr::select(sample,cluster)
phenotype_data = plyr::join(
  df1,
  df2,
  by="sample",
  type="left")
met$cluster = phenotype_data$cluster
identical(phenotype_data$sample,met$sample)
rm(df1,df2,phenotype_data)
met_df <- data.frame(
  "Sample.mean" = colMeans(assay(met), na.rm = TRUE),
  "cluster" = met$cluster
) %>% 
  drop_na() %>%
  filter(cluster %in% c("high_keratinization","low_keratinization"))

ggpubr::ggboxplot(
  data = met_df,
  y = "Sample.mean",
  x = "cluster",
  color = "cluster",
  add = "jitter",
  ylab = expression(paste("Mean DNA methylation (", beta, "-values)")),
  xlab = ""
) + stat_compare_means() 
#sample_annotation_col = colData(met)[c("sample_type","cluster")] %>% as.data.frame()

# ELMER analysis-----------------------------------------

distal.probes <- get.feature.probe(
  genome = genome,
  met.platform = "450K",
  rm.chr = paste0("chr", c("X", "Y"))
)


mae <- createMAE(
  exp = rna,
  met = met,
  save = FALSE,
  linearize.exp = TRUE,
  filter.probes = distal.probes,
  met.platform = "450K",
  genome = genome,
  TCGA = TRUE
)

df1= as_tibble_col(mae$sample,column_name =c("sample"))
df2=as_tibble(colData(rna)) %>% 
  dplyr::select(sample,cluster)
phenotype_data = plyr::join(
  df1,
  df2,
  by="sample",
  type="left")
identical(phenotype_data$sample,mae$sample)

mae$cluster = phenotype_data$cluster
rm(df1,df2,phenotype_data)

sig.diff.hypo <- get.diff.meth(
  data = mae,
  group.col = "cluster",
  group1 =  "high_keratinization",
  group2 = "low_keratinization",
  save = TRUE,
  mode = "supervised",
  diff.dir = "hypo",
  cores = 8,
  dir.out = "result",
  pvalue = 0.01
)

nearGenes.hypo <- GetNearGenes(data = mae,
                               probes = sig.diff.hypo$probe,
                               numFlankingGenes = 20)

gene.pair.hypo <- get.pair(
  data = mae,
  group.col = "cluster",
  group1 =  "high_keratinization",
  group2 = "low_keratinization",
  nearGenes = nearGenes.hypo,
  mode = "supervised",
  diff.dir = "hypo",
  dir.out = "result",
  label = "hypo",
  cores = 8,
  permu.dir = "result/permu"
)

enriched.motif.hypo <- get.enriched.motif(
  data = mae,
  label = "hypo",
  probes = gene.pair.hypo$Probe,
  dir.out = "result",
  min.incidence=5
)

TF.hypo <- get.TFs(
  data = mae,
  group.col = "cluster",
  group1 =  "high_keratinization",
  group2 = "low_keratinization",
  mode = "supervised",
  enriched.motif = enriched.motif.hypo,
  dir.out = "result",
  diff.dir = "hypo",
  label = "hypo",
  cores = 8
)
TF.target.hypo <- getTFtargets(gene.pair.hypo,
                               enriched.motif.hypo,
                               TF.hypo,
                               dir.out = "result")

sig.diff.hyper <- get.diff.meth(
  data = mae,
  group.col = "cluster",
  group1 =  "high_keratinization",
  group2 = "low_keratinization",
  save = TRUE,
  mode = "supervised",
  diff.dir = "hyper",
  cores = 8,
  dir.out = "result",
  pvalue = 0.01
)

nearGenes.hyper <- GetNearGenes(data = mae,
                                probes = sig.diff.hyper$probe,
                                numFlankingGenes = 20)

gene.pair.hyper <- get.pair(
  data = mae,
  group.col = "cluster",
  group1 =  "high_keratinization",
  group2 = "low_keratinization",
  nearGenes = nearGenes.hyper,
  mode = "supervised",
  diff.dir = "hyper",
  dir.out = "result",
  label = "hypo",
  cores = 8,
  permu.dir = "result/permu"
)

enriched.motif.hyper <- get.enriched.motif(
  data = mae,
  label = "hyper",
  probes = gene.pair.hyper$Probe,
  dir.out = "result",
  min.incidence=5
)

TF.hyper <- get.TFs(
  data = mae,
  group.col = "cluster",
  group1 =  "high_keratinization",
  group2 = "low_keratinization",
  mode = "supervised",
  enriched.motif = enriched.motif.hyper,
  dir.out = "result",
  diff.dir = "hyper",
  label = "hyper",
  cores = 8
)
TF.target.hyper <- getTFtargets(gene.pair.hyper,
                                enriched.motif.hyper,
                                TF.hyper,
                                dir.out = "result")