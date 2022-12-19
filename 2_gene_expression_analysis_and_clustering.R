# load library ------------------------------------------------------------

library(MultiAssayExperiment)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)
library(msigdbr)
library(GSEABase)
library(edgeR)
library(singscore)
library(pheatmap)
library(NbClust)
library(miRBaseConverter)

working_dir = "/home/yusri/Documents/project/lusc_keratinization/"
setwd(working_dir)
genome <- "hg38"

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# load data -------------------------------------------------------------
rna_all <- loadRData("LUSC_rna_hg38.rda")

lusc_deg = DGEList(counts = assay(rna_all,'tpm_unstrand'), genes = rowData(rna_all))
prop_expressed = rowMeans(assay(rna_all,'tpm_unstrand') > 1)

# filter low expression gene
keep = prop_expressed > 0.5
rna = rna_all[keep,]
lusc_deg = lusc_deg[keep, , keep.lib.sizes = FALSE] 

assay(rna, 'logTPM') = log(assay(rna,'tpm_unstrand'))
all_genes = rowData(rna)$gene_name

# filter duplicate gene name
dups = unique(all_genes[duplicated(all_genes)])
rna = rna[!all_genes %in% dups, ]
rownames(rna) = rowData(rna)$gene_name

ensmbl_id_name_convert = rowData(rna)[,c("gene_id","gene_name")] %>% 
  as_tibble() %>%
  separate(col=gene_id,remove=FALSE,sep="\\.",into = c("id","id2")) %>%
  mutate(id2=NULL)

rm(dups,keep)

# get gene sets------------------------------------------------------------
keratinization_gene_set_name = c(
  "GOBP_KERATINIZATION",
  "GOBP_CORNIFICATION",
  "GOBP_KERATINOCYTE_APOPTOTIC_PROCESS",
  "GOBP_KERATINOCYTE_DEVELOPMENT",
  "GOBP_KERATINOCYTE_DIFFERENTIATION",
  "GOCC_KERATIN_FILAMENT",
  "GOBP_NEGATIVE_REGULATION_OF_KERATINOCYTE_DIFFERENTIATION",
  "GOBP_POSITIVE_REGULATION_OF_KERATINOCYTE_DIFFERENTIATION",
  "GOBP_REGULATION_OF_KERATINOCYTE_DIFFERENTIATION",
  "GOMF_KERATIN_FILAMENT_BINDING"
)
keratinization_gene_df = msigdbr(species = "Homo sapiens",category="C5") %>%
  dplyr::filter(gs_name %in% keratinization_gene_set_name) %>%
  dplyr::select(gs_name,gene_symbol) %>%
  dplyr::rename(set_name = gs_name,gene_name=gene_symbol) %>%
  distinct()


keratinization_gene_names = keratinization_gene_df$gene_name %>% unique
gene_set_name <- keratinization_gene_df$set_name %>% unique
gene_sets = list()
gsc = list()
for (i in gene_set_name) {
  print(i)
  genes = keratinization_gene_df %>%
    filter(set_name == i) %>%
    pull(gene_name)
  gsc[[i]] = GeneSet(genes,setName=i,geneIdType=SymbolIdentifier())
  gene_sets[[i]] = genes
  rm(genes)
  rm(i)
}
gsc <- GeneSetCollection(gsc)


# calculate Single sample scoring of molecular phenotypes
eranks = rankGenes(assay(rna, 'logTPM'))
lusc_multiscore =  multiScore(eranks, gsc)
lusc_keratinization_scores = lusc_multiscore$Scores

# Calculate best value of K. K=3 was selected from the result of NbClust
# res_nblust <- NbClust(as.tibble(t(lusc_keratinization_scores)),method ="complete")
# rna$cluster = as.factor(res_nblust$Best.partition)

# clustering the sample
distance_mat <- dist(t(lusc_keratinization_scores), method = 'euclidean')
Hierar_cl <- hclust(distance_mat, method = "complete")
fit <- cutree(Hierar_cl, k = 3 )
ordered_sample =Hierar_cl$labels[Hierar_cl$order]
rna$cluster = factor(fit,levels=c("2","1","3"),labels = c("high_keratinization","medium_keratinization","low_keratinization"))


identical(names(fit),rna$barcode)

sample_annotation_col = colData(rna)[c("sample_type","cluster")] %>% as.data.frame()
pheatmap(lusc_keratinization_scores,
         annotation_col = sample_annotation_col,
         main = "singscore",
         # scale= "row",
         cluster_rows = F,
         cluster_cols = T,
         cutree_cols=3,
         clustering_method = "complete",
         # treeheight_col=2,
         show_colnames = FALSE,
         fontsize_row=5)

# test------------------------------------------------------------------------
sample_annotation_col = colData(rna)[c("sample_type","cluster","paper_Expression.Subtype")] %>% 
  as.data.frame() %>%
  arrange(cluster) %>%
  drop_na()

pheatmap(lusc_keratinization_scores[,rownames(sample_annotation_col)],
         annotation_col = sample_annotation_col,
         main = "singscore",
         # scale= "row",
         cluster_rows = F,
         cluster_cols = F,
         cutree_cols=3,
         clustering_method = "complete",
         # treeheight_col=2,
         show_colnames = FALSE,
         fontsize_row=5)

# differential expression --------------------------------------------------------------

rna_wo_normal = rna[,rna$sample_type =="Primary Tumor"]
rownames(rna_wo_normal) = rowData(rna_wo_normal)$gene_id

# save transcriptomic profiling of tumor with cluster data if needed
# save(rna_wo_normal,file="rna_with_cluster.rda")

# TCGABiolink workflow of DEGs analysis
dataPrep <- TCGAanalyze_Preprocessing(
  object = rna_wo_normal, 
  cor.cut = 0.6
)                  

dataNorm <- TCGAanalyze_Normalization(
  tabDF = dataPrep,
  geneInfo = geneInfoHT,
  method = "gcContent"
)                

dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
)   

high_keratinization_sample = rna_wo_normal[,rna_wo_normal$cluster == "high_keratinization"]$barcode
low_keratinization_sample = rna_wo_normal[,rna_wo_normal$cluster == "low_keratinization"]$barcode

dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,high_keratinization_sample],
  mat2 = dataFilt[,low_keratinization_sample],
  Cond1type = "high_keratinization",
  Cond2type = "low_keratinization",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT",
  pipeline = "edgeR"
)

logdataFiltDE <- subset(dataFilt, subset = rownames(dataFilt) %in% rownames(dataDEGs)) %>% log1p()
logdataFiltDE_rownames= tibble(id = rownames(logdataFiltDE)) %>%
  dplyr::left_join(ensmbl_id_name_convert %>% dplyr::select(id,gene_name),by="id") %>%
  distinct()

if (identical(rownames(logdataFiltDE),logdataFiltDE_rownames$id)) {
  rownames(logdataFiltDE) =logdataFiltDE_rownames$gene_name
}

# save DEG analysis result and log(x+1) rna matrix of DEGs
write.table(logdataFiltDE,"lusc_keratinization_deg.tsv",col.names=NA)
write.csv(dataDEGs,"DEG_analysis_result.csv")
