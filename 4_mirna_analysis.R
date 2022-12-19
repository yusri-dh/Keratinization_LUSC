# load library ------------------------------------------------------------

library(MultiAssayExperiment)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)
library(edgeR)
library(mirTarRnaSeq)
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
rna <- loadRData("rna_with_cluster.rda")
dataDEGs <- read.csv("DEG_analysis_result.csv",row.names=1)
  

gene_id_name_pairs = dataDEGs %>%
  rownames_to_column("gene_id") %>%
  dplyr::select(c("gene_id","gene_name"))

# Only take mature miRNA
isoform <-loadRData("LUSC_isoform_hg38.rda") %>%
  dplyr::select(read_count,miRNA_region,barcode) %>%
  mutate(sample=substr(barcode,1,16))  %>%
  filter(grepl("mature",miRNA_region)) %>%
  mutate(acc_id=substring(miRNA_region,8)) %>%
  mutate(mirna_id = miRNA_AccessionToName(acc_id)$TargetName)
sum(grepl("^MIMAT",isoform$acc_id)) == dim(isoform)[1]

# create mirna expression matrix
mirna= aggregate(read_count~mirna_id+sample,data=isoform,FUN=sum) %>%
  spread(sample,read_count) %>%
  column_to_rownames("mirna_id")

# differential mirna analysis-----------------------------------------------------------------
df1 = tibble(sample=colnames(mirna))
df2=as_tibble(colData(rna)) %>% 
  dplyr::select(sample,cluster,sample_type,barcode)
phenotype_data = plyr::join(
  df1,
  df2,
  by="sample",
  type="left") %>%
  distinct(sample,.keep_all = TRUE) %>%
  drop_na()

high_keratinization_sample = phenotype_data %>% filter(cluster == "high_keratinization" ) %>% pull(sample)
low_keratinization_sample = phenotype_data %>% filter(cluster == "low_keratinization" ) %>% pull(sample)
all_sample = phenotype_data$sample


dataFiltrna <- TCGAanalyze_Filtering(
  tabDF = assay(rna[,phenotype_data$barcode],"unstranded"),
  method = "quantile", 
  qnt.cut =  0.25
) 

colnames(dataFiltrna) <- colnames(dataFiltrna)  %>% substr(1,16)
dataFiltmirna <- TCGAanalyze_Filtering(
  tabDF = mirna[,phenotype_data$sample],
  method = "quantile", 
  qnt.cut =  0.25
)

identical(colnames(dataFiltrna),colnames(dataFiltmirna))


dataDEGsmirna <- TCGAanalyze_DEA(
  mat1 = as.matrix(dataFiltmirna[,high_keratinization_sample]),
  mat2 = as.matrix(dataFiltmirna[,low_keratinization_sample]),
  metadata=FALSE,
  Cond1type = "high_keratinization",
  Cond2type = "low_keratinization",
  fdr.cut = 0.01 ,
  logFC.cut= 1,
  method = "glmLRT"
) 
rm(df1,df2,phenotype_data,high_keratinization_sample,low_keratinization_sample)

# find mirna-mrna significant relationship
rna_fc = dataDEGs %>%
  dplyr::select(logFC) %>%
  dplyr::rename(FC1=logFC)
mirna_fc = dataDEGsmirna["logFC"] %>% rename(FC1=logFC)
inter0 <- twoTimePoint(rna_fc, mirna_fc)
outs <- twoTimePointSamp(rna_fc, mirna_fc)

miRanda <- getInputSpecies("Human1",threshold=140)
#Identify miRNA mRNA relationships bellow a P value threshold, default is 0.05
sig_InterR <- threshSigInter(inter0, outs)
#Intersect the mirRanda file with your output results
results <- mirandaIntersectInter(sig_InterR, outs, rna_fc, mirna_fc, miRanda)
gene_id_name_pairs = dataDEGs%>% 
  rownames_to_column("gene_id") %>%
  select(gene_id,gene_name) %>%
  rename(V2=gene_id)
sig_InterR <- sig_InterR%>%
  left_join(gene_id_name_pairs,by="V2")
CorRes<-results$corrs %>%
  left_join(gene_id_name_pairs,by="V2") %>%
  mutate(gene_id=V2)%>%
  mutate(V2=gene_name)
mirRnaHeatmapDiff(CorRes,main= "miRNA-mRNA relationship", legend_breaks = c(0,2,4,6,8,10), 
                  legend_labels = c("0", "2", "4", "6", "8", "Absolute Log\nfold difference")) 
write_csv(CorRes,"mirtarrnaseq_result.csv")
#------------------------------------------------------------------------------ 
cat(unique(sig_InterR$gene_name),sep=",")
cat(unique(sig_InterR$V1),sep=",")
cat(rownames(mirna_fc),sep=",")
cat(dataDEGs$gene_name,sep=",")
cat(unique(CorRes$V2),sep=" ")
cat(unique(CorRes$gene_id),sep=" ")
