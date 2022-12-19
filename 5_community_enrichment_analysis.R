# initialization -----------------------------------------------------------

library(tidyverse)
library(gprofiler2)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(org.Hs.eg.db)

working_dir = "/home/yusri/Documents/project/lusc_keratinization/"
setwd(working_dir)
genome <- "hg38"
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
# data load---------------------------------------------------------------------------
rna <- loadRData("rna_with_cluster.rda")
universe_gene <- bitr(rowData(rna)$gene_name,fromType = "SYMBOL",
                      toType="ENTREZID", 
                      OrgDb="org.Hs.eg.db",
                      drop=F) %>%
  pull(ENTREZID)
dataDEGs <- read.csv("DEG_analysis_result.csv",row.names=1)


gene_id_name_pairs = dataDEGs %>%
  rownames_to_column("gene_id") %>%
  dplyr::select(c("gene_id","gene_name")) %>%
  mutate(entrez_id = bitr(gene_name,fromType = "SYMBOL",
                          toType="ENTREZID", 
                          OrgDb="org.Hs.eg.db",
                          drop=F)$ENTREZID)

df = read_csv("arboreto_node_df.csv") %>%
  dplyr::rename(gene_name=gene_id) %>%
  dplyr::left_join(gene_id_name_pairs,by="gene_name")

# community enrichment analysis----------------------------------------------
count_membership = dplyr::count(df,membership)

selected_membership = c(0,1,2,3,4,5)
cumsum(count_membership[selected_membership,]$n)/sum(count_membership$n)

query_gene_member <- function(df,member){
  res <-df %>%
    filter(membership == member) %>%
    dplyr::select(gene_id,gene_name,entrez_id)
  
  return(res)
}
query_gene_entrez_id = list()
query_gene_name = list()
query_gene_ensembl_id = list()

for (i in selected_membership) {
  print(i)
  query <- query_gene_member(df,i)
  identifier = paste0("Community",as.character(i))
  query_gene_entrez_id[[identifier]] = query$entrez_id
  query_gene_name[[identifier]] = query$gene_name
  query_gene_ensembl_id[[identifier]] = query$gene_id
}

query_gene_entrez_id = query_gene_entrez_id[order(lengths(query_gene_entrez_id),decreasing=TRUE)]
query_gene_name = query_gene_name[order(lengths(query_gene_name),decreasing=TRUE)]
query_gene_ensembl_id = query_gene_ensembl_id[order(lengths(query_gene_ensembl_id),decreasing=TRUE)]


# functional enrichment using g:profiler----------------------------------------
gostres <- gost(query = query_gene_ensembl_id,organism="hsapiens",
                sources = c("GO","KEGG","REAC","WP"),
                as_short_link = F,
                evcodes=T)


# Save all significant result
write_csv(gostres$result,"significant_gsea.csv")

# functional enrichment using clusterprofiler----------------------------------------
create.enrichResult <-function(g){
  gp_mod = g$result[,c("query", "source", "term_id",
                             "term_name", "p_value", "query_size", 
                             "intersection_size", "term_size", 
                             "effective_domain_size", "intersection")]
  gp_mod$GeneRatio = paste0(gp_mod$intersection_size,  "/", gp_mod$query_size)
  
  gp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)
  names(gp_mod) = c("Cluster", "Category", "ID", "Description", "p.adjust", 
                    "query_size", "Count", "term_size", "effective_domain_size", 
                    "geneID", "GeneRatio", "BgRatio")
  gp_mod$geneID = gsub(",", "/", gp_mod$geneID)
  #row.names(gp_mod) = gp_mod$ID
  

  # define as enrichResult object
  gp_mod_enrich  = new("enrichResult", result = gp_mod)
  return(gp_mod_enrich)
}

create.compareClusterResult <-function(g){
  gp_mod = g$result[,c("query", "source", "term_id",
                       "term_name", "p_value", "query_size", 
                       "intersection_size", "term_size", 
                       "effective_domain_size", "intersection")]
  gp_mod$GeneRatio = paste0(gp_mod$intersection_size,  "/", gp_mod$query_size)
  
  gp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)
  names(gp_mod) = c("Cluster", "Category", "ID", "Description", "p.adjust", 
                    "query_size", "Count", "term_size", "effective_domain_size", 
                    "geneID", "GeneRatio", "BgRatio")
  gp_mod$geneID = gsub(",", "/", gp_mod$geneID)
  #row.names(gp_mod) = gp_mod$ID
  
  
  # define as enrichResult object
  gp_mod_cluster = new("compareClusterResult", compareClusterResult = gp_mod)
  return(gp_mod_cluster)
}



# Reactome pathway enrichment using g:profiler----------------------------------------

gostres_reac <- gost(query = query_gene_ensembl_id,organism="hsapiens",
                    sources = c("REAC"),
                    as_short_link = F,
                    evcodes=T)
write_csv(gostres_reac$result,"significant_reactome_enrichment.csv")

ego_reac = create.compareClusterResult(gostres_reac)


dotplot(ego_reac)+
  ggtitle("Reactome enrichment analysis") +
  theme(axis.text.y = element_text(size = 10)) +
  scale_x_discrete(name=NULL,
                   labels=c("Community0\n(145)"= "Community0",
                            "Community1\n(117)"= "Community1",
                            "Community2\n(64)"= "Community2",
                            "Community3\n(67)"= "Community3",
                            "Community4\n(75)"= "Community4",
                            "Community5\n(57)"= "Community5"))
f
ggplot_build(f)$layout$panel_params[[1]]$x$get_labels()
