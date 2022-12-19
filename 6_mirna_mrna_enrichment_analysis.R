library(tidyverse)
library(gprofiler2)

df = read_csv("mirtarrnaseq_result.csv")

mirnas = unique(df$V1)

query_gene_name = list()

for (i in mirnas) {
  print(i)
  temp_df = df %>% filter(V1==i)
  identifier = i
  query_gene_name[[identifier]] = unique(temp_df$gene_name)
  rm(temp_df)
}

gostres <- gost(query = query_gene_name,organism="hsapiens",
                sources = c("GO","KEGG","REAC","WP"),
                as_short_link = F,
                evcodes=T)

gostreslink <- gost(query = query_gene_name,organism="hsapiens",
                sources = c("GO","KEGG","REAC","WP"),
                as_short_link = T,
                evcodes=T)
browseURL(gostreslink)

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