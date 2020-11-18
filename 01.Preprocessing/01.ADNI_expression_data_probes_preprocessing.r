rm(list = ls())
library("AnnotationDbi"); library("hgu219.db"); library(data.table)
setwd("d:/data/github/")

#ADNI_Gene_Expression_Profile.csv: expression data
#ADNI_Gene_Expression_Profile.csv must be needed and publicly available data at ADNI site.
data1 <- fread("ADNI_Gene_Expression_Profile.csv")
data1 <- data.frame(data1)
data1 <- data1[, -748]
probe_ID = unlist(data1[-(1:9), 1])
platform_out = select(hgu219.db, keys= probe_ID, columns=c("SYMBOL", "ENTREZID"),keytype="PROBEID")

######## protein coding gene filter
#Homo_sapiens.GRCh38.94.gtf could be downloaded at http://asia.ensembl.org/Homo_sapiens/Info/Index.
info = rtracklayer::import("Homo_sapiens.GRCh38.94.gtf")
index = which(info$gene_biotype=="protein_coding")
length(index)

info <- info[index]
coding_gene_list <- unique(info$gene_name)
coding_gene_id <- unique(info$gene_id)

length(coding_gene_list)
length(coding_gene_id)

coding_gene_list <- coding_gene_list[!coding_gene_list==""]

coding_platform_out <- NULL
for (i in 1:nrow(platform_out)){
  if (i %% 2000 == 0){
    print(i)
  }
  if (platform_out$SYMBOL[i] %in% coding_gene_list)
    coding_platform_out <- rbind(coding_platform_out,platform_out[i,])
}
nrow(coding_platform_out)

write.table(coding_platform_out, "affymetrix_u219_protein_coding.txt", col.names = T, row.names = F, sep = "\t")
