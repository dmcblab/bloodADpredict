rm(list = ls());options(stringsAsFactors = F)
#######################################################
library(dplyr); library(data.table); library(plyr); library(limma)
setwd("d:/data/github/")

####
dat1 = fread("ANM2_corrected.txt", header = T, stringsAsFactors = F)
dat1 = data.frame(dat1)
rownames(dat1) = dat1[, 1]
dat1 <- dat1[, -1]

#expr value to log2(Fold change)
mean = rowMeans(dat1[, -1])
dat1 <- log2(dat1 / mean)


####
load("AD_sets_demo.rdata")
dat1_demographic_matrix <- model.matrix( ~ anm2$Dx + anm2$age + anm2$sex)
colnames(dat1_demographic_matrix) <- c("Intercept","Disease", "Age", "Gender")

limma_fit_cv1 <- lmFit(dat1, dat1_demographic_matrix)
limma_fit_cv1 <- eBayes(limma_fit_cv1)
deg_cv1 <- topTable(limma_fit_cv1, coef=2 , number=3000, p.value=1, sort.by = "logFC")
write.table(deg_cv1, "DEG_ANM2.txt", row.names = F, col.names = T, sep = "\t")

