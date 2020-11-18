rm(list = ls());options(stringsAsFactors = F)
library(data.table); library(readxl); library(samr)

#
dat1 <- fread("ADNI_corrected.txt"); dat1 = data.frame(dat1)

#
symbol <- dat1[, 1]
dat1 = dat1[, -1]

#
load("AD_sets_demo.rdata")
dx = adni$Dx+1

samfit01 <- SAM(as.matrix(dat1), dx, resp.type="Two class unpaired", geneid = symbol, 
                nperms = 500, logged2 = FALSE, fdr.output=1)
a01_up <- samfit01$siggenes.table$genes.up
a01_down <- samfit01$siggenes.table$genes.lo
a001_all <- rbind(a01_up, a01_down)
write.table(a001_all, "DEG_ADNI.txt", col.names = T, row.names = F, sep = "\t")

