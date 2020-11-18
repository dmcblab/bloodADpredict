rm(list = ls())
#################
library(sva); library(data.table) ; library(plyr)
setwd("d:/data/github/")

dat1 <- fread("GSE63060.txt", header = F)
dat2 <- fread("GSE63061.txt", header = F)
dat3 <- fread("ADNI.txt", header = F)
dat1 <- data.frame(dat1)
dat2 <- data.frame(dat2)
dat3 <- data.frame(dat3)

#
colnames(dat1)[1] = "symbol"
colnames(dat2)[1] = "symbol"
colnames(dat3)[1] = "symbol"
dat4 <- join(dat1, dat2, by = "symbol")
dat4 <- join(dat4, dat3, by = "symbol")

#removing_NA_values
for(i in 2:ncol(dat4)){
  dat4 = dat4[!is.na(dat4[, i]), ]
}

#Extract symbol
symbol <- dat4[-1, 1]
dat4 = dat4[, -1]

#Extract disease status
Alzheimer = unlist(dat4[1, ])
dat4 = dat4[-1, ]

#Setting_Batch
a = rep(1, (ncol(dat1) -1))
b = rep(2, (ncol(dat2) -1))
c = rep(3, (ncol(dat3) -1))
batch <- c(a, b, c)

#
for(i in 1:ncol(dat4)){
  dat4[, i] <- as.numeric(dat4[, i])
  print(i)
}

#
combat_edata = ComBat(dat = as.matrix(dat4), batch = batch, par.prior=TRUE, prior.plots=FALSE)

#dx <- as.numeric(alzheimer)
dat1 = combat_edata[, batch == 1]
dat2 = combat_edata[, batch == 2]
dat3 = combat_edata[, batch == 3]
dat1 = cbind(symbol, dat1)
dat2 = cbind(symbol, dat2)
dat3 = cbind(symbol, dat3)
colnames(dat1) = c("symbol", paste0("ADNI.", 1:10))
colnames(dat2) = c("symbol", paste0("ANM1.", 1:10))
colnames(dat3) = c("symbol", paste0("ANM2.", 1:10))

demo.adni = data.frame(ID = colnames(dat1)[-1], Dx = unlist(Alzheimer[batch == 1]))
demo.anm1 = data.frame(ID = colnames(dat2)[-1], Dx = unlist(Alzheimer[batch == 2]))
demo.anm2 = data.frame(ID = colnames(dat3)[-1], Dx = unlist(Alzheimer[batch == 3]))

#
write.table(dat1, "ADNI_corrected.txt", sep = '\t', row.names = F)
write.table(dat2, "ANM1_corrected.txt", sep = '\t', row.names = F)
write.table(dat3, "ANM2_corrected.txt", sep = '\t', row.names = F)

#
write.table(demo.adni, "ADNI_demo.txt", sep = '\t', row.names = F)
write.table(demo.anm1, "ANM1_demo.txt", sep = '\t', row.names = F)
write.table(demo.anm2, "ANM2_demo.txt", sep = '\t', row.names = F)

