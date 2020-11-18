rm(list = ls())
library(data.table); library(plyr)
setwd("d:/data/github")

#Note that following data-sets were not real-data, but just sample data-sets for running source codes.
load("expr_sample_GSE63060.rdata")
idx1 <- dat_dx$RID[dat_dx$final_diagnosis == "NL"]
idx2 <- dat_dx$RID[dat_dx$final_diagnosis == "AD"]
ctrl = length(idx1)
ad = length(idx2)
dat1 <- dat1[, c("ID", idx1, idx2)]

#expression value: chr to num
for(i in 2:ncol(dat1)){
  dat1[, i] <- as.numeric(dat1[, i])
}

#
load("probes_symbol_gse63060.rdata")
dat1 = join(p2s, dat1, by = "ID")
idx1 = which(dat1$ILMN_Gene != "")
dat1 = dat1[idx1, ]
idx1 = which(!is.na(dat1[, 5]))
dat1 = dat1[idx1, ]
dat1 = dat1[, -1]

#selecting a probe with median value.(multiple probes annotating a gene)
dat1$mean = rowMeans(dat1[, 2:ncol(dat1)])
dat1 = dat1[c(order(dat1$mean)), ]

list.tmp = split(dat1, dat1[, 1], drop = T)
dat2 <- NULL
for(i in 1:length(list.tmp)){
  temp <- data.frame((list.tmp[i]))
  row_num <- nrow(temp)
  if(row_num < 2){
    temp <- temp
  }else if(row_num %% 2 == 0){
    row_num2 <- row_num / 2
    temp_sub1 <- temp[(row_num2:(row_num2+1)), ]
    temp <- temp[1, ]
    temp[, 2:ncol(temp_sub1)] <- colMeans(temp_sub1[, 2:ncol(temp_sub1)], na.rm = FALSE, dims = 1)
  }else if(row_num %% 2 != 0){
    row_num2 <- row_num / 2
    row_num2 <- ceiling(row_num2)
    temp <- temp[row_num2, ]
  }
  colnames(temp) <- NA
  dat2 <- rbind(dat2, temp)
  if(i %% 1000 == 0){
    print(i)
  }
  
}
dat2 <- dat2[, -ncol(dat2)]
colnames(dat2) = names(dat1)[-ncol(dat1)]

#integrating disease status label (AD or CN) to expression dataset
sample1 <- rep(0, ctrl)
sample2 <- rep(1, ad)
dx <- c("dx", sample1, sample2)
dat2 <- rbind(dx, dat2)
write.table(dat2, "GSE63060.txt", sep = '\t', col.names = T, row.names = F)


