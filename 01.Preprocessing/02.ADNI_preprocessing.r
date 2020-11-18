rm(list = ls());options(stringsAsFactors = F)
library(data.table); library(plyr)
setwd("d:/data/github/")

#Note that following data-sets were not real-data, but just sample data-sets for running source codes.
load("expr_sample_ADNI.rdata")
idx1 <- dat_dx$RID[dat_dx$final_diagnosis == "NL"]
idx2 <- dat_dx$RID[dat_dx$final_diagnosis == "AD"]
ctrl = length(idx1)
ad = length(idx2)

data1 <- data1[, c("PROBEID", idx1, idx2)]

#ADNI, expr, chr to numeric
for(i in 2:ncol(data1)){
  data1[, i] <- as.numeric(data1[, i])
}

#removing probes with low-expression
value = unlist(data1[, -1])
median = median(value)
idx1 <- NULL
for(i in 1:nrow(data1)){
  temp <- which(data1[i, ] > median)
  if (length(temp) > round((ncol(data1)-1)/2)){
    idx1 <- c(idx1, i)
  }
  if (i %% 1000 == 0){
    print(i)
  }
}
data1 <- data1[idx1, ]

#probes to gene
#The coding_platform_out.rdata could be made by running 01.ADNI_expression_data_probes_preprocessing.r.
load("coding_platform_out.rdata")
data1 <- join(coding_platform_out, data1, by = "PROBEID")
data1 = data1[!is.na(data1[, 5]), ]
data1 = data1[, -1]
data1$SYMBOL = as.character(data1$SYMBOL)

#
data1$mean = rowMeans(data1[, 2:ncol(data1)])
data1 = data1[c(order(data1$mean)), ]

list.tmp = split(data1, data1[, 1], drop = T)

dat2 = NULL
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
dat2 = dat2[, -(ncol(dat2))]
dat2 = data.frame(dat2)
colnames(dat2) = names(data1)[-ncol(data1)]

#integrating disease status label (AD or CN) to expression dataset
sample1 <- rep(0, ctrl)
sample2 <- rep(1, ad)
dx <- c("dx", sample1, sample2)
dat2 <- rbind(dx, dat2)
write.table(dat2, "ADNI.txt", sep = '\t', col.names = T, row.names = F)



