rm(list = ls()); options(stringsAsFactors = F)
library(data.table)
#setting your directory consisting of data.
setwd("d:/data/github/")

#ADNI_Gene_Expression_Profile.csv: expression data
#ADNI_Gene_Expression_Profile.csv must be needed and publicly available data at ADNI site.
#DXSUM_PDXCONV_ADNIALL.csv: data for patient with diagnosis
#DXSUM_PDXCONV_ADNIALL.csv must be needed and publicly available data at ADNI site.
dat1 = fread("ADNI_Gene_Expression_Profile.csv")
dat1 = dat1[, -748]
pts = fread("DXSUM_PDXCONV_ADNIALL.csv")

#
sample = t(dat1[c(2, 3), -(1:3)])
date = sample[, 1]
sample = sample[, 2]
sample = strsplit(sample, split = "_")
sample1 = NULL
for(i in 1:length(sample)){
  sample1 = c(sample1, sample[[i]][3])
}
sample1 = as.numeric(sample1)

#
pts2 = pts[, c("RID", "VISCODE", "VISCODE2", "DXCHANGE", "DXCURREN")]
idx = which(!is.na(match(pts2$RID, sample1)))
pts2 = pts2[idx, ]

#
for(i in 1:length(sample1)){
  idx1 = which(pts2$RID == sample1[i])
  pts2[idx1, "expression_date"] = date[i]
}

#pts3_manual is data for patient without diagnosis at the time of examination for expression data.
pts3_auto = NULL
pts3_manual = NULL
for(i in 1:length(sample1)){
  temp = pts2[pts2$RID == sample1[i], ]
  idx1 = which(temp$VISCODE == temp$expression_date)
  idx2 = which(temp$VISCODE2 == temp$expression_date)
  idx = union(idx1, idx2)
  if(length(idx) > 0){
    dx = c(temp$DXCHANGE[idx], temp$DXCURREN[idx])
    dx = dx[!is.na(dx)]
    temp[idx, "dx_at_expr"] = dx
    pts3_auto = rbind(pts3_auto, temp)
  } else {
    pts3_manual = rbind(pts3_manual, temp)
  }
}

pts3_auto = pts3_auto[!is.na(pts3_auto$dx_at_expr), ]

#If you mannually annotate diagnosis for pts3_manual, I recommend you "pts4 = rbind(pts3_auto, pts3_manual)".
pts4 = pts3_auto
pts4$final_diagnosis = pts4$dx_at_expr
pts4$final_diagnosis[pts4$dx_at_expr == 1] = "NL"
pts4$final_diagnosis[pts4$dx_at_expr == 2] = "MCI"
pts4$final_diagnosis[pts4$dx_at_expr == 3] = "AD"
pts4$final_diagnosis[pts4$dx_at_expr == 4] = "MCI"
pts4$final_diagnosis[pts4$dx_at_expr == 5] = "AD"
pts4$final_diagnosis[pts4$dx_at_expr == 6] = "AD"
pts4$final_diagnosis[pts4$dx_at_expr == 7] = "NL"
pts4$final_diagnosis[pts4$dx_at_expr == 8] = "MCI"
pts4$final_diagnosis[pts4$dx_at_expr == 9] = "NL"
pts4 = pts4[, c("RID", "final_diagnosis")]
write.table(pts4, "adni_pts_info.txt", sep = "\t", row.names = F)
