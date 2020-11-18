rm(list = ls()); options(stringsAsFactors = F)

library(data.table); library(caret); library(ROCR)
library(glmnet); library(PRROC); library(randomForest)
library(mlbench); library(e1071)
setwd("d:/data/github/")

#Preprocessing Train And Test set (PTAT)
PTAT = function(expression.set, demographic.set, feature){
  rownames(expression.set) = expression.set[, 1]
  expression.set = expression.set[, -1]
  dat = cbind(demographic.set$Dx, t(expression.set[feature, ]))
  dat = data.frame(dat, stringsAsFactors = F)
  colnames(dat)[1] = "dx"
  return(dat)
}

#Features
load("feature_anm1.rdata")

#Loading samples with diagnosis (label)
#"AD_sets_demo.rdata" is randomly made, not real data.
load("AD_sets_demo.rdata")

#Preparing_dataset
train = fread("ANM1_corrected.txt")
test = fread("ADNI_corrected.txt")
train = data.frame(train)
test = data.frame(test)
rownames(train) = train$sample
rownames(test) = test$sample

#Explanation for "expr_sample_ADNI_ANM1_ANM2.rdata"
#We made another sample data that could be runned in this code.
#Note that these made sets were not real sets.
load("expr_sample_ADNI_ANM1_ANM2.rdata")

#DEG (SVM, RF)
train = expr.anm1
test = expr.adni
train1 = PTAT(expression.set = train, demographic.set = adni, deg_anm1)
test1 = PTAT(expression.set = test, demographic.set = anm1, deg_anm1)
SVM_DEG = svm(as.factor(dx) ~ ., data=train1, kernel="radial", probability = TRUE)
pred = predict(SVM_DEG, test1,  decision.values = TRUE, probability = TRUE)
prob = attr(pred, "probabilities")[,2]
fg = prob[test1$dx == 1]
bg = prob[test1$dx == 0]
roc = roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
AUC_SVM = roc$auc

#
RF_DEG = randomForest(as.factor(dx) ~ ., data=train1, mtry = floor(sqrt(ncol(train1)-1)), ntree = 500, importance = T)
prob = predict(RF_DEG, test1, type = "prob")[, 2]
fg <- prob[test1$dx == 1]
bg <- prob[test1$dx == 0]
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
AUC_RF = roc$auc
performance = c("DEG", AUC_SVM, AUC_RF)

#TF (SVM, RF)
train1 = PTAT(expression.set = train, demographic.set = anm1, feature = tf_anm1)
test1 = PTAT(expression.set = test, demographic.set = adni, feature = tf_anm1)
SVM_TF = svm(as.factor(dx) ~ ., data=train1, kernel="radial", probability = TRUE)
pred = predict(SVM_TF, test1,  decision.values = TRUE, probability = TRUE)
prob = attr(pred, "probabilities")[,2]
fg = prob[test1$dx == 1]
bg = prob[test1$dx == 0]
roc = roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
AUC_SVM = roc$auc

#
RF_TF = randomForest(as.factor(dx) ~ ., data=train1, mtry = floor(sqrt(ncol(train1)-1)), ntree = 500, importance = T)
prob = predict(RF_TF, test1, type = "prob")[, 2]
fg <- prob[test1$dx == 1]
bg <- prob[test1$dx == 0]
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
AUC_RF = roc$auc
temp = c("DEG+TF", AUC_SVM, AUC_RF)
performance = rbind(performance, temp)


#CFG (SVM, RF)
train1 = PTAT(expression.set = train, demographic.set = anm1, feature = cfg_anm1)
test1 = PTAT(expression.set = test, demographic.set = adni, feature = cfg_anm1)
SVM_CFG = svm(as.factor(dx) ~ ., data=train1, kernel="radial", probability = TRUE)
pred = predict(SVM_CFG, test1,  decision.values = TRUE, probability = TRUE)
prob = attr(pred, "probabilities")[,2]
fg = prob[test1$dx == 1]
bg = prob[test1$dx == 0]
roc = roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
AUC_SVM = roc$auc

#
RF_CFG = randomForest(as.factor(dx) ~ ., data=train1, mtry = floor(sqrt(ncol(train1)-1)), ntree = 500, importance = T)
prob = predict(RF_CFG, test1, type = "prob")[, 2]
fg <- prob[test1$dx == 1]
bg <- prob[test1$dx == 0]
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
AUC_RF = roc$auc
temp = c("DEG+CFG", AUC_SVM, AUC_RF)
performance = rbind(performance, temp)

colnames(performance) = c("Features", "SVM", "RF")
performance = data.frame(performance)
performance$Test_set = "ADNI"
performance1 = performance

#Preparing_dataset
test = fread("ANM2_corrected.txt"); test = data.frame(test, stringsAsFactors = F)
rownames(test) = test$sample

#Loading sample datasets
load("expr_sample_ADNI_ANM1_ANM2.rdata")

#DEG (SVM, RF)
test1 = PTAT(expression.set = expr.anm2, demographic.set = anm2, feature = deg_anm1)
pred = predict(SVM_DEG, test1,  decision.values = TRUE, probability = TRUE)
prob = attr(pred, "probabilities")[,2]
fg = prob[test1$dx == 1]
bg = prob[test1$dx == 0]
roc = roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
AUC_SVM = roc$auc

#
prob = predict(RF_DEG, test1, type = "prob")[, 2]
fg <- prob[test1$dx == 1]
bg <- prob[test1$dx == 0]
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
AUC_RF = roc$auc
performance = c("DEG", AUC_SVM, AUC_RF)

#TF (SVM, RF)
test1 = PTAT(expression.set = expr.anm2, demographic.set = anm2, feature = tf_anm1)
pred = predict(SVM_TF, test1,  decision.values = TRUE, probability = TRUE)
prob = attr(pred, "probabilities")[,2]
fg = prob[test1$dx == 1]
bg = prob[test1$dx == 0]
roc = roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
AUC_SVM = roc$auc

#
prob = predict(RF_TF, test1, type = "prob")[, 2]
fg <- prob[test1$dx == 1]
bg <- prob[test1$dx == 0]
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
AUC_RF = roc$auc
temp = c("DEG+TF", AUC_SVM, AUC_RF)
performance = rbind(performance, temp)


#CFG (SVM, RF)
test1 = PTAT(expression.set = expr.anm2, demographic.set = anm2, feature = cfg_anm1)
SVM_CFG = svm(as.factor(dx) ~ ., data=train1, kernel="radial", probability = TRUE)
pred = predict(SVM_CFG, test1,  decision.values = TRUE, probability = TRUE)
prob = attr(pred, "probabilities")[,2]
fg = prob[test1$dx == 1]
bg = prob[test1$dx == 0]
roc = roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
AUC_SVM = roc$auc

#
prob = predict(RF_CFG, test1, type = "prob")[, 2]
fg <- prob[test1$dx == 1]
bg <- prob[test1$dx == 0]
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
AUC_RF = roc$auc
temp = c("DEG+CFG", AUC_SVM, AUC_RF)
performance = rbind(performance, temp)

colnames(performance) = c("Features", "SVM", "RF")
performance = data.frame(performance)
performance$Test_set = "ANM2"

final_AUC = rbind(performance1, performance)
