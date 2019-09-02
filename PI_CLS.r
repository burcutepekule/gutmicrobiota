#! /usr/bin/env Rscript
rm(list=ls())

args <- commandArgs(TRUE)
print(args)

sampIdx <- as.integer(args[1])
numTree <- as.integer(args[2])

library(readr)
library(party)

DATA_REG_TRAIN_CLS_B            <- read_csv(paste0("./R/PI_DATA_5E",sampIdx,"/DATA_REG_TRAIN_CLS_BLNCD.csv"))
DATA_REG_TRAIN_CLS_B$prevelance <- as.factor(DATA_REG_TRAIN_CLS_B$prevelance)
DATA_REG_TRAIN_CLS_B            <- DATA_REG_TRAIN_CLS_B

tree_cls_b.cf                  <- cforest(prevelance ~ ., data = DATA_REG_TRAIN_CLS_B, control = cforest_unbiased(mtry = 4, ntree = numTree))
varImpNormal_cls_b             <- varimp(tree_cls_b.cf)
tableNormal                    <- rbind(varImpNormal_cls_b)

write.table(tableNormal,paste0("./R/PI_RESULTS/tableNormal_T",numTree,"_CLS_S_5E",sampIdx,".txt"), sep="\t")

varImpCond_cls_b   <- varimp(tree_cls_b.cf, conditional = TRUE)
tableCond          <- rbind(varImpCond_cls_b)

write.table(tableCond,paste0("./R/PI_RESULTS/tableCond_T",numTree,"_CLS_S_5E",sampIdx,".txt"), sep="\t")





