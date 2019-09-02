#! /usr/bin/env Rscript
rm(list=ls())

args <- commandArgs(TRUE)
print(args)

sampIdx <- as.integer(args[1])
numTree <- as.integer(args[2])

library(readr)
library(party)

DATA_REG_TRAIN_POSONLY <- read_csv(paste0("./R/PI_DATA_5E",sampIdx,"/DATA_REG_TRAIN_POSONLY.csv"))

tree_posonly.cf        <- cforest(prevelance ~ ., data = DATA_REG_TRAIN_POSONLY, control = cforest_unbiased(mtry = 4, ntree = numTree))
varImpNormal_posonly   <- varimp(tree_posonly.cf)
tableNormal            <- rbind(varImpNormal_posonly)

write.table(tableNormal,paste0("./R/PI_RESULTS/tableNormal_T",numTree,"_REG_S_5E",sampIdx,".txt"), sep="\t")

varImpCond_posonly <- varimp(tree_posonly.cf ,conditional = TRUE)
tableCond          <- rbind(varImpCond_posonly)

write.table(tableCond,paste0("./R/PI_RESULTS/tableCond_T",numTree,"_REG_S_5E",sampIdx,".txt"), sep="\t")



