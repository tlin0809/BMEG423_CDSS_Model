# Team: Team 3

# Author: Traye Lin 

# Version: RF Model 1

# set working directory (for myself)
setwd("C:\\Users\\Traye Lin's Laptop\\OneDrive - UBC\\Desktop\\BMEG423\\Traye_R_model")

## Clear workspace
rm(list = ls())

## TeamName
team <- "Team 3"

## Load helper function to read data
source("load_data.R")

# load datasets
SEPSISdat_train <- load_data("train")
SEPSISdat_test  <- load_data("test")

# print mortality rate - should be the same as training_gbm.R (nothing changed)
cat("Training mortality rate:", mean(SEPSISdat_train$inhospital_mortality), "\n")
cat("Testing mortality rate:",  mean(SEPSISdat_test$inhospital_mortality), "\n")

# data cleaning - following training_gbm.R (nothing changed)
library(mice)

SEPSISdat_train <- complete(mice(SEPSISdat_train, method = "pmm", m = 1))
SEPSISdat_test  <- complete(mice(SEPSISdat_test,  method = "pmm", m = 1))


# add derived variables (Shock Index); nothing changed
SEPSISdat_train$SI <- SEPSISdat_train$hr_bpm_adm / SEPSISdat_train$sysbp_mmhg_adm
SEPSISdat_test$SI  <- SEPSISdat_test$hr_bpm_adm  / SEPSISdat_test$sysbp_mmhg_adm


# RF with class weights
## weight_minority = N_total / (2 * N_minority)
wt_minority <- nrow(SEPSISdat_train) / (2 * mean(SEPSISdat_train$inhospital_mortality)*nrow(SEPSISdat_train))
wt_majority <- nrow(SEPSISdat_train) / (2 * (1-mean(SEPSISdat_train$inhospital_mortality))*nrow(SEPSISdat_train))

## imbalance ratio 
ratio <- wt_minority / wt_majority # 23
cw <- c("0" = 1, "1" = 23) # potential class weights; fine-tune later

# Random Forest model using ranger
library(ranger)

## increase num. of trees and decrease node size to compensate for variances and adjust for specificity
rf_model2 <- ranger(
  formula = as.factor(inhospital_mortality) ~ .,
  data = SEPSISdat_train,
  probability = TRUE,
  num.trees = 1500,
  mtry = floor(sqrt(ncol(SEPSISdat_train))),
  min.node.size = 5,           # allow deeper trees
  importance = "impurity",
  class.weights = cw
)


# generate predictions, choose classification threshold
train_pred2 <- rf_model2$predictions[, 2]   # column 2 is P(mortality=1)
test_pred2  <- predict(rf_model2, SEPSISdat_test)$predictions[, 2]

summary(train_pred2) # check distribution 

# ROC, determine threshold by Youden index (same as provided gbm model)
library(pROC)

roc_RF_train2 <- roc(SEPSISdat_train$inhospital_mortality, train_pred2) # training AUC
roc_RF_test2 <- roc(SEPSISdat_test$inhospital_mortality, test_pred2) # testing AUC 

# plot ROC curvce
plot(roc_RF_train2, main = paste0("RF Training AUC = ", round(roc_RF_train2$auc, 3)))
plot(roc_RF_test2, main = paste0("RF Testing AUC = ", round(roc_RF_test2$auc, 3)))

threshold2 <- coords(
  roc_RF_train2, 
  "best", 
  best.method = "youden", # J = sensitivity + specificity - 1
  ret = "threshold"
)

cat("Threshold:", as.numeric(threshold2), "\n")


# save model #1
saveRDS(rf_model2, paste0(team, "_RF_model.2.rds"))

# values for load_sepsis_model()
myModel2 <- list(
  thresh = round(as.numeric(threshold2), 4),
  model_path = paste0(team, "_RF_model.2.rds")
)

dput(myModel2) # print model threshold & path

# model performance 
source(file.path("..","scoring","evaluate_performance.R"))
res<-NULL

res <- rbind(
  res,
  evaluate_model(SEPSISdat_train$inhospital_mortality,
                 train_pred2,
                 as.numeric(threshold),
                 "Training",
                 0) # training
)
res <- rbind(
  res,
  evaluate_model(SEPSISdat_test$inhospital_mortality,
                 test_pred2,
                 as.numeric(threshold),
                 "Testing",
                 0) # testing
)


print(res)