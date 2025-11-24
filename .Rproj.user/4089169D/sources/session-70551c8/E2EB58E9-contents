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

# Random Forest model using ranger
library(ranger) # may need to install ranger if you haven't

rf_model <- ranger(
  formula = as.factor(inhospital_mortality) ~ .,
  data = SEPSISdat_train,
  probability = TRUE,        # gives class probabilities
  num.trees = 800,           # higher trees improve AUC
  mtry = floor(sqrt(ncol(SEPSISdat_train))),  # RF rule
  min.node.size = 10,        # prevents overfitting
  importance = "impurity"
) # model training 


# generate predictions, choose classification threshold
train_pred <- rf_model$predictions[, 2]   # column 2 is P(mortality=1)
test_pred  <- predict(rf_model, SEPSISdat_test)$predictions[, 2]

summary(train_pred) # check distribution 

# ROC, determine threshold by Youden index (same as provided gbm model)
library(pROC)

roc_RF_train <- roc(SEPSISdat_train$inhospital_mortality, train_pred) # training AUC
roc_RF_test <- roc(SEPSISdat_test$inhospital_mortality, test_pred) # testing AUC 

# plot ROC curvce
plot(roc_RF_train, main = paste0("RF Training AUC = ", round(roc_RF_train$auc, 3)))
plot(roc_RF_test, main = paste0("RF Testing AUC = ", round(roc_RF_test$auc, 3)))

threshold <- coords(
  roc_RF_train, 
  "best", 
  best.method = "youden", # J = sensitivity + specificity - 1
  ret = "threshold"
)

cat("Threshold:", as.numeric(threshold), "\n")


# save model #1
saveRDS(rf_model, paste0(team, "_RF_model.1.rds"))

# values for load_sepsis_model()
myModel <- list(
  thresh = round(as.numeric(threshold), 4),
  model_path = paste0(team, "_RF_model.1.rds")
)

dput(myModel) # print model threshold & path

# model performance 
source(file.path("..","scoring","evaluate_performance.R"))
res<-NULL

res <- rbind(
  res,
  evaluate_model(SEPSISdat_train$inhospital_mortality,
                 train_pred,
                 as.numeric(threshold),
                 "Training",
                 0) # training
)
res <- rbind(
  res,
  evaluate_model(SEPSISdat_test$inhospital_mortality,
                 test_pred,
                 as.numeric(threshold),
                 "Testing",
                 0) # testing
)
  

print(res)
