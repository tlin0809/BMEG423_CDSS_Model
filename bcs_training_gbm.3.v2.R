# Team: Team 3

# Version: GBM model 3 v2 (BCS + manual drop features + hyperparameter search)

# set working directory (for myself; can skip)
setwd("C:\\Users\\Traye Lin's Laptop\\OneDrive - UBC\\Desktop\\BMEG423\\Traye_R_model")

rm(list = ls())

team <- "Team 3"

# data preparation ----
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

# BCS ---- 
## Blantyre coma scale (BCS) scores

bcs_scores <- function(data) { 
  data$bcs_eye_score <- ifelse(
    data$bcseye_adm == "Watches or follows", 2,
    ifelse(data$bcseye_adm == "Eyes open to voice", 1, 0)
  )
  
  data$bcs_motor_score <- ifelse(
    data$bcsmotor_adm == "Localizes painful stimulus", 2,
    ifelse(data$bcsmotor_adm == "Withdraws limb from pain", 1, 0)
  )
  
  data$bcs_verbal_score <- ifelse(
    data$bcsverbal_adm == "Cries appropriately with pain, or, if verbal, speaks", 2,
    ifelse(data$bcsverbal_adm == "Moans or cries inappropriately to pain", 1, 0)
  )
  
  data$bcs_total <- data$bcs_eye_score +
    data$bcs_motor_score +
    data$bcs_verbal_score
  
  data$coma_flag <- ifelse(data$bcs_total <= 2, 1, 0)
  
  return(data)
  
}

SEPSISdat_train <- bcs_scores(SEPSISdat_train)
SEPSISdat_test  <- bcs_scores(SEPSISdat_test)

# feature pruning - drop features 
## feature engineering to be continued
dropFeats <- c(
  "inhospital_mortality",
  "travelmethod_adm", "traveldist_adm",
  "waterpure_adm", "cookloc_adm", 
  "lightfuel_adm", "tobacco_adm",
  "bednet_adm", "momagefirstpreg_adm", 
  "householdsize_adm",
  "alivechildren_adm", "deadchildren_adm"
)

trainX <- SEPSISdat_train[, !(names(SEPSISdat_train) %in% dropFeats)]
testX  <- SEPSISdat_test[,  !(names(SEPSISdat_test)  %in% dropFeats)]

trainY <- SEPSISdat_train$inhospital_mortality
testY  <- SEPSISdat_test$inhospital_mortality

# convert to matrix for lightgbm
train_mat <- as.matrix(trainX)
test_mat  <- as.matrix(testX) 

# Hyperparameter search ---- 
## grid search for AUC optimization

library(lightgbm)
library(PRROC)

## weight_minority = N_total / (2 * N_minority) (skip this; for myself)
wt_minority <- nrow(SEPSISdat_train) / (2 * mean(SEPSISdat_train$inhospital_mortality)*nrow(SEPSISdat_train))
wt_majority <- nrow(SEPSISdat_train) / (2 * (1-mean(SEPSISdat_train$inhospital_mortality))*nrow(SEPSISdat_train))
ratio <- wt_minority / wt_majority # imbalance ratio ~= 1:23

## balance class
pos_weight_true <- sum(trainY == 0) / sum(trainY == 1)   # ~23

scale_pos_grid <- c(10, 15, 18, 20, 23) # 23 is true imbalance ratio; decrease weighting to control false positives
learning_rate_grid <- c(0.05, 0.1) # from smooth to faster
num_leaves_grid <- c(20, 31, 40) # explore smaller and larger trees 
# max_depth_grid <- c(-1, 6) # compare free vs. controlled tree depth
nrounds_grid <- c(150, 300, 500)

best_model <- NULL
best_score <- -Inf
best_params <- NULL

for (w in scale_pos_grid) {
  for (lr in learning_rate_grid) {
    for (leaves in num_leaves_grid) {
      for (depth in max_depth_grid) {
        for (nr in nrounds_grid) {
          
          cat(
            "Testing with: ",
            "scale_pos_weight =", w, 
            "| learning_rate =", lr,
            "| num_leaves =", leaves,
            "| max_depth =", depth,
            "| nrounds =", nr, 
            "\n"
          )
          
          params <- list(
            objective = "binary",
            metric = "binary_logloss",
            learning_rate = lr,
            num_leaves = leaves,
            max_depth = depth,
            scale_pos_weight = w,
            feature_pre_filter = FALSE
          )
          
          dtrain <- lgb.Dataset(train_mat, label = trainY)
          
          gbm_tmp <- lightgbm(
            data = dtrain,
            params = params,
            nrounds = nr,
            verbose = -1
          )
          
          pred_tmp <- predict(gbm_tmp, train_mat)
          
          aucpr_tmp <- pr.curve(
            scores.class0 = pred_tmp,
            weights.class0 = trainY
          )$auc.integral
          
          if (aucpr_tmp > best_score) {
            best_score <- aucpr_tmp
            best_model <- gbm_tmp
            best_params <- list(
              scale_pos_weight = w,
              learning_rate = lr,
              num_leaves = leaves,
              depth = depth,
              nrounds = nr
            )
          }
        }
      }
    }
  }
}

cat("Best parameters found:\n")
print(best_params)
cat("Best AUPRC:", best_score, "\n")

gbm_model3 <- best_model



########################
# manual fix overfitting ----
params <- list(
  objective = "binary",
  metric = "auc",
  
  # parameters found by grid search
  scale_pos_weight = 10,
  learning_rate    = 0.05,
  num_leaves       = 20,
  max_depth        = 6,        # replaced -1; prevent infinite depth
  nthread          = 4,
  
  # prevent overfitting 
  lambda_l1        = 0.3,
  lambda_l2        = 0.3,
  min_data_in_leaf = 30,
  feature_fraction = 0.8,
  bagging_fraction = 0.8,
  bagging_freq     = 5
)

trainX <- as.matrix(trainX)
testX  <- as.matrix(testX)

gbm_model3 <- lightgbm(
  data  = trainX,
  label = SEPSISdat_train$inhospital_mortality,
  params = params,
  nrounds = 150,
  verbose = 1
)
########################

# plot AUC ----
library(pROC)

train_pred3 <- predict(gbm_model3, train_mat)
test_pred3  <- predict(gbm_model3, test_mat)

roc_train3 <- roc(trainY, train_pred3)
roc_test3  <- roc(testY,  test_pred3)

plot(roc_train3, main=paste0("GBM Model 3 Train AUC=", round(roc_train3$auc,3)))
plot(roc_test3, col="red", add=TRUE)
text(0.4,0.4,paste0("Test AUC=", round(roc_test3$auc,3)), col="red")

## unchanged threshold determination (by youden index)
threshold3 <- coords(
  roc_train3,
  "best",
  best.method = "youden",
  ret = "threshold"
)

cat("Selected threshold:", as.numeric(threshold3), "\n")


# save model ----

lgb.save(gbm_model3, paste0(team, "_lightgbm_model3.txt"))

myModel3 <- list(
  thresh = round(as.numeric(threshold3), 4),
  model_path = paste0(team, "_lightgbm_model3.txt"),
  dropped_features = dropFeats
)

dput(myModel3)


# evaluate performance ----

source(file.path("..","scoring","evaluate_performance.R"))

res <- NULL

res <- rbind(
  res,
  evaluate_model(trainY, train_pred3, as.numeric(threshold3), "Training", 0)
)
res <- rbind(
  res,
  evaluate_model(testY, test_pred3, as.numeric(threshold3), "Testing", 0)
)

print(res)