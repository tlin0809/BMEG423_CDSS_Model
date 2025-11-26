# Team: Team 3

# Version: GBM model 3 v1 (BCS + manual drop features)

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

# drop features ----
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

# GBM model with reduced features ----
library('lightgbm')
gbm_model3 <- lightgbm(
  data = as.matrix(subset(SEPSISdat_train,select=-c(inhospital_mortality, 
                                                    travelmethod_adm, traveldist_adm,
                                                    waterpure_adm, cookloc_adm, 
                                                    lightfuel_adm, tobacco_adm,
                                                    bednet_adm, momagefirstpreg_adm, 
                                                    householdsize_adm,
                                                    alivechildren_adm, 
                                                    deadchildren_adm)))
  , label = SEPSISdat_train$inhospital_mortality
  , objective = "binary"
  , nrounds=10
)

# plot AUC ----
train_pred3 <- predict(gbm_model3, train_mat)
test_pred3  <- predict(gbm_model3, test_mat)

library('pROC')
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
lgb.save(gbm_model3, paste0(team, "_lightgbm_model3.model"))

myModel <- list(
  thresh = round(as.numeric(threshold3), 4),
  model_path = paste0(team, "_lightgbm_model3.model"),
  dropped_features = dropFeats
)

dput(myModel)


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
