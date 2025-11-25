# Team: Team 3

# Author: Traye Lin 

# Version: RF Model 3

# set working directory (change this on personal laptop)
# setwd("C:\\Users\\Traye Lin's Laptop\\OneDrive - UBC\\Desktop\\BMEG423\\Traye_R_model")
setwd("C:\\Users\\traye.lin\\Downloads\\BMEG423_CDSS_Assignment_2025\\r_modelBuilding")

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



# grid search (on full list of features)

## weight_minority = N_total / (2 * N_minority) (skip this)
wt_minority <- nrow(SEPSISdat_train) / (2 * mean(SEPSISdat_train$inhospital_mortality)*nrow(SEPSISdat_train))
wt_majority <- nrow(SEPSISdat_train) / (2 * (1-mean(SEPSISdat_train$inhospital_mortality))*nrow(SEPSISdat_train))

ratio <- wt_minority / wt_majority # imbalance ratio ~= 1:23

################################################
## class weights (skip this after running it once)

cw_grid <- list(
  c("0" = 1, "1" = 10),
  c("0" = 1, "1" = 15),
  c("0" = 1, "1" = 18),
  c("0" = 1, "1" = 20),
  c("0" = 1, "1" = 23)   # true imbalance ratio
)

mtry_grid <- c(
  floor(sqrt(length(imp_features))),
  floor(length(imp_features) / 3)
)

min_node_grid <- c(3, 4, 5, 6, 7, 8, 9, 10)

trees_grid <- c(800, 900, 1000, 1100, 1200, 1300, 1400, 1500)

best_model <- NULL
best_score <- -Inf   # AUPRC is the objective
best_params <- NULL

library(pROC)
library(PRROC)

source(file.path("..", "scoring", "evaluate_performance.R"))

for (cw in cw_grid) {
  for (mtry_val in mtry_grid) {
    for (node_val in min_node_grid) {
      for (ntree in trees_grid) {
        
        cat("Testing:",
            "cw=", cw["1"],
            "mtry=", mtry_val,
            "node=", node_val,
            "trees=", ntree, "\n")
        
        rf_tmp <- ranger(
          formula = as.factor(inhospital_mortality) ~ .,
          data = SEPSIS_train_reduced,
          probability = TRUE,
          num.trees = ntree,
          mtry = mtry_val,
          min.node.size = node_val,
          importance = "impurity",
          class.weights = cw
        )
        
        pred_tmp <- rf_tmp$predictions[,2]
        aucpr_tmp <- pr.curve(scores.class0 = pred_tmp,
                              weights.class0 = SEPSIS_train_reduced$inhospital_mortality)$auc.integral
        
        if (aucpr_tmp > best_score) {
          best_score <- aucpr_tmp
          best_model <- rf_tmp
          best_params <- list(cw=cw, mtry=mtry_val, node=node_val, trees=ntree)
        }
      }
    }
  }
}

cat("Best parameters: \n")
print(best_params) # cw: 1:15; mtry: 21; min_node: 8; n_trees: 1100
cat("Best AUPRC: ", best_score, "\n") # 0.2977

rf_model3 <- best_model
################################################


################################################
## continue here if: obtained optimal parameters already  
## skip here if: ran above grid search for the first time
cw <- c("0" = 1, "1" = 15) 

mtry <- 21

min_node <- 8

n_trees <- 1100

rf_model3 <- ranger(
  formula = as.factor(inhospital_mortality) ~ .,
  data = SEPSIS_train_reduced,
  probability = TRUE,
  num.trees = n_trees,
  mtry = mtry,
  min.node.size = min_node,
  importance = "impurity",
  class.weights = cw
)
################################################

# predictions + threshold 
train_pred3 <- rf_model3$predictions[,2]
test_pred3  <- predict(rf_model3, SEPSISdat_test)$predictions[,2]

roc_train3 <- roc(SEPSISdat_train$inhospital_mortality, train_pred3)
roc_test3  <- roc(SEPSISdat_test$inhospital_mortality,  test_pred3)

plot(roc_train3, main=paste0("RF Model 3 — Train AUC=", round(roc_train3$auc,3)))
plot(roc_test3, col="red", add=TRUE)
text(0.4,0.4,paste0("Test AUC=", round(roc_test3$auc,3)), col="red")

threshold3 <- coords(
  roc_train3,
  "best",
  best.method = "youden",
  ret = "threshold"
)

cat("Selected threshold:", as.numeric(threshold3), "\n")


# save model
saveRDS(rf_model3, paste0(team, "_RF_model.3.rds"))

myModel3 <- list(
  thresh = round(as.numeric(threshold3), 4),
  model_path = paste0(team, "_RF_model.3.rds"),
  selected_features = imp_features
)

dput(myModel3)

# model performance 
source(file.path("..","scoring","evaluate_performance.R"))
res <- NULL

res <- rbind(
  res,
  evaluate_model(
    SEPSISdat_train$inhospital_mortality,
    train_pred3,
    as.numeric(threshold3),
    "Training", 0
  )
)

res <- rbind(
  res,
  evaluate_model(
    SEPSISdat_test$inhospital_mortality,
    test_pred3,
    as.numeric(threshold3),
    "Testing", 0
  )
)

print(res)