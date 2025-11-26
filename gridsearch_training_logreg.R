## Clear workspace
rm(list=ls())

## TeamName
team<-"Example"

## Load all the data so we can quickly combine it and explore it. 
source("load_data.R")
SEPSISdat_train<-load_data("train")
sum(SEPSISdat_train$inhospital_mortality)/nrow(SEPSISdat_train)
SEPSISdat_test<-load_data("test")
sum(SEPSISdat_test$inhospital_mortality)/nrow(SEPSISdat_test)

## Do some data cleaning, e.g. imputation for missing values
library(mice)
library(pROC)
library(glmnet)
library(PRROC)
SEPSISdat_train<-complete(mice(SEPSISdat_train, method = "pmm",m=1))
SEPSISdat_test<-complete(mice(SEPSISdat_test, method = "pmm",m=1))

## Add some derived variables, e.g. shock index= HR/SBP
SEPSISdat_train$SI<-SEPSISdat_train$hr_bpm_adm/SEPSISdat_train$sysbp_mmhg_adm
SEPSISdat_test$SI<-SEPSISdat_test$hr_bpm_adm/SEPSISdat_test$sysbp_mmhg_adm

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


cols <- colnames(SEPSISdat_train)
cols_to_drop <- c("studyid_adm","inhospital_mortality","tobacco_adm",
                  "travelmethod_adm", "traveldist_adm", "waterpure_adm",
                  "cookloc_adm", "lightfuel_adm", "bednet_adm", 
                  "momagefirstpreg_adm", "householdsize_adm",
                  "alivechildren_adm", "deadchildren_adm")

features <- setdiff(cols, cols_to_drop)

# Create Matrix (remove Intercept -1 because glmnet adds its own)
f <- as.formula(paste("~", paste(features, collapse = "+")))

X_train <- model.matrix(f, SEPSISdat_train)[, -1]
y_train <- as.factor(SEPSISdat_train$inhospital_mortality)

X_test  <- model.matrix(f, SEPSISdat_test)[, -1]
y_test  <- SEPSISdat_test$inhospital_mortality

# ------------------------------------------------------------------------------
# 3. Grid Search for Penalized Logistic Regression
# ------------------------------------------------------------------------------

alpha_grid <- seq(0, 1, by = 0.2) 
cw_grid <- list(
  c("0" = 1, "1" = 1),
  c("0" = 1, "1" = 5),
  c("0" = 1, "1" = 10),
  c("0" = 1, "1" = 15)
)

best_model <- NULL
best_score <- -Inf 
best_params <- NULL

cat("Starting Grid Search...\n")

for (cw in cw_grid) {
  weights_vec <- ifelse(y_train == "1", cw["1"], cw["0"])
  
  for (a in alpha_grid) {
    set.seed(123)
    cv_fit <- cv.glmnet(
      x = X_train, y = y_train, family = "binomial",
      alpha = a, weights = weights_vec,
      type.measure = "auc", nfolds = 5
    )
    
    current_score <- max(cv_fit$cvm)
    
    cat("Testing: Alpha =", a, "| Class Weight 1 =", cw["1"], "| CV AUC =", round(current_score, 4), "\n")
    
    if (current_score > best_score) {
      best_score <- current_score
      best_model <- cv_fit
      best_params <- list(alpha = a, weights = cw, lambda = cv_fit$lambda.min)
    }
  }
}

cat("Best CV AUC:", best_score, "\n")
print(best_params)

# ------------------------------------------------------------------------------
# 4. Generate Predictions for Plotting
# ------------------------------------------------------------------------------

# Generate probabilities (need as.vector for pROC)
prob_train <- as.vector(predict(best_model, newx = X_train, s = "lambda.min", type = "response"))
prob_test  <- as.vector(predict(best_model, newx = X_test,  s = "lambda.min", type = "response"))

# Add back to dataframe for convenience (optional, but matches your style)
SEPSISdat_train$probSepsisLR <- prob_train
SEPSISdat_test$probSepsisLR  <- prob_test

# ------------------------------------------------------------------------------
# 5. Plot the AUC & Find Threshold
# ------------------------------------------------------------------------------

# --- Train ROC ---
roc_LR <- roc(inhospital_mortality ~ probSepsisLR, data=SEPSISdat_train)

# Plot Train
plot(roc_LR, main=paste0('Train AUC=', round(roc_LR$auc,3)))

# Find Threshold (Youden)
thresh_full <- coords(roc_LR, "b", best.method="youden", input = "threshold", transpose = T,
                      ret = c("threshold", "sensitivity","specificity","ppv","npv","fp","tp","fn","tn"))

# Pick the threshold value
threshold <- thresh_full[1] # Extracts just the numeric threshold

# --- Test ROC ---
roc_LR_test <- roc(inhospital_mortality ~ probSepsisLR, data=SEPSISdat_test)

# Add Test to Plot
plot(roc_LR_test, add=T, col='red')
text(0.3, 0.3, paste0('Test AUC=', round(roc_LR_test$auc,3)), col="red")


# ------------------------------------------------------------------------------
# 6. Submission Preparation
# ------------------------------------------------------------------------------

final_coefs <- coef(best_model, s = "lambda.min")
final_coefs_matrix <- as.matrix(final_coefs)

myModel <- NULL
myModel$thresh <- round(as.numeric(threshold), 3)
myModel$const  <- final_coefs_matrix[1, 1]
myModel$coeffs <- final_coefs_matrix[-1, 1]
# Remove zero coefficients (features dropped by Lasso)
myModel$coeffs <- myModel$coeffs[myModel$coeffs != 0] 

dput(myModel)

save(best_model, best_params, file=paste0(team, "_Best_LR_Grid.RData"))

# ------------------------------------------------------------------------------
# 7. Quick Performance Evaluation
# ------------------------------------------------------------------------------
source(file.path("..","scoring","evaluate_performance.R"))
res <- NULL

res <- rbind(res, evaluate_model(SEPSISdat_train$inhospital_mortality, SEPSISdat_train$probSepsisLR, threshold, "Training", 0))
res <- rbind(res, evaluate_model(SEPSISdat_test$inhospital_mortality, SEPSISdat_test$probSepsisLR, threshold, "Testing", 0))
print(res)