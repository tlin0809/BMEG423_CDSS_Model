# Team: Team 3

# Version: GBM model 3 v3 (BCS + feature engineering; no hyperparameter search)

# set working directory (for myself; can skip)
setwd("C:\\Users\\Traye Lin's Laptop\\OneDrive - UBC\\Desktop\\BMEG423\\Traye_R_model")
setwd("C:\\Users\\traye.lin\\Downloads\\BMEG423_CDSS_Assignment_2025\\r_modelBuilding")

rm(list = ls())

team <- "Team 3"

# data preparation ----
source("load_data.R")

# load datasets
SEPSISdat_train <- load_data("train")
SEPSISdat_test  <- load_data("test")

# print mortality rate
cat("Training mortality rate:", mean(SEPSISdat_train$inhospital_mortality), "\n")
cat("Testing mortality rate:",  mean(SEPSISdat_test$inhospital_mortality), "\n")

# load all required packages 
library(mice)
library(dplyr)
library(lightgbm)
library(pROC)

# imputation
SEPSISdat_train <- complete(mice(SEPSISdat_train, method = "pmm", m = 1))
SEPSISdat_test  <- complete(mice(SEPSISdat_test,  method = "pmm", m = 1))

# BCS score ----
## Blantyre Coma Scale 
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

# feature engineering (shock index, strong features, etc.) ----
add_clinical_features <- function(data) {
  data %>%
    # shock index + SIPA
    mutate(
      shock_index = hr_bpm_adm / sysbp_mmhg_adm,
      sipa_month_flag = case_when(
        agecalc_adm < 4  & shock_index > 2.0 ~ 1,   # 0–3 months
        agecalc_adm >= 4 & agecalc_adm < 7  & shock_index > 1.5 ~ 1,  # 4–6 months
        agecalc_adm >= 7 & agecalc_adm < 12 & shock_index > 1.4 ~ 1,  # 7–11 months
        agecalc_adm >= 12 & agecalc_adm < 84  & shock_index > 1.2 ~ 1, # 1–6 years
        agecalc_adm >= 84 & agecalc_adm < 156 & shock_index > 1.0 ~ 1, # 7–12 years
        agecalc_adm >= 156 & shock_index > 0.9 ~ 1,                    # 13+ years
        TRUE ~ 0
      )
    ) %>%
    # lactate
    mutate(
      lactate_level = lactate_mmolpl_adm,
      lactate_month_flag = case_when(
        agecalc_adm <= 2  & lactate_level > 3.5 ~ 1,   # 0–2 months
        agecalc_adm >= 3 & agecalc_adm <= 24 & lactate_level > 3.3 ~ 1, # 3–24 months
        agecalc_adm > 24 & lactate_level > 2.4 ~ 1,    # >24 months
        TRUE ~ 0
      )
    ) %>%
    # respiratory rate
    mutate(
      respiratory_rate = rr_brpm_app_adm,
      RR_flag = case_when(
        agecalc_adm < 12 & respiratory_rate < 20 ~ 1,  # slow breathing
        agecalc_adm < 2  & respiratory_rate >= 60 ~ 1, # fast 0–2 months
        agecalc_adm >= 2  & agecalc_adm < 12  & respiratory_rate >= 50 ~ 1, # 2–11 months
        agecalc_adm >= 12 & agecalc_adm < 60 & respiratory_rate >= 40 ~ 1,  # 1–5 years
        agecalc_adm >= 60 & agecalc_adm < 144 & respiratory_rate >= 30 ~ 1, # 5–12 years
        agecalc_adm >= 144 & respiratory_rate >= 18 ~ 1,                    # 12+ years
        respiratory_rate >= 70 ~ 1,                                         # very fast
        TRUE ~ 0
      )
    ) %>%
    # temperature
    mutate(
      temperature = temp_c_adm,
      temp_flag = case_when(
        temperature < 36.5 ~ 1,
        temperature > 37.5 ~ 1,
        TRUE ~ 0
      )
    ) %>%
    # glucose
    mutate(
      glucose = glucose_mmolpl_adm,
      gluc_flag = case_when(
        agecalc_adm <= 1 & glucose < 2.5 ~ 1,   # 0–1 month
        agecalc_adm > 1  & glucose > 3.33 ~ 1,  # >1 month
        TRUE ~ 0
      )
    ) %>%
    # capillary refill (--> numeric)
   mutate(
      caprefill_sec = as.numeric(caprefill_adm),   # ensure numeric seconds
      
      caprefill_flag = case_when(
        # neonates: up to 7 days old (agecalc_adm is in months; 7 days ~= 0.23 months)
        agecalc_adm <= 0.23 & caprefill_sec > 7 ~ 1,     # neonate abnormal
        agecalc_adm <= 0.23 & caprefill_sec <= 7 ~ 0,    # neonate normal
        
        # non-neonates:
        caprefill_sec >= 3 ~ 1,   # abnormal refill (>3 sec)
        TRUE ~ 0
      )
    ) %>%
    # respiratory distress
    mutate(
      respdistress_adm = as.numeric(respdistress_adm == "TRUE")
    ) %>%
    # Malaria
    mutate(
      malariastatuspos_adm = as.numeric(malariastatuspos_adm == "TRUE")
    ) %>%
    # HIV
    mutate(
      hivstatus_adm = as.numeric(hivstatus_adm == "TRUE")
    ) %>%
    # Comorbidity score (0–5)
    mutate(
      comorbid_score =
        (comorbidity_adm_airway == "TRUE") +
        (comorbidity_adm_cardiac == "TRUE") +
        (comorbidity_adm_sicklecell == "TRUE") +
        (comorbidity_adm_tuberculosis == "TRUE") +
        (comorbidity_adm_disability == "TRUE")
    ) %>%
    # MAP + shock flag
    mutate(
      mean_arterial_pressure = (sysbp_mmhg_adm + 2 * diasbp_mmhg_adm) / 3,
      min_MAP = 40 + 1.5 * (agecalc_adm / 12),
      MAP_flag = case_when(
        mean_arterial_pressure < min_MAP ~ 1,
        TRUE ~ 0
      )
    )
}

SEPSISdat_train <- add_clinical_features(SEPSISdat_train)
SEPSISdat_test  <- add_clinical_features(SEPSISdat_test)

# drop features ----
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
gbm_model3 <- lightgbm(
  data = train_mat,
  label = trainY,
  objective = "binary",
  nrounds = 10
)

# plot AUC ----
train_pred3 <- predict(gbm_model3, train_mat)
test_pred3  <- predict(gbm_model3, test_mat)

roc_train3 <- roc(trainY, train_pred3)
roc_test3  <- roc(testY,  test_pred3)

plot(roc_train3, main = paste0("GBM Model 3 Train AUC = ", round(roc_train3$auc, 3)))
plot(roc_test3, col = "red", add = TRUE)
text(0.4, 0.4, paste0("Test AUC = ", round(roc_test3$auc, 3)), col = "red")

# threshold (Youden) ----
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
source(file.path("..", "scoring", "evaluate_performance.R"))

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