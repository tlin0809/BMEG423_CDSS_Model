# Team 3 

# Version: 1st GBM submission  


# recompute BCS scores like druing training
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


# called by score_submission.R
get_sepsis_score = function(data, myModel) {
  
  library(mice)
  data <- complete(mice(data, method = "pmm", m = 1)) # impute missing values 

  data$SI <- data$hr_bpm_adm / data$sysbp_mmhg_adm # shock index
  
  data <- bcs_scores(data)   # BCS scores
  
  dropFeats <- myModel$drop_feats # drop features (removed during training)
  
  data_model <- data[, !(names(data) %in% dropFeats)]
  
  data_mat <- as.matrix(data_model)   # convert to matrix
  
  probSepsis <- predict(myModel$bst, data_mat) # predict with GBM 
  
  label <- probSepsis >= myModel$thresh  # convert to binary label 
  
  return(data.frame(probSepsis = probSepsis, label = label)) # return dataframe 
}

# load model
load_sepsis_model <- function() {
  library(lightgbm)

  bst <- lgb.load("Team 3_lightgbm_model3.model") # load GBM model 3 file saved in Model 3 training
  
  myModel <- list(
    thresh = 0.1016,   # determined from GBM model 3 v1
    bst = bst,
    drop_feats = c(
      "inhospital_mortality",
      "travelmethod_adm", "traveldist_adm",
      "waterpure_adm", "cookloc_adm", 
      "lightfuel_adm", "tobacco_adm",
      "bednet_adm", "momagefirstpreg_adm", 
      "householdsize_adm",
      "alivechildren_adm", "deadchildren_adm"
    ) # dropped features
  )
  
  return(myModel)
}
