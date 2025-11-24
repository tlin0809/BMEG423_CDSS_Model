###############################################################################
# get_sepsis_score.R
# Generates sepsis risk probability and binary label for each patient.
###############################################################################

get_sepsis_score <- function(data, model) {
  
  # --------------------------------------------------------------------------
  # 1. Reproduce preprocessing used in training
  # --------------------------------------------------------------------------
  
  # Add shock index (must match training preprocessing exactly)
  data$SI <- data$hr_bpm_adm / data$sysbp_mmhg_adm
  
  # --------------------------------------------------------------------------
  # 2. Generate probability predictions using the ranger model
  # --------------------------------------------------------------------------
  preds <- predict(model$model, data = data)$predictions[, 2]
  
  # --------------------------------------------------------------------------
  # 3. Convert probabilities into binary labels using stored threshold
  # --------------------------------------------------------------------------
  labels <- ifelse(preds > model$thresh, 1, 0)
  
  # --------------------------------------------------------------------------
  # 4. Return output in required format
  # --------------------------------------------------------------------------
  return(list(prob = preds, label = labels))
}
