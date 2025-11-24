# Team: Team 3

# Author: Traye Lin 

# Version: get_sepsis_score v1


# compute sepsis risk probability and labels for patients
get_sepsis_score <- function(data, model) {
  
  data$SI <- data$hr_bpm_adm / data$sysbp_mmhg_adm # shock index
  
  # probability predictions using selected model (e.g. RF, GBM)
  prediction <- predict(model$model, data = data)$predictions[, 2]
  
  # convert probabilities into labels (binary) using the stored threshold
  labels <- ifelse(prediction > model$threshold, 1, 0)
  
  return(list(prob = prediction, label = labels)) # return output
}
