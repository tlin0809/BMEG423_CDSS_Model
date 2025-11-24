# Team: Team 3

# Author: Traye Lin 

# Version: load_sepsis_model v1

# load model and threshold for sepsis prediction
load_sepsis_model <- function() {
  model <- readRDS("Team 3_RF_model.1.rds") # replace with model selected 
  model_info <- list(
    threshold = 0.0652, # replace with threshold determined 
    model = model
  )
  
  return(model_info) # return model + threshold
}
