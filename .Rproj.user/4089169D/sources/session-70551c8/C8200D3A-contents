###############################################################################
# load_sepsis_model.R
# Loads the Random Forest model and threshold for sepsis prediction.
###############################################################################

load_sepsis_model <- function() {
  
  # Load the trained model (ensure file name matches your saved file)
  model <- readRDS("Team 3_RF_model.rds")
  
  # These values are generated during training and printed via dput(myModel)
  model_info <- list(
    thresh = 0.XXX,                # <-- REPLACE with your exact printed threshold
    model = model
  )
  
  return(model_info)
}
