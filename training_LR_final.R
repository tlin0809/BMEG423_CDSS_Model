setwd("C:\\Users\\Traye Lin's Laptop\\OneDrive - UBC\\Desktop\\BMEG423\\Traye_R_model")


#clear work space
rm(list = ls())
team <- "Team 3"

#all libs

library(mice)

# load data ----
source("load_data.R")
SEPSISdat_train <- load_data("train")
SEPSISdat_test <- load_data("test")
#missing value with mice
SEPSISdat_train <- complete(mice(SEPSISdat_train, method = "pmm", m = 1))
SEPSISdat_test <- complete(mice(SEPSISdat_test, method = "pmm", m = 1))


# process neuro measures (BCS scores) ----
## higher score = more likely to die
process_neuro <- function(df) {
  df$bcseye_adm<- gsub("Watches or follows", "1", df$bcseye_adm)
  df$bcseye_adm <- gsub("Fails to watch or follow", "2", df$bcseye_adm)
  
  df$bcsmotor_adm<- gsub("Localizes painful stimulus", "1", df$bcsmotor_adm)
  df$bcsmotor_adm <- gsub("Withdraws limb from painful stimulus", "2", df$bcsmotor_adm)
  df$bcsmotor_adm<- gsub("No response or inappropriate response", "3", df$bcsmotor_adm)
  
  df$bcsverbal_adm<- gsub("Cries appropriately with pain, or, if verbal, speaks", "1", df$bcsverbal_adm)
  df$bcsverbal_adm <- gsub("Moan or abnormal cry with pain", "2", df$bcsverbal_adm)
  df$bcsverbal_adm <- gsub("No vocal response to pain", "3", df$bcsverbal_adm)
  
  return(df)
}
SEPSISdat_train <- process_neuro(SEPSISdat_train)
SEPSISdat_test <- process_neuro(SEPSISdat_test)

#replace categorical with numerical----
encode_all_categoricals <- function(df) {
  cat_cols <- names(df)[sapply(df, function(x) is.factor(x) || is.character(x))]
  for (col in cat_cols) {
    df[[col]] <- as.numeric(factor(df[[col]]))
  }
  return(df)
}
SEPSISdat_train <- encode_all_categoricals (SEPSISdat_train)
SEPSISdat_test  <- encode_all_categoricals(SEPSISdat_test)

#feature engineering ----
phoenix_score <- function(df) {
  #respiratory
  df$"phoenix_score_respiratory" <- 0
  #some data are not available
  #score of 1 if higher than normal range (https://www.ncbi.nlm.nih.gov/books/NBK596717/table/ch1survey.T.normal_respiratory_rate_by_a/)
  #score of 1 if has Severe respiratory distress (respdistress_adm)
  #score of 1 if any oxygen saturation is 1 if < than 97%
  get_rr_upper <- function(age_months) {
    if (age_months < 1)  return(60)     
    if (age_months < 12) return(60)    
    if (age_months < 120) return(50)   
    if (age_months < 12*18) return(22)   
    return(20)                         
  }
  rr_upper <- sapply(df$agecalc_adm, get_rr_upper)
  df$phoenix_score_respiratory <- df$phoenix_score_respiratory + (df$rr_brpm_app_adm > rr_upper)
  df$phoenix_score_respiratory <- df$phoenix_score_respiratory + (df$respdistress_adm == 1)
  df$phoenix_score_respiratory <- df$phoenix_score_respiratory + ((df$spo2site1_pc_oxi_adm < 97)+(df$spo2site2_pc_oxi_adm < 97))
  
  
  #cardiovascular
  df$"phoenix_score_cardiovascular" <- 0
  lactate_pts <- ifelse(df$lactate_mmolpl_adm >= 11, 2, ifelse(df$lactate_mmolpl_adm >= 5, 1, 0))
  #map alternative
  df$map <- (df$sysbp_mmhg_adm + 2*df$diasbp_mmhg_adm) / 3
  #low 2 points
  df$map_low  <- with(df, ifelse(agecalc_adm < 1,   17,
                                 ifelse(agecalc_adm < 12,  25,
                                        ifelse(agecalc_adm < 24,  31,
                                               ifelse(agecalc_adm < 60,  32,
                                                      ifelse(agecalc_adm < 12*12, 36,
                                                             38))))))
  # in range 1 point                                            
  df$map_high <- with(df, ifelse(agecalc_adm < 1,   30,
                                 ifelse(agecalc_adm < 12,  38,
                                        ifelse(agecalc_adm < 24,  43,
                                               ifelse(agecalc_adm < 60,  44,
                                                      ifelse(agecalc_adm < 12*12, 48,
                                                             51))))))
  
  df$phoenix_score_cardiovascular <- with(df,ifelse(map < map_low, 2, ifelse(map <= map_high, 1, 0)))
  df$phoenix_score_cardiovascular <- df$phoenix_score_cardiovascular + lactate_pts 
  
  df$map_low <- NULL
  df$map_high <- NULL
  
  
  #coagulation
  df$phoenix_hematocrit_score <- 0
  for (i in 1:nrow(df)) {
    age_m <- df$agecalc_adm[i]
    hct   <- df$hematocrit_gpdl_adm[i]
    
    # +2 if <33
    if (hct < 33) {
      df$phoenix_hematocrit_score[i] <- 2
      next
    }
    #1 if below normal
    lower_lim <- ifelse(age_m < 1,   31,     
                        ifelse(age_m < 2,   28,    
                               ifelse(age_m < 6,   29,    
                                      ifelse(age_m < 24,  33,     
                                             ifelse(age_m < 72,  34,     
                                                    ifelse(age_m < 144, 35,     
                                                           ifelse(age_m < 216, 36,     
                                                                  41)))))))
    if (hct < lower_lim) {
      df$phoenix_hematocrit_score[i] <- 1
    }
  }
  
  #neuro
  #lower response has more points
  df$neuro <- df$bcseye_adm+df$bcsmotor_adm+df$bcsverbal_adm
  
  df$curved_neuro<-0
  get_cruved_neuro <- function(neuro_score) {
    if (age_months <= 3)  return(0)     
    if (age_months <=5) return(1)    
    if (age_months <=7) return(2)   
    return(20)                         
  }
  df$curved_neuro <- sapply(df$neuro, get_rr_upper)
  
  
  #total
  df$"phoenix_score_total" <- df$"phoenix_score_respiratory" + df$"phoenix_score_cardiovascular" + df$phoenix_hematocrit_score+df$curved_neuro 
  return(df)
  
}
## get score
SEPSISdat_train <- phoenix_score(SEPSISdat_train)
SEPSISdat_test  <- phoenix_score(SEPSISdat_test)


# PRISM-III APS (another mortality measure) ----
prism_score <- function(df) {
  df$prism_res<-0
  #function for get ages
  get_rr_upper <- function(age_months) {
    if (age_months < 1)  return(100)     
    if (age_months < 12) return(100)    
    if (age_months < 144) return(80)   
    return(60)                         
  }
  rr_upper <- sapply(df$agecalc_adm, get_rr_upper)
  df$prism_res <- df$prism_res+ 2.501*(df$rr_brpm_app_adm > rr_upper)
  
  df$prism_bp_sys<-0
  df$prism_bp_dia<-0
  df$prism_hr<-0
  df$prism_glu<-0
  df$prism_temp<-0
  
  #for loop everything 
  for (i in 1:nrow(df)) {
    age <- df$agecalc_adm[i]
    bp_s <-df$sysbp_mmhg_adm[i]
    bp_a <-df$diasbp_mmhg_adm[i]
    hr <-df$hr_bpm_adm[i]
    glu <-df$glucose_mmolpl_adm[i]
    temp <-df$temp_c_adm[i]
    
    
    #newborn
    if (age < 1) {
      
      #bp sys
      if (bp_s >= 51 & bp_s <= 55) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 5.205} 
      else if (bp_s >= 40 & bp_s <= 50) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 14.010} 
      else if (bp_s < 40) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 60.407} 
      else if (bp_s > 125) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 2.789}
      
      #bp dia
      if (bp_a >80) {df$prism_bp_dia[i] <- df$prism_bp_dia[i] + 4.915} 
      
      #hr
      if (hr >= 195 & hr <= 214) {df$prism_hr[i] <- df$prism_hr[i] + 2.915} 
      else if (hr >= 215 & hr <= 225) {df$prism_hr[i] <- df$prism_hr[i] + 5.214} 
      else if (hr < 75) {df$prism_hr[i] <- df$prism_hr[i] + 3.493} 
      else if (hr > 225) {df$prism_hr[i] <- df$prism_hr[i] + 10.306}
      
      #glucose
      if (glu >= 50 & glu <= 59) {df$prism_glu[i] <- df$prism_glu[i] + 3.332} 
      else if (glu >= 40 & glu <= 49) {df$prism_glu[i] <- df$prism_glu[i] + 11.928} 
      else if (glu >= 30 & glu <= 39) {df$prism_glu[i] <- df$prism_glu[i] + 12.640} 
      else if (glu >= 160 & glu <= 200) {df$prism_glu[i] <- df$prism_glu[i] + 2.312} 
      else if (glu >= 201 & glu <= 250) {df$prism_glu[i] <- df$prism_glu[i] + 4.497} 
      else if (glu >= 251 & glu <= 400) {df$prism_glu[i] <- df$prism_glu[i] + 12.787} 
      else if (glu < 30) {df$prism_glu[i] <- df$prism_glu[i] + 17.445} 
      else if (glu > 400) {df$prism_glu[i] <- df$prism_glu[i] + 12.787}
      
      #temp
      if (temp < 33) {df$prism_temp[i] <- df$prism_temp[i] + 30.940} 
      else if (temp > 40) {df$prism_temp[i] <- df$prism_temp[i] + 5.805}  
      
      
      #infant
    } else if (age < 12) {
      
      #bp sys
      if (bp_s >= 56 & bp_s <= 65) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 5.205} 
      else if (bp_s >= 45 & bp_s <= 50) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 14.010} 
      else if (bp_s < 45) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 60.407} 
      else if (bp_s > 135) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 2.789}
      
      #bp dia
      if (bp_a >95) {df$prism_bp_dia[i] <- df$prism_bp_dia[i] + 4.915}
      
      #hr
      if (hr >= 195 & hr <= 214) {df$prism_hr[i] <- df$prism_hr[i] + 2.915} 
      else if (hr >= 215 & hr <= 225) {df$prism_hr[i] <- df$prism_hr[i] + 5.214} 
      else if (hr < 75) {df$prism_hr[i] <- df$prism_hr[i] + 3.493} 
      else if (hr > 225) {df$prism_hr[i] <- df$prism_hr[i] + 10.306}
      
      #glucose
      if (glu >= 50 & glu <= 59) {df$prism_glu[i] <- df$prism_glu[i] + 3.332} 
      else if (glu >= 40 & glu <= 49) {df$prism_glu[i] <- df$prism_glu[i] + 11.928} 
      else if (glu >= 30 & glu <= 39) {df$prism_glu[i] <- df$prism_glu[i] + 12.640} 
      else if (glu >= 160 & glu <= 200) {df$prism_glu[i] <- df$prism_glu[i] + 2.312} 
      else if (glu >= 201 & glu <= 250) {df$prism_glu[i] <- df$prism_glu[i] + 4.497} 
      else if (glu >= 251 & glu <= 400) {df$prism_glu[i] <- df$prism_glu[i] + 12.787} 
      else if (glu < 30) {df$prism_glu[i] <- df$prism_glu[i] + 17.445} 
      else if (glu > 400) {df$prism_glu[i] <- df$prism_glu[i] + 12.787}
      
      #temp
      if (temp < 33) {df$prism_temp[i] <- df$prism_temp[i] + 30.940} 
      else if (temp > 40) {df$prism_temp[i] <- df$prism_temp[i] + 5.805}  
      
      
      
      #child
    } else if (age < 144) {
      #bp sys
      if (bp_s >= 66 & bp_s <= 75) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 5.205} 
      else if (bp_s >= 55 & bp_s <= 65) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 14.010} 
      else if (bp_s < 55) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 60.407} 
      else if (bp_s > 150) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 2.789}
      
      #bp dia
      if (bp_a >100) {df$prism_bp_dia[i] <- df$prism_bp_dia[i] + 4.915} 
      
      #hr
      if (hr >= 165 & hr <= 184) {df$prism_hr[i] <- df$prism_hr[i] + 2.915} 
      else if (hr >= 185 & hr <= 205) {df$prism_hr[i] <- df$prism_hr[i] + 5.214} 
      else if (hr < 55) {df$prism_hr[i] <- df$prism_hr[i] + 3.493} 
      else if (hr > 205) {df$prism_hr[i] <- df$prism_hr[i] + 10.306}
      
      #glucose
      if (glu >= 50 & glu <= 59) {df$prism_glu[i] <- df$prism_glu[i] + 3.332} 
      else if (glu >= 40 & glu <= 49) {df$prism_glu[i] <- df$prism_glu[i] + 11.928} 
      else if (glu >= 30 & glu <= 39) {df$prism_glu[i] <- df$prism_glu[i] + 12.640} 
      else if (glu >= 160 & glu <= 200) {df$prism_glu[i] <- df$prism_glu[i] + 2.312} 
      else if (glu >= 201 & glu <= 250) {df$prism_glu[i] <- df$prism_glu[i] + 4.497} 
      else if (glu >= 251 & glu <= 400) {df$prism_glu[i] <- df$prism_glu[i] + 12.787} 
      else if (glu < 30) {df$prism_glu[i] <- df$prism_glu[i] + 17.445} 
      else if (glu > 400) {df$prism_glu[i] <- df$prism_glu[i] + 12.787}
      
      #temp
      if (temp < 33) {df$prism_temp[i] <- df$prism_temp[i] + 30.940} 
      else if (temp > 40) {df$prism_temp[i] <- df$prism_temp[i] + 5.805}  
      
      
      
      #adolescent
    } else {
      #bp sys
      if (bp_s >= 76 & bp_s <= 85) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 5.205} 
      else if (bp_s >= 65 & bp_s <= 75) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 14.010} 
      else if (bp_s < 65) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 60.407} 
      else if (bp_s > 190) {df$prism_bp_sys[i] <- df$prism_bp_sys[i] + 2.789}
      
      #bp dia
      if (bp_a >110) {df$prism_bp_dia[i] <- df$prism_bp_dia[i] + 4.915} 
      
      #hr
      if (hr >= 135 & hr <= 144) {df$prism_hr[i] <- df$prism_hr[i] + 2.915} 
      else if (hr >= 145 & hr <= 155) {df$prism_hr[i] <- df$prism_hr[i] + 5.214} 
      else if (hr < 55) {df$prism_hr[i] <- df$prism_hr[i] + 3.493} 
      else if (hr > 155) {df$prism_hr[i] <- df$prism_hr[i] + 10.306}
      
      #glucose
      if (glu >= 50 & glu <= 59) {df$prism_glu[i] <- df$prism_glu[i] + 3.332} 
      else if (glu >= 40 & glu <= 49) {df$prism_glu[i] <- df$prism_glu[i] + 11.928} 
      else if (glu >= 30 & glu <= 39) {df$prism_glu[i] <- df$prism_glu[i] + 12.640} 
      else if (glu >= 160 & glu <= 200) {df$prism_glu[i] <- df$prism_glu[i] + 2.312} 
      else if (glu >= 201 & glu <= 250) {df$prism_glu[i] <- df$prism_glu[i] + 4.497} 
      else if (glu >= 251 & glu <= 400) {df$prism_glu[i] <- df$prism_glu[i] + 12.787} 
      else if (glu < 30) {df$prism_glu[i] <- df$prism_glu[i] + 17.445} 
      else if (glu > 400) {df$prism_glu[i] <- df$prism_glu[i] + 12.787}
      
      #temp
      if (temp < 33) {df$prism_temp[i] <- df$prism_temp[i] + 30.940} 
      else if (temp > 40) {df$prism_temp[i] <- df$prism_temp[i] + 5.805}  
      
    }
    
    
  }
  df$prism_total=df$prism_res+df$prism_bp_sys+
    df$prism_bp_dia+
    df$prism_hr+
    df$prism_glu+
    df$prism_temp
  
  return (df)
  
  
  
}

SEPSISdat_train=prism_score(SEPSISdat_train)
SEPSISdat_test=prism_score(SEPSISdat_test)

# PELOD-2 (another mortality measure)
## An Update of the PEdiatric Logistic Organ Dysfunction Score
PELOD_score <- function(df) {
  df$pelod_cardio<-0
  df$pelod_lac<-0
  
  #for loop everything 
  for (i in 1:nrow(df)) {
    age <- df$agecalc_adm[i]
    map <-df$map[i]
    lac <-df$lactate_mmolpl_adm[i]
    
    #newborn
    if (age < 1) {
      #cardio vescular
      #lac
      if (lac >= 5.0 & lac <= 10.9) {df$pelod_lac[i] <- df$pelod_lac[i] + 1} 
      else if (lac >= 11) {df$pelod_lac[i] <- df$pelod_lac[i] + 4}
      
      #MAP
      if (map >= 31 & map <= 45) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 2} 
      else if (map >= 17 & map <= 30) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 3} 
      else if (map <=16) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 6}
      
      
      
    } else if (age < 11) {
      #cardio vescular
      #lac
      if (lac >= 5.0 & lac <= 10.9) {df$pelod_lac[i] <- df$pelod_lac[i] + 1} 
      else if (lac >= 11) {df$pelod_lac[i] <- df$pelod_lac[i] + 4}
      
      #MAP
      if (map >= 39 & map <= 54) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 2} 
      else if (map>= 25 & map <= 38) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 3} 
      else if (map <=24) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 6}
      
      
      
      
      
    } else if (age < 23) {
      #cardio vescular
      #lac
      if (lac >= 5.0 & lac <= 10.9) {df$pelod_lac[i] <- df$pelod_lac[i] + 1} 
      else if (lac >= 11) {df$pelod_lac[i] <- df$pelod_lac[i] + 4}
      
      #MAP
      if (map >= 44 & map <= 59) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 2} 
      else if (map >= 31 & map <= 43) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 3} 
      else if (map <=30) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 6}
      
      
      
      
      
    }else if (age < 59) {
      #cardio vescular
      #lac
      if (lac >= 5.0 & lac <= 10.9) {df$pelod_lac[i] <- df$pelod_lac[i] + 1} 
      else if (lac >= 11) {df$pelod_lac[i] <- df$pelod_lac[i] + 4}
      
      #MAP
      if (map >= 46 & map <= 61) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 2} 
      else if (map >= 32 & map <= 44) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 3} 
      else if (map <=31) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 6}
      
      
    }else if (age < 143) {
      #cardio vescular
      #lac
      if (lac >= 5.0 & lac <= 10.9) {df$pelod_lac[i] <- df$pelod_lac[i] + 1} 
      else if (lac >= 11) {df$pelod_lac[i] <- df$pelod_lac[i] + 4}
      
      #MAP
      if (map >= 49 & map <= 64) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 2} 
      else if (map >= 36 & map <= 48) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 3} 
      else if (map <=35) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 6}
      
      
    } else {
      #cardio vescular
      #lac
      if (lac >= 5.0 & lac <= 10.9) {df$pelod_lac[i] <- df$pelod_lac[i] + 1} 
      else if (lac >= 11) {df$pelod_lac[i] <- df$pelod_lac[i] + 4}
      
      #MAP
      if (map >= 52 & map <= 66) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 2} 
      else if (map >= 38 & map <= 51) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 3} 
      else if (map <=37) {df$pelod_cardio[i] <- df$pelod_cardio[i] + 6}
      
    }
  }
  df$pelod_total<-df$pelod_cardio+df$pelod_lac+df$pelod_lac
  
  
  return (df)
}
SEPSISdat_train=PELOD_score(SEPSISdat_train)
SEPSISdat_test=PELOD_score(SEPSISdat_test)


# keep only group 1

feature_group_1 <- c(
  "agecalc_adm", 
  "lactate_mmolpl_adm","hr_bpm_adm","rr_brpm_app_adm","sysbp_mmhg_adm","diasbp_mmhg_adm",
  "temp_c_adm","spo2site1_pc_oxi_adm","spo2site2_pc_oxi_adm","spo2onoxy_adm","glucose_mmolpl_adm","hematocrit_gpdl_adm",
  "hctpretransfusion_adm","caprefill_adm","respdistress_adm","malariastatuspos_adm","hivstatus_adm","bcseye_adm",
  "bcsmotor_adm","bcsverbal_adm", 
  "inhospital_mortality",
  "phoenix_score_respiratory","phoenix_score_cardiovascular","phoenix_hematocrit_score", "map",
  "curved_neuro", "phoenix_score_total", "neuro",
  "prism_res","prism_sys","prism_dia","prism_hr","prism_glu","prism_temp", "prism_total",
  "pelod_cardio","pelod_lac" ,"pelod_total"
)

SEPSISdat_train <- SEPSISdat_train  [ , intersect(feature_group_1, names(SEPSISdat_train ))]
SEPSISdat_test<- SEPSISdat_test [ , intersect(feature_group_1, names(SEPSISdat_train ))]


library(dplyr)

## Build a regression model
cols <- colnames(SEPSISdat_train)
# Drop the tobacco_adm variable as the observations change which cause issues
cols <- cols[!(cols %in% c("studyid_adm","inhospital_mortality","tobacco_adm"))]

formula <- as.formula(paste0("inhospital_mortality ~ ",paste0(cols,sep="",collapse="+")))
logReg <- glm(formula,data=SEPSISdat_train,family=binomial(link='logit'))
summary(logReg)


source(file.path("..","scoring","evaluate_performance.R"))

# threshold sweep ----
cutoffs <- seq(-6,-4,0.001)
accuracy_train <- NULL
accuracy_test<- NULL
x<-NULL
## Quick but not necessarily great way to find a threshold
SEPSISdat_train$probSepsisLR <- predict(logReg,newdata=subset(SEPSISdat_train,select=-c(inhospital_mortality)))
SEPSISdat_test$probSepsisLR <- predict(logReg,newdata=subset(SEPSISdat_test,select=-c(inhospital_mortality)))


for (i in cutoffs){
  res_train <- try(
    evaluate_model(SEPSISdat_train$inhospital_mortality,
                   SEPSISdat_train$probSepsisLR,
                   i, "Training", 0),
    silent = TRUE
  )
  res_test <- try(
    evaluate_model(SEPSISdat_test$inhospital_mortality,
                   SEPSISdat_test$probSepsisLR,
                   i, "Testing", 0),
    silent = TRUE
  )
  if (!inherits(res_train, "try-error") && !inherits(res_test, "try-error")) {
    accuracy_train <- c(accuracy_train, res_train$weighted_score)
    accuracy_test  <- c(accuracy_test,  res_test$weighted_score)
    x <- c(x,  i)
  } else {
    #message("skipping cutoff = ", i, "some error")
  }
  
}

# see if the 2 are the same, if so we are not over fitting
#find index of bet fit 
max(accuracy_train)
index=match(max(accuracy_train),accuracy_train)
x[index]
index=match(max(accuracy_test),accuracy_test)
x[index]


#plot to see weighted score vs threshold 
plot(x, accuracy_train, pch =19,type='b',col= "steelblue",
     main ="LR Train", xlab="Cutoff Level", ylab = "Accuracy %")

plot(x, accuracy_test, pch =19,type='b',col= "steelblue",
     main ="LR Test", xlab="Cutoff Level", ylab = "Accuracy %")


# Plot the AUC
library('pROC')
roc_LR <- roc(inhospital_mortality ~ probSepsisLR,data=SEPSISdat_train)
plot(roc_LR,main=paste0('LR Train AUC=',round(roc_LR$auc,3)))


thresh<-coords(roc_LR, "b", best.method="youden", input = "threshold", transpose = T,
               ret = c("threshold", "sensitivity","specificity","ppv","npv","fp","tp","fn","tn"))

roc_LR_test <- roc(inhospital_mortality ~ probSepsisLR,data=SEPSISdat_test)
plot(roc_LR_test,add=T,col='red')
text(0.3,0.3,paste0('AUC_test=',round(roc_LR_test$auc,3)),col="red")
threshold<-thresh[1]


## Report the values to put into my get_sepsis_score's load_sepsis_model function
myModel<- NULL
myModel$thresh <- round(x[index],3)
myModel$const <- logReg$coefficients[1]
myModel$coeffs <- logReg$coefficients[2:length(logReg$coefficients)]
dput(myModel)

# save model ----
save(logReg, file=paste0(team,"_","LR.RData"))

# evaluate performance ----
source(file.path("..","scoring","evaluate_performance.R"))
res<-NULL
best_threshold <-cutoffs[index]
res<-rbind(res,evaluate_model(SEPSISdat_train$inhospital_mortality,SEPSISdat_train$probSepsisLR,best_threshold ,"Training",0))
res<-rbind(res,evaluate_model(SEPSISdat_test$inhospital_mortality,SEPSISdat_test$probSepsisLR,best_threshold ,"Testing",0))
print(res)


