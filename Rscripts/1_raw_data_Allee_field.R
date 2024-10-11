# raw data list 


Allee_data <- list()

#Fig1c
load("./Rdata/aten_insitu_sperm_conc.RData")
Allee_data[["Fig1c_top_raw"]] <- aten_insitu_sperm_conc 
#Fig1c
load("./Rdata/adig_insitu_sperm_conc.RData")
Allee_data[["Fig1c_bot_raw"]] <- adig_insitu_sperm_conc 
#Fig1d
load("./Rdata/Fig1d_raw.RData") # fert_data
Allee_data[["Fig1d_raw"]] <- fert_data 
#Fig1e
load('./Rdata/ahyac_insitu_sperm_conc.RData') #ahyac_insitu_sperm_conc
Allee_data[["Fig1e_raw"]] <- ahyac_insitu_sperm_conc 
#Controls
load("./Rdata/controls1.RData") #controls
Allee_data[["Controls"]] <- controls 
#Colony density surveys reefflat
load("./Rdata/coral_den_reefflat.RData") #coral_den_reefflat
Allee_data[["Col_density_surveys_reefflat"]] <- coral_den_reefflat 
##Colony spacing surveys reefslope and flat
load("C:/Users/gerar/Documents/1_R_projects/2023-palau/Rdata/intercol_dist_natural.RData") #intercol_dist_natural
Allee_data[["intercol_dist_natural"]] <- intercol_dist_natural 
#Egg concentration experiment
load("./Rdata/egg_conc_exp.RData") #egg_conc_exp
Allee_data[["egg_conc_exp"]] <- egg_conc_exp 

###
