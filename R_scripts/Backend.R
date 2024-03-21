# Teo and Todd modified


# IMMEDIATE PRIORITIES:
#could try a high res model for the intial patch config and then convert back to a 1 m model.
# check time_k should  = dt instead of dx
# ATTEMPTING TO ADD STOKES DRIFT FOR X,Y VELOCITIES
# you made sperm sinking as a porotion/s. can this be done for advection and ascent?
# Slow decrease in cons mass after an initial increase. check
# the model can still give output greater than 1!? Can it still?
# run for edge, and weight based on position


# Run all models
# 1) Run bunya species
# 2) run uncertainity for all species
# - a hyac
#- a ten
# - a dig
# 3) run sens_ranges
# 4) int_dia. Run at 5 x 5?
# 5) intermediate sizes


# FIX SETPOINT CODE
# - LOTS OF MODEL BASED ON INITLA NO_TIMESTEP. tRYING USING END NO_TIMESTEP, AND CHANGING INITAL TIMESTEP FOR MODEL RUNS
# test setpoint at 500, end  = 1000.    Results = 500 = 0.073, 1000 = 0.2845.     Compare to running whole = 0.2840



# Secondary priorities:
# see python code for vectorised way of dispersal function (from 12 steps to navier). See dispersal5_diff_vec (code at bottom for outside)
# Possible could try the analystic solution (see gaussian model) but not sure if this allow cells without discretisation.
# ADD SOME TYPE OF BROWIAN MOTION OF SOMETHING LIKE stochastic_difference_x = np.sqrt(2 * D_L * dt) * np.random.normal()  https://pylag.readthedocs.io/en/latest/examples/lateral_adv_diff_analytic.html
# how can add asynchronous in? Move the model to a individual grids
# - one way could be to run a individual sperm map, transpose it x times for the patch and save these as lists. Then add together. Then runt he egg/embry map
# could use a extrapolated logarithmic model to predict to 2 hrs
# add in sperm sinking Kono or similar
#- the only way i can think to do this is to do a proportion from each advectarray (but this is tricky because slick is changing in the interim).  So not sure.
# could data.table speed up processing
# add in wall (topgraphical front), just set the length as not far? Sounds too easy


# Issues/fixed
#  - eggs should not be distributed in 1 m^3, so conversion to uL (/10^9) is too low i.e the egg conc will be much high.Pending
# - E0 will produce higher fert than usual because of (dilution) and likely not an issue anyway. Think this should be set to minimum because it functions
# to remove sperm but  ... (changed to 0.001 serm /ul = 1 sperm/mL. If E0 is set to 1, then it means Fe must be <0.01 to match other sfert curves.
# possible self fert effect - adjusted by using eggmap and running sequentially. # - self fer could be controlled using 1 egg and multiple sperm colonies.
# Then repeated for each position in the patch by moving the egg colony.
# I think I see what he is getting at now with gamete viability. Styan model uses tau.time to indicate when decay starts, but since each
# increment is 1 sec, then it will never kick in. So time needs to sum and then be fed back into the model.Add a gradual decay into the model like Teo
# need to check that the fertilisation increase at each time point is linear by comparing the Styan model for 1 hr in a group vs iteratively. Checked. This is ok
# Bundle breaking time data
# - he only takes every 59th second in the fert count vastly under estimating the fertilisation
# max fert shouldn't be 100% but should be adjusted for lab experiment maxes. Egg viability has been added in.
# could make tau a function of sperm concentration (although not sure we have data on that or if it is a thing)
# fecundity, sensense, polyps/branches, and infertility regions need to be added. Alvarez has this at 0.87 (SOM). (This has been adjusted)
# Added egg rad monte carlo in to model loop (confirmed working)
# Adjusted spe_speed to match mean Morita 2006.
# Added sperm swim monte carlo
# advection WC?

# egg radius might be shrunk. Maybe use Wallace

## Couple of weird things:
# I think maybe that the diffusion should be divided by the  number of receiving cells? But this would significantly reduce his mixing
# Maybe this is because of depth diffusion seems too strong
# is the nested loop acting instantaneous or iterative? I think you'd want each calculation independent of other calculations in the grid each
# for each time point. NOTE. Is is acting instantaneous - checked at one iteration
# I pulled fertcounter out from the last loop because it was updating for each cell rather than each time point??

# everything currently in metres (m) and seconds (s)

run_model <- function(scenario = "default", debug = F, store = F) {
  if(scenario == "default") {
    params <- settings_default
  } else if(scenario == "test") {
    params <- settings_test
  } else if(scenario == "custom") {
    params <- settings_custom
  } else {
    stop("Invalid scenario. Choose 'default', 'test', or provide custom settings.")
  }

  # load packages -----------------------------------------------------------

  # .rs.restartR()
  # library(foreach)
  # library(doParallel)
  library(abind)
  library(dplyr)
  library(lattice)
  library(gridExtra)
  library(ggplot2)
  library(RColorBrewer)
  library(truncnorm)
  library(profvis)
  #library(magick)
  #library(gganimate)
  library(Rcpp)

  source("./R_scripts/boundary_0_func.R")
  source("./R_scripts/config_func.R")
  source("./R_scripts/plot_fun1.R") # sperm
  source("./R_scripts/plot_fun2.R") # eggs
  source("./R_scripts/dispersal_func5_diff.R") # dispersal4
  #source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code



  # DEFAULTS. no_timestep =  3000, bundlebreak = 1050 , bundlebreak_sd = 289 , trans_dim = 46
  # a hya = colony_diam = 27, patch = 2022ahyac , depth = 4,  flow_rate = 0.10.  Total eggs 178315
  # tenuis = colony_diam = 17, patch = 2021tenuis, depth = 3,  flow_rate = 0.15
  # digit = colony_diam = 20, patch = 2021digit, depth = 3,  flow_rate = 0.11

  # TEST no_timestep =  500,bundlebreak =  400 , bundlebreak_sd = 100 , trans_dim = 15, sub_steps = 1



  # Input variables -------------------------------------------------------
  ## space/time parameters
  (no_timestep <- params$no_timestep) # (s) CHANGE THIS. min = 370 The total exposure time = 10800 (s) or 180 (min) or 3 hrs. Needs to meet min array
  time <- params$time # (s) intialise time
  polarbody <- params$polarbody # (s)  Teo uses seconds. Ive left it at 0 so things fert immediately after breakage
  (bundlebreak <- params$bundlebreak) # (s) 1350 From Wolstone 2004 with humulis removed. Teo  = 900 (s) = 15 (min). For 1 min iteration
  bundlebreak_sd <- params$bundlebreak_sd # (s)289 is mean indiv SD from Wolstone 2004 with humulis removed, see rscript

  ## patch parameters
  patch_type <- params$patch_type ## config types are patch_bot, patch_center, patch_outer, 2022ahyac, 2021tenuis, 2021digit, petal
  (int_col_spac <- params$int_col_spac) # (m) intercolonial spacing
  egg_viab <- params$egg_viab # (prop) egg viability.  Buccheri
  colony_diam <- params$colony_diam # (cm) colony diameter. A, hyac = 27. a. tenuis = 17.  digit = 20.
  polyp_den <- params$polyp_den # (polyps/cm^2) polyp density   #Wallace 1985 A. hya 182   #80.6 Alvarez SOM A. hyacinthus  #72-88 = 80 A digitifera in Doropoulos 2018
  fecun_zone <- params$fecun_zone # (prop) fecund of colony. Alvarez SOM for A. dig. A. hyac Palau = 0.70,
  bund_asc <- params$bund_asc # (m/s)  #125 s/m
  bund_asc_egg <- params$bund_asc_egg # (m/s)
  patch_density <- params$patch_density

  ## physical parameters
  Ks <- params$Ks # (m) is the characteristic roughness length-scale, using values from Baird and Atkinson’s (1997) flume experiments with coral skeletons. G - Think should be 0.27
  # but see Lentz paper that reports lower
  flow_rate <- params$flow_rate # (m s-1).  Probably close 300m/hr  or 5m/minute or 0.083 m/s s not far off FOR oNE tREE?
  # note that in the log shear equation u should be u(z) i,e flow at height above seabed
  van_K_c <- params$van_K_c # van kormans constant (
  long_const <- params$long_const # 0.20  (from Wistari n=4), 5.93  = NOTE THEE ARE DIFF CONSTANT NOT THE COEFFICENT. Elder. Teo =  0.15 based on shear flow ignoring
  trans_const <- params$trans_const # 0.135 (from Wistari, n = 4),  0.23 = Elder. Fischer 1979 = 0.15 based on experiments average
  vert_const <- params$vert_const # Fischer 1979 from Jobson and Sayre (1970)
  depth_m <- params$depth_m # (m). Depth at site. Also 6 m to match median Madin 2016/Muir 2015?
  # depth <- depth_m * down_scale
  (trans_dim <- params$trans_dim) # (m) These are the no. of columns in teh grid

  ## fertilisation kinetics parameters
  tb <- params$tb # 0.26 from estimation (dimensionless but 0 -> inf) polyspermy strength
  fe <- params$fe # 0.0037     0.0000952 / 0.00186  (G changed to match styans). Vogel had 0.01
  E0 <- params$E0 # 0.0364 eggs/ul = 1 egg/mL = from optimized . 10^6 eggs/m^3  :G changed)      #(eggmap5[i,j,k]) / 10^9
  spe_speed <- params$spe_speed # (mm/s)  #Teo used max 0.35. Webplot us 0.272 +/- 0.049 from Morita 2006 with egg?
  spe_speed_sd <- params$spe_speed_sd # (mms/s)
  egg_dia <- params$egg_dia #mm  #0.517   Waalce stahorns. 0.356 my Acropora tenuis data. SD = 0.011 mine dia = 0.0055 rad (added into loop rnorm)
  egg_dia_sd <- params$egg_dia_sd #
  eggperb <- params$eggperb # (eggs/bundle) from Wallace 1985. Add range 4-11. 3-10 for A. digit (Dorpopoulos 2018)
  speperb <- params$speperb # (sperm/bundle) From Furukawa 2020 fig 3 Acroora diviricara 2.1 x 10^6. sperm per bundle. Chris hyacinthus data 5.1 * 10^6.



  # calculations ------------------------------------------------------------

  ## size/time adjusted parameters, overrides orginals
  grid_size <- 1
  sub_steps <- 1 #
  down_scale <- 1 / grid_size
  (no_timestep <- no_timestep * sub_steps)
  polarbody <- polarbody * sub_steps
  (bundlebreak <- bundlebreak * sub_steps)
  bundlebreak_sd <- bundlebreak_sd * sub_steps
  (int_col_spac <- int_col_spac * down_scale)
  bund_asc <- bund_asc * down_scale / sub_steps
  bund_asc_egg <- bund_asc_egg / sub_steps # (m/s)
  flow_rate <- flow_rate * down_scale
  (trans_dim <- trans_dim * down_scale)


  ## secondary calculations
  #time/space
  cell_vol <- grid_size^3 # (m)
  cell_vol_uL <- cell_vol * 10^9 # (ul)
  down_scale <- 1 / grid_size
  dt <- 1 / sub_steps # timestep
  bund_probs <- dtruncnorm(1:no_timestep, a = 0, mean = bundlebreak, sd = bundlebreak_sd) # (prop)?
  ## patch params
  E0 <- E0 * cell_vol
  no_height <- sqrt(patch_density)   #4  #species config ignores these inputs ( ith ink)
  no_width <- sqrt(patch_density)   #5
  long_dim <- as.integer((flow_rate * no_timestep + 60) * down_scale) # (m) These are the no. of rows. Used these dimension before they are the max  eggs spread over 2 hrs.
  egg_rad <- egg_dia / 2 #
  egg_rad_sd <- egg_dia_sd / 2 #

  # create grids ------------------------------------------------------------
  dx <- 1 * grid_size # grid size
  dy <- 1 * grid_size # grid size
  dz <- 1 * grid_size
  z_dim <- (depth_m * down_scale) + 2 # dummy top and bottom layer (so 7 for 5 m depth)
  surface <- 2
  benthos <- z_dim - 1
  blank_spemap <- array(rep(0, long_dim * trans_dim * z_dim), c(long_dim, trans_dim, z_dim)) # original is 253 x 253 with 7 dimensions. Creats z number of matrices
  # fill the top and bottom layer with NA to serve as surface and bottom. Fill edges with NA

  # blank_spemap = boundary(blank_spemap = blank_spemap, z_dim = z_dim)
  blank_spemap <- boundary_0(x = blank_spemap, z_dim = z_dim)
  eggmap <- spemap <- blank_bo <- blank_unfert <- blank_phi_mono <- embmap <- blank_eggmap <- blank_spemap # empty array to fill fert data onto


  # configuration -----------------------------------------------------------
  # Use patch_cent when flow is <0.05

  #########################################
  #testing fine grid size for patch at 0.01m (not working yet)
  # grid_size_fine <- 0.01 # New grid size in meters
  # down_scale_fine <- 1 / grid_size_fine
  # long_dim_fine <- long_dim * 100
  # trans_dim_fine <- trans_dim * 100
  # spemap_fine <- array(0, dim = c(long_dim_fine, trans_dim_fine, z_dim))
  # eggmap_fine = spemap_fine
  # source("./R/config_fine.R")
  # config_lst <- config_fine(spemap_fine = spemap_fine, spemap = spemap, eggmap = eggmap, eggmap_fine = eggmap_fine, long_dim_fine = long_dim_fine,
  #                           trans_dim_fine = trans_dim_fine, z_dim = z_dim, int_col_spac = int_col_spac, config1 = patch_type, benthos = benthos,
  #                           no_height = no_height, no_width = no_width, patch_density = patch_density, colony_diam = colony_diam,
  #                           polyp_den = polyp_den, fecun_zone = fecun_zone) # type, length, width
  #
  # spemap <- config_lst$spemap
  # eggmap <- config_lst$eggmap
  # dim(spemap)
  # config_lst$small_mat_coarse

  ####################
  #original patch scale
  config_lst <- config(spemap = spemap, eggmap = eggmap, long_dim = long_dim, trans_dim = trans_dim, z_dim = z_dim, int_col_spac = int_col_spac,
                       config1 = patch_type, benthos = benthos, no_height = no_height, no_width = no_width, patch_density = patch_density,
                       colony_diam = colony_diam, fecun_zone = fecun_zone, polyp_den = polyp_den) # type, length, width
  spemap0 <- config_lst$spemap00
  eggmap0 <- config_lst$eggmap00
  ##############

  smallmat <- config_lst$smallmat
  sum(smallmat)
  spe_config <- config_lst$plot1
  #spe_config
  egg_config <- config_lst$plot2
  plot_fun1(spemap0, benthos)
  plot_fun2(eggmap0, benthos)

  # coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
  # levelplot(t(apply(smallmat[], 2, rev)),
  #   col.regions = coul, xlab = "Transverse",
  #   ylab = "Longitudinal", main = "Conc. (cells/m^3)"
  # )

  ##adults to bundles####this is now in config
  ## Create a patch and embed in the grid
  #col_area <- pi * (colony_diam / 2)^2 # (cm^2)
  #fecund <- col_area * polyp_den * fecun_zone # total bundles per colony
  #fecund <-  polyp_den * fecun_zone # total bundles per colony
  # fecund * 8 *0.08 *25  bundles * eggs * fertsuc * colonies
  #spemap0 <- spemap * fecund # bundles/adult
  #eggmap0 <- eggmap * fecund
  #p0 <- plot_fun1(spemap0, benthos)


  # Stokes drift#############
  # Input constants
  # wp = 4.5  #wave period (Fab). Time between waves
  # wave_height = 0.5   #(Fab)
  # z = -0.5  # Depth from surface (m) - negative because it's below the surface. Indepednat of model
  # k = 0.5  # Wave number (radians/m)  estimate
  # d = depth_m   #total depth to seabed
  #
  # #Derived constants
  # w = (2 * pi)/wp  # Angular frequency (radians/s)
  # lambda = 2 * pi / k  # Wavelength (m)
  # a = wave_height/2  # Wave amplitude (m)
  #
  # (u_stokes_drift_shallow = (w * k * a^2 * cosh(2 * k * (z + d))) / (2 * sinh(k * d)^2)) #Equaton 4 https://link.springer.com/article/10.1007/s10652-021-09811-8
  # direction_angle = 50  # degrees
  #
  # # Convert angle to radians for calculation
  # direction_angle_radians <- direction_angle * pi / 180
  #
  # # Calculating x and y components of the velocity
  # x_velocity <- u_stokes_drift_shallow * cos(direction_angle_radians)
  # y_velocity <- u_stokes_drift_shallow * sin(direction_angle_radians)


  ## Advection
  h_s <- rev(0:(z_dim - 3) + 0.5) / down_scale # create vec of off h_s using mid point from top to bot
  # Drag coefs
  drag_coef <- data.frame(depth_m = 1:10, dc = c(0.033, 0.019, 0.015, 0.012, 0.011, 0.010, 0.009, 0.009, 0.008, 0.008)) # from Lentz
  drag_coef$shear <- sqrt(drag_coef$dc * (flow_rate / down_scale)^2) # sometime ppl include a 1/2 but Lentz does. This is by subbing Eq(1) into u* in the text
  # shear_surf = 0.1 * flow_rate
  shear_surf <- drag_coef$shear[depth_m] # this extracts u* based on depth of benthos
  z0 <- 0.05 # hydrodynamic roughness scale from mean in Lentz Fig 9
  flow_rate_z <- ((shear_surf / van_K_c) * log(h_s / z0) * down_scale) # 0.05
  # advect_array = seq(1/flow_rate, no_timestep, 1/flow_rate) ##every advect_array step advect 1 cell longitudinal (up)
  advect_arrays <- lapply(flow_rate_z, function(fr) {
    as.integer(seq(1 / fr, no_timestep, 1 / fr))
  })

  ## Sperm sinking
  sperm_sink_cm_h <- 60 # from Kono
  sperm_sink_m_h <- sperm_sink_cm_h / 10
  sperm_sink_m_s <- sperm_sink_m_h / 3600


  ## Diffusion coefs
  D_L <- long_const * depth_m * shear_surf
  D_T <- trans_const * depth_m * shear_surf
  D_Z <- vert_const * depth_m * shear_surf

  ## Ascent array bund
  ascent_array <- as.integer(seq(1 / bund_asc, no_timestep, 1 / bund_asc)[1:(z_dim - 3)]) ## seconds to ascend each cell. ignore last move and 2 boundaries = 3.
  ascent_array_egg <- as.integer(seq(1 / bund_asc_egg, no_timestep, 1 / bund_asc_egg)[1:(z_dim - 3)])[!is.na(as.integer(seq(1 / bund_asc_egg, no_timestep, 1 / bund_asc_egg)[1:(z_dim - 3)]))]

  # model prep-----------------------------------------------------------------
  ## blank maps###
  spemap3_3 <- spemap3_2 <- spemap3_1 <- spemap3_0 <- spemap3 <- spemap1_3 <- spemap2 <- spemap2_1 <- spemap1_2 <- blank_spemap # create a temp mat to store popped bundles
  fertmap <- polymap <- embmap <- eggmap3_4 <- eggmap3_3 <- eggmap3_2 <- eggmap3_1 <- eggmap3_0 <- eggmap3 <- eggmap2 <- eggmap1_3 <- eggmap2_1 <- eggmap1_2 <- blank_eggmap # create a temp mat to store popped bundles

  ## Create 10 empty lists for sperm and egg temp store
  labels_temp <- c("1", "1_0", "1_1", "1_2", "1_3", "2", "2_0", "2_1", "2_2", "3", "3_0", "3_1", "3_2", "3_3", "3_4")
  for (i in 1:length(labels_temp)) {
    list_name <- paste0("spemap_temp", labels_temp[i])
    list_name1 <- paste0("eggmap_temp", labels_temp[i])
    assign(list_name, list())
    assign(list_name1, list())
  }

  embmap_store <- list() # set list space  1
  phi_mono_store <- list() # set list space  1
  polymap_store <- list() # set list space 1
  fert_store <- list() # set list space 1
  spe_count_mat_store <- list() # set list space
  egg_count_mat_store <- list() # set list space
  fertcounter <- rep(0, (no_timestep) / 1) # dummy array of total exposure time to be filled with fert values
  polycounter <- rep(0, (no_timestep) / 1) # dummy array of total exposure time to be filled with poly values

  spemap1 <- spemap0 # looping maps
  eggmap1 <- eggmap0 * egg_viab
  tot_spe_bund <- sum(spemap1, na.rm = T) # total bundles in spemap
  tot_egg_bund <- sum(eggmap1, na.rm = T) # total bundles in eggmap
  tot_sperm <- tot_spe_bund * speperb
  tot_eggs <- tot_egg_bund * eggperb
  bunds_sep_spe <- tot_spe_bund * bund_probs # number of bundles to break at each time point
  bunds_sep_egg <- tot_egg_bund * bund_probs # number of bundles to break at each time point
  bund_store <- numeric() # initialise bundle sep


  # spemap1 <- readRDS(file = "spemap_temp1_1_output_500.rds")[[500]]
  # eggmap1 <- readRDS(file = "eggmap_temp1_1_output_500.rds")[[500]]
  # spemap3_0 <- readRDS(file = "spemap_temp3_0_output_500.rds")[[500]]
  # eggmap3_0 <- readRDS(file = "eggmap_temp3_0_output_500.rds")[[500]]
  # time = 500


  # Start model -------------------------------------------------------------
  strt <- Sys.time()
  store_output <- store
  debug_output <- debug

  while (time <= no_timestep) # when time is less than 15 (min) run bundles
  {
    # bundle ascent (1's)
    spemap1_1 <- spemap1 # Loops from midway
    eggmap1_1 <- eggmap1

    # advection surface sperm bundles
    # spemap1_2a <- blank_spemap  #blanking ok because whole model is being replaced
    # for (i in 1:(long_dim))
    # {
    #   for (j in 1:(trans_dim))
    #   {
    #     for (k in 1:(z_dim))
    #     {
    #       if (is.na(spemap1_1[i,j,k])==FALSE )
    #       {
    #         if (time %in% advect_array)
    #           {
    #           spemap1_2a[(i-1),j,k] = spemap1_1[i,j,k]
    #       } else {
    #          spemap1_2a[i,j,k] = spemap1_1[i,j,k]
    #           }
    #       }
    #     }
    #   }
    # }
    ## Advection loop
    spemap1_2a <- blank_spemap # Blank map for the new model
    for (i in 2:(long_dim - 1)) {
      for (j in 2:(trans_dim - 1)) {
        for (k in 2:(z_dim - 1)) {
          if (!is.na(spemap1_1[i, j, k])) {
            if (time %in% advect_arrays[[k - 1]]) { # at time for the depth
              # new_i <- max(1, i - k) # Ensure new index is within bounds
              spemap1_2a[(i - 1), j, k] <- spemap1_1[i, j, k]
            } else {
              spemap1_2a[i, j, k] <- spemap1_1[i, j, k]
            }
          }
        }
      }
    }


    ## reset boundaries to zero
    spemap1_2a[1, , ] <- 0
    spemap1_2a[long_dim, , ] <- 0
    spemap1_2a[, 1, ] <- 0
    spemap1_2a[, trans_dim, ] <- 0
    spemap1_2a[, , 1] <- 0
    spemap1_2a[, , z_dim] <- 0


    ## 1_2 disperse
    # disperse (horizontal sperm bundles in 1_1) for each cell and add to previous cell to 1_2
    spemap1_2b <- blank_spemap # blanking ok because whole model is being replaced

    for (i in 1:(long_dim))
    {
      for (j in 1:(trans_dim))
      {
        for (k in 1:(z_dim))
        {
          if (is.na(spemap1_2a[i, j, k]) == FALSE) {
            dispdata <- dispersal5_diff(x = spemap1_2a, i = i, j = j, k = k, dx = dx, dy = dy, dz = dz, dt = dt, D_L = D_L, D_T = D_T, D_Z = D_Z, shear_surf = shear_surf, long_dim = long_dim, trans_dim = trans_dim, z_dim = z_dim, buoyant = TRUE) # this runs on the original to stop doubling up. disp2 true is horizontal any depth
            spemap1_2b[i, j, k] <- spemap1_2a[i, j, k] + dispdata # note this has to be adding (technically subtracting) to the previous map i.e previous count minus diffused
          }
        }
      }
    }

    # vectorised approach
    # dispersion_data1 <- dispersal5_diff_vec(spemap1_2a, dx = dx, dy = dy, dz = dz, dt = dt, D_L = D_L, D_T = D_T, D_Z = D_Z, shear_surf = shear_surf, buoyant = TRUE)
    # spemap1_2b <- spemap1_2a + dispersion_data1


    ## no-flux boundary condition for i = 1
    spemap1_2b[2, , ] <- spemap1_2b[2, , ] + spemap1_2b[1, , ]
    spemap1_2b[1, , ] <- 0

    # no-flux boundary condition for i = long_dim
    spemap1_2b[dim(spemap1_2b)[1] - 1, , ] <- spemap1_2b[dim(spemap1_2b)[1] - 1, , ] + spemap1_2b[dim(spemap1_2b)[1], , ]
    spemap1_2b[dim(spemap1_2b)[1], , ] <- 0

    # no-flux boundary condition for j = 1
    spemap1_2b[, 2, ] <- spemap1_2b[, 2, ] + spemap1_2b[, 1, ]
    spemap1_2b[, 1, ] <- 0

    # no-flux boundary condition for j = trans_dim
    spemap1_2b[, dim(spemap1_2b)[2] - 1, ] <- spemap1_2b[, dim(spemap1_2b)[2] - 1, ] + spemap1_2b[, dim(spemap1_2b)[2], ]
    spemap1_2b[, dim(spemap1_2b)[2], ] <- 0

    # no-flux boundary condition for k = 1
    spemap1_2b[, , 2] <- spemap1_2b[, , 2] + spemap1_2b[, , 1]
    spemap1_2b[, , 1] <- 0

    # no-flux boundary condition for k = z_dim
    spemap1_2b[, , dim(spemap1_2b)[3] - 1] <- spemap1_2b[, , dim(spemap1_2b)[3] - 1] + spemap1_2b[, , dim(spemap1_2b)[3]]
    spemap1_2b[, , dim(spemap1_2b)[3]] <- 0


    ## Conservation of mass correction
    spemap1_2b <- spemap1_2b * (sum(spemap1_1, na.rm = TRUE) / sum(spemap1_2b, na.rm = TRUE)) # Orig = new * Orig/new

    if (store_output) {
      spemap_temp1_1[[time]] <- spemap1_1 # pre-disperse
      spemap_temp1_2[[time]] <- spemap1_2b # afer disperse
    }

    # advection surface egg bundles
    # eggmap1_2a <- blank_spemap  #blanking ok because whole model is being replaced
    # for (i in 1:(long_dim))
    # {
    #   for (j in 1:(trans_dim))
    #   {
    #     for (k in 1:(z_dim))
    #     {
    #       if (is.na(eggmap1_1[i,j,k])==FALSE )
    #       {
    #         if (time %in% advect_array){
    #           eggmap1_2a[(i-1),j,k] = eggmap1_1[i,j,k]
    #         } else {
    #           eggmap1_2a[i,j,k] = eggmap1_1[i,j,k]
    #         }
    #       }
    #     }
    #   }
    # }

    eggmap1_2a <- blank_eggmap # Blank map for the new model
    for (i in 2:(long_dim - 1)) {
      for (j in 2:(trans_dim - 1)) {
        for (k in 2:(z_dim - 1)) {
          if (!is.na(eggmap1_1[i, j, k])) {
            if (time %in% advect_arrays[[k - 1]]) { # at time for the depth
              # new_i <- max(1, i - k) # Ensure new index is within bounds
              eggmap1_2a[(i - 1), j, k] <- eggmap1_1[i, j, k]
            } else {
              eggmap1_2a[i, j, k] <- eggmap1_1[i, j, k]
            }
          }
        }
      }
    }
    ## reset boundaries to zero
    eggmap1_2a[1, , ] <- 0
    eggmap1_2a[long_dim, , ] <- 0
    eggmap1_2a[, 1, ] <- 0
    eggmap1_2a[, trans_dim, ] <- 0
    eggmap1_2a[, , 1] <- 0
    eggmap1_2a[, , z_dim] <- 0



    ## disperse (horizontal sperm bundles in 1_1) for each cell and add to previous cell to 1_2
    eggmap1_2b <- blank_spemap # blanking ok because whole model is being replaced
    for (i in 1:(long_dim))
    {
      for (j in 1:(trans_dim))
      {
        for (k in 1:(z_dim))
        {
          if (is.na(eggmap1_2a[i, j, k]) == FALSE) {
            dispdata <- dispersal5_diff(x = eggmap1_2a, i = i, j = j, k = k, dx = dx, dy = dy, dz = dz, dt = dt, D_L = D_L, D_T = D_T, D_Z = D_Z, shear_surf = shear_surf, long_dim = long_dim, trans_dim = trans_dim, z_dim = z_dim, buoyant = TRUE) # this runs on the original to stop doubling up. disp2 true is horizontal any depth
            eggmap1_2b[i, j, k] <- eggmap1_2a[i, j, k] + dispdata # note this has to be adding (technically subtracting) to the previous map i.e previous count minus diffused
          }
        }
      }
    }

    # vectorised approach
    # dispersion_data2 <- dispersal5_diff_vec(eggmap1_2a, dx = dx, dy = dy, dz = dz, dt = dt, D_L = D_L, D_T = D_T, D_Z = D_Z, shear_surf = shear_surf, buoyant = TRUE)
    # eggmap1_2b <- eggmap1_2a + dispersion_data2

    # # no-flux boundary condition for i = 1
    # if (i == 1 && eggmap1_2b[i, j, k] != 0) {
    #   eggmap1_2b[2, j, k] <- eggmap1_2b[2, j, k] + eggmap1_2b[i, j, k]
    #   eggmap1_2b[i, j, k] <- 0
    # }
    # # no-flux boundary condition for i = long_dim
    # if (i == long_dim && eggmap1_2b[i, j, k] != 0) {
    #   eggmap1_2b[i-1, j, k] <- eggmap1_2b[i-1, j, k] + eggmap1_2b[i, j, k]
    #   eggmap1_2b[i, j, k] <- 0
    # }
    # # no-flux boundary condition for j = 1
    # if (j == 1 && eggmap1_2b[i, j, k] != 0) {
    #   eggmap1_2b[i, 2, k] <- eggmap1_2b[i, 2, k] + eggmap1_2b[i, j, k]
    #   eggmap1_2b[i, j, k] <- 0
    # }
    # # no-flux boundary condition for j = trans_dim
    # if (j == trans_dim && eggmap1_2b[i, j, k] != 0) {
    #   eggmap1_2b[i, j-1, k] <- eggmap1_2b[i, j-1, k] + eggmap1_2b[i, j, k]
    #   eggmap1_2b[i, j, k] <- 0
    # }
    # # no-flux boundary condition for k = 1
    # if (k == 1 && eggmap1_2b[i, j, k] != 0) {
    #   eggmap1_2b[i, j, 2] <- eggmap1_2b[i, j, 2] + eggmap1_2b[i, j, k]
    #   eggmap1_2b[i, j, k] <- 0
    # }
    # # no-flux boundary condition for k = z_dim
    # if (k == z_dim && eggmap1_2b[i, j, k] != 0) {
    #   eggmap1_2b[i, j, k-1] <- eggmap1_2b[i, j, k-1] + eggmap1_2b[i, j, k]
    #   eggmap1_2b[i, j, k] <- 0
    # }

    ## no-flux boundary condition for i = 1
    eggmap1_2b[2, , ] <- eggmap1_2b[2, , ] + eggmap1_2b[1, , ]
    eggmap1_2b[1, , ] <- 0

    # no-flux boundary condition for i = long_dim
    eggmap1_2b[dim(eggmap1_2b)[1] - 1, , ] <- eggmap1_2b[dim(eggmap1_2b)[1] - 1, , ] + eggmap1_2b[dim(eggmap1_2b)[1], , ]
    eggmap1_2b[dim(eggmap1_2b)[1], , ] <- 0

    # no-flux boundary condition for j = 1
    eggmap1_2b[, 2, ] <- eggmap1_2b[, 2, ] + eggmap1_2b[, 1, ]
    eggmap1_2b[, 1, ] <- 0

    # no-flux boundary condition for j = trans_dim
    eggmap1_2b[, dim(eggmap1_2b)[2] - 1, ] <- eggmap1_2b[, dim(eggmap1_2b)[2] - 1, ] + eggmap1_2b[, dim(eggmap1_2b)[2], ]
    eggmap1_2b[, dim(eggmap1_2b)[2], ] <- 0

    # no-flux boundary condition for k = 1
    eggmap1_2b[, , 2] <- eggmap1_2b[, , 2] + eggmap1_2b[, , 1]
    eggmap1_2b[, , 1] <- 0

    # no-flux boundary condition for k = z_dim
    eggmap1_2b[, , dim(eggmap1_2b)[3] - 1] <- eggmap1_2b[, , dim(eggmap1_2b)[3] - 1] + eggmap1_2b[, , dim(eggmap1_2b)[3]]
    eggmap1_2b[, , dim(eggmap1_2b)[3]] <- 0


    ## Conservation of mass correction
    eggmap1_2b <- eggmap1_2b * (sum(spemap1_1, na.rm = TRUE) / sum(eggmap1_2b, na.rm = TRUE))


    if (store_output) {
      eggmap_temp1_1[[time]] <- eggmap1_1 # before disperse
      eggmap_temp1_2[[time]] <- eggmap1_2b # bundles
    }

    ## 1_3 ascent bund####
    ## ascent water column sperm bundles at various time points (1_3)
    spemap1_3 <- blank_spemap
    if (time %in% ascent_array) { # Move up at ascent array time points
      # Moving everything up except for k boundaries
      for (i in 1:long_dim) {
        for (j in 1:trans_dim) {
          for (k in 2:(z_dim - 1)) {
            spemap1_3[i, j, (k - 1)] <- spemap1_2b[i, j, k]
          }
        }
      }
      # no-flux boundary condition at k=2, by moving k = 1 to 2
      spemap1_3[, , 2] <- spemap1_3[, , 1] + spemap1_3[, , 2]
      # Reset boundaries to 0
      spemap1_3[, , 1] <- 0
      spemap1_3[, , z_dim] <- 0
    } else { # If not in ascent_array time, just copy the values
      spemap1_3[1:long_dim, 1:trans_dim, 2:(z_dim - 1)] <- spemap1_2b[1:long_dim, 1:trans_dim, 2:(z_dim - 1)]
    }

    if (store_output) {
      spemap_temp1_3[[time]] <- spemap1_3 # ascent
      # spemap_temp3[[125]][,,6]
    }

    ## ascent water column egg bundles at various time points (1_3)
    eggmap1_3 <- blank_eggmap
    if (time %in% ascent_array) { # Move up at ascent array time points
      # Moving everything up except for k boundaries
      for (i in 1:long_dim) {
        for (j in 1:trans_dim) {
          for (k in 2:(z_dim - 1)) {
            eggmap1_3[i, j, (k - 1)] <- eggmap1_2b[i, j, k]
          }
        }
      }
      # no-flux boundary condition at k=2, by moving k = 1 to 2
      eggmap1_3[, , 2] <- eggmap1_3[, , 1] + eggmap1_3[, , 2]
      # Reset boundaries to 0
      eggmap1_3[, , 1] <- 0
      eggmap1_3[, , z_dim] <- 0
    } else { # If not in ascent_array time, just copy the values
      eggmap1_3[1:long_dim, 1:trans_dim, 2:(z_dim - 1)] <- eggmap1_2b[1:long_dim, 1:trans_dim, 2:(z_dim - 1)]
    }

    if (store_output) {
      eggmap_temp1_3[[time]] <- eggmap1_3 # ascent
    }

    spemap2_1 <- spemap1_3
    eggmap2_1 <- eggmap1_3

    if (store_output) {
      spemap_temp2_1[[time]] <- spemap2_1 # advection
    }

    if (store_output) {
      eggmap_temp2_1[[time]] <- eggmap2_1
    }

    ##dissociated separated bundles####
    ## Prepare maps for bundle dissociation. Separate popped bundles to map 3s
    spemap3 <- blank_spemap # blanks because emptied to 3_1 at each run
    tot_bund_spemap2_1 <- sum(spemap2_1, na.rm = T) # total bundles at no_timestep
    # bundle_frac_spe = bunds_sep_spe[time] / tot_bund_spemap2_1  #fraction of bundles in each cell as prop of recent
    bundle_frac_spe <- bunds_sep_spe[time] / tot_spe_bund # fraction of bundles in each cell  as prop of original

    remove_spe <- spemap2_1 * bundle_frac_spe # probs * no per cell
    bund_store[time] <- sum(remove_spe, na.rm = T) * speperb # store no of bundles removed for each timepoint
    # Subtract the calculated value from each cell in spemap2_1 and add to spemap3
    spemap2_1 <- spemap2_1 - remove_spe #
    spemap3 <- spemap3 + remove_spe

    spemap1 <- spemap2_1 #  LOOPS back to original and overrides
    spemap2_2 <- spemap2_1 # no purpose, can adjust later

    if (store_output) {
      spemap_temp2_2[[time]] <- spemap2_2 # sperm bundles after popped removed
      spemap_temp3[[time]] <- spemap3 # popped bundles
    }

    ## Prepare maps for bundle dissociation. Move popped bundles to map 3s
    eggmap3 <- blank_spemap # blanks because emptied to 3_1 at each run
    tot_bund_eggmap2_1 <- sum(eggmap2_1, na.rm = T) # total bundles at no_timestep
    bundle_frac_egg <- bunds_sep_egg[time] / tot_bund_eggmap2_1 # fraction of bundles in each cell
    remove_egg <- eggmap2_1 * bundle_frac_egg # probs * no per cell
    # Subtract the calculated value from each cell in eggmap2_1 and add to eggmap3
    eggmap2_1 <- eggmap2_1 - remove_egg # ensure no negative values
    eggmap3 <- eggmap3 + remove_egg

    eggmap1 <- eggmap2_1 #  LOOPS back to original and overrides
    eggmap2_2 <- eggmap2_1 # no purpose, can adjust later

    if (store_output) {
      eggmap_temp2_2[[time]] <- eggmap2_2 # egg bundles after popped removed
      eggmap_temp3[[time]] <- eggmap3 # popped bundles
    }

    ##free + existing gametes####
    ## dissociated sperm bundles
    # takes new bundles(3), breaks them, adds them to previous count (3_0)
    spemap3_1 <- spemap3_0 # initialize with blank, LOOPS to here from end code.
    spe_count_mat <- spemap3 * speperb
    spemap3_1 <- spemap3_0 + spe_count_mat # existing sperm + popped sperm

    if (store_output) {
      spemap_temp3_0[[time]] <- spemap3_0 # looped mat (existing sperm)
      spe_count_mat_store[[time]] <- spe_count_mat
      spemap_temp3_1[[time]] <- spemap3_1 # combined sperm (popped + existing
    }

    ## dissociation egg bundles
    eggmap3_1 <- eggmap3_0
    egg_count_mat <- eggmap3 * eggperb # 2 dim mat
    eggmap3_1 <- eggmap3_0 + egg_count_mat

    if (store_output) {
      eggmap_temp3_0[[time]] <- eggmap3_0
      eggmap_temp3_1[[time]] <- eggmap3_1
    }


    ## advect
    spemap3_2a <- blank_spemap # Blank map for the new model
    for (i in 2:(long_dim - 1)) {
      for (j in 2:(trans_dim - 1)) {
        for (k in 2:(z_dim - 1)) {
          if (!is.na(spemap3_1[i, j, k])) {
            if (time %in% advect_arrays[[k - 1]]) { # at time for the depth
              # new_i <- max(1, i - k) # Ensure new index is within bounds
              spemap3_2a[(i - 1), j, k] <- spemap3_1[i, j, k]
            } else {
              spemap3_2a[i, j, k] <- spemap3_1[i, j, k]
            }
          }
        }
      }
    }
    ## reset boundaries to zero
    spemap3_2a[1, , ] <- 0
    spemap3_2a[long_dim, , ] <- 0
    spemap3_2a[, 1, ] <- 0
    spemap3_2a[, trans_dim, ] <- 0
    spemap3_2a[, , 1] <- 0
    spemap3_2a[, , z_dim] <- 0



    ## 3_2 dispersed gametes####
    spemap3_2b <- blank_spemap
    for (i in 1:(long_dim))
    {
      for (j in 1:(trans_dim))
      {
        for (k in 1:(z_dim)) # don't want top layer i.e [,,2]
        {
          if (is.na(spemap3_2a[i, j, k]) == FALSE) {
            dispdata <- dispersal5_diff(x = spemap3_2a, i = i, j = j, k = k, dx = dx, dy = dy, dz = dz, dt = dt, D_L = D_L, D_T = D_T, D_Z = D_Z, shear_surf = shear_surf, long_dim = long_dim, trans_dim = trans_dim, z_dim = z_dim, buoyant = FALSE) # this runs on the original to stop doubling up. disp2 true is horizontal any depth
            spemap3_2b[i, j, k] <- spemap3_2a[i, j, k] + dispdata # note this has to be adding (technically subtracting) to the previous map i.e previous count minus diffused
          }
        }
      }
    }

    # vectorised approach
    # dispersion_data3 <- dispersal5_diff_vec(spemap3_2a, dx = dx, dy = dy, dz = dz, dt = dt, D_L = D_L, D_T = D_T, D_Z = D_Z, shear_surf = shear_surf, buoyant = TRUE)
    # spemap3_2b <- spemap3_2a + dispersion_data3


    ## no-flux boundary condition for i = 1
    spemap3_2b[2, , ] <- spemap3_2b[2, , ] + spemap3_2b[1, , ]
    spemap3_2b[1, , ] <- 0

    # no-flux boundary condition for i = long_dim
    spemap3_2b[dim(spemap3_2b)[1] - 1, , ] <- spemap3_2b[dim(spemap3_2b)[1] - 1, , ] + spemap3_2b[dim(spemap3_2b)[1], , ]
    spemap3_2b[dim(spemap3_2b)[1], , ] <- 0

    # no-flux boundary condition for j = 1
    spemap3_2b[, 2, ] <- spemap3_2b[, 2, ] + spemap3_2b[, 1, ]
    spemap3_2b[, 1, ] <- 0

    # no-flux boundary condition for j = trans_dim
    spemap3_2b[, dim(spemap3_2b)[2] - 1, ] <- spemap3_2b[, dim(spemap3_2b)[2] - 1, ] + spemap3_2b[, dim(spemap3_2b)[2], ]
    spemap3_2b[, dim(spemap3_2b)[2], ] <- 0

    # no-flux boundary condition for k = 1
    spemap3_2b[, , 2] <- spemap3_2b[, , 2] + spemap3_2b[, , 1]
    spemap3_2b[, , 1] <- 0

    # no-flux boundary condition for k = z_dim
    spemap3_2b[, , dim(spemap3_2b)[3] - 1] <- spemap3_2b[, , dim(spemap3_2b)[3] - 1] + spemap3_2b[, , dim(spemap3_2b)[3]]
    spemap3_2b[, , dim(spemap3_2b)[3]] <- 0


    if (store_output) {
      spemap_temp3_2[[time]] <- spemap3_2b
    }

    spemap3_3 <- spemap3_2b

    ## sinking sperm

    moving_sperm <- sperm_sink_m_s * spemap3_2b[, , 2]
    spemap3_3[, , 3] <- spemap3_2b[, , 3] + moving_sperm
    spemap3_3[, , 2] <- spemap3_2b[, , 2] - moving_sperm


    ## advection  eggs
    eggmap3_2a <- blank_eggmap # Blank map for the new model
    for (i in 2:(long_dim - 1)) {
      for (j in 2:(trans_dim - 1)) {
        for (k in 2:(z_dim - 1)) {
          if (!is.na(eggmap3_1[i, j, k])) {
            if (time %in% advect_arrays[[k - 1]]) { # at time for the depth
              # new_i <- max(1, i - k) # Ensure new index is within bounds
              eggmap3_2a[(i - 1), j, k] <- eggmap3_1[i, j, k]
            } else {
              eggmap3_2a[i, j, k] <- eggmap3_1[i, j, k]
            }
          }
        }
      }
    }
    # reset boundaries to zero
    eggmap3_2a[1, , ] <- 0
    eggmap3_2a[long_dim, , ] <- 0
    eggmap3_2a[, 1, ] <- 0
    eggmap3_2a[, trans_dim, ] <- 0
    eggmap3_2a[, , 1] <- 0
    eggmap3_2a[, , z_dim] <- 0


    ## dispersion
    eggmap3_2b <- blank_eggmap
    for (i in 1:(long_dim))
    {
      for (j in 1:(trans_dim))
      {
        for (k in 1:(z_dim)) # don't want top layer i.e [,,2]
        {
          if (is.na(eggmap3_2a[i, j, k]) == FALSE) {
            dispdata <- dispersal5_diff(x = eggmap3_2a, i = i, j = j, k = k, dx = dx, dy = dy, dz = dz, dt = dt, D_L = D_L, D_T = D_T, D_Z = D_Z, shear_surf = shear_surf, long_dim = long_dim, trans_dim = trans_dim, z_dim = z_dim, buoyant = TRUE) # this runs on the original to stop doubling up. disp2 true is horizontal any depth
            eggmap3_2b[i, j, k] <- eggmap3_2a[i, j, k] + dispdata # note this has to be adding (technically subtracting) to the previous map i.e previous count minus diffused
          }
        }
      }
    }

    # vectorised approach
    # dispersion_data4 <- dispersal5_diff_vec(eggmap3_2a, dx = dx, dy = dy, dz = dz, dt = dt, D_L = D_L, D_T = D_T, D_Z = D_Z, shear_surf = shear_surf, buoyant = TRUE)
    # eggmap3_2b <- eggmap3_2a + dispersion_data4

    # # no-flux boundary condition for i = 1
    # if (i == 1 && eggmap3_2b[i, j, k] != 0) {
    #   eggmap3_2b[2, j, k] <- eggmap3_2b[2, j, k] + eggmap3_2b[i, j, k]
    #   eggmap3_2b[i, j, k] <- 0
    # }
    # # no-flux boundary condition for i = long_dim
    # if (i == long_dim && eggmap3_2b[i, j, k] != 0) {
    #   eggmap3_2b[i-1, j, k] <- eggmap3_2b[i-1, j, k] + eggmap3_2b[i, j, k]
    #   eggmap3_2b[i, j, k] <- 0
    # }
    # # no-flux boundary condition for j = 1
    # if (j == 1 && eggmap3_2b[i, j, k] != 0) {
    #   eggmap3_2b[i, 2, k] <- eggmap3_2b[i, 2, k] + eggmap3_2b[i, j, k]
    #   eggmap3_2b[i, j, k] <- 0
    # }
    # # no-flux boundary condition for j = trans_dim
    # if (j == trans_dim && eggmap3_2b[i, j, k] != 0) {
    #   eggmap3_2b[i, j-1, k] <- eggmap3_2b[i, j-1, k] + eggmap3_2b[i, j, k]
    #   eggmap3_2b[i, j, k] <- 0
    # }
    # # no-flux boundary condition for k = 1
    # if (k == 1 && eggmap3_2b[i, j, k] != 0) {
    #   eggmap3_2b[i, j, 2] <- eggmap3_2b[i, j, 2] + eggmap3_2b[i, j, k]
    #   eggmap3_2b[i, j, k] <- 0
    # }
    # # no-flux boundary condition for k = z_dim
    # if (k == z_dim && eggmap3_2b[i, j, k] != 0) {
    #   eggmap3_2b[i, j, k-1] <- eggmap3_2b[i, j, k-1] + eggmap3_2b[i, j, k]
    #   eggmap3_2b[i, j, k] <- 0
    # }

    ## no-flux boundary condition for i = 1
    eggmap3_2b[2, , ] <- eggmap3_2b[2, , ] + eggmap3_2b[1, , ]
    eggmap3_2b[1, , ] <- 0

    # no-flux boundary condition for i = long_dim
    eggmap3_2b[dim(eggmap3_2b)[1] - 1, , ] <- eggmap3_2b[dim(eggmap3_2b)[1] - 1, , ] + eggmap3_2b[dim(eggmap3_2b)[1], , ]
    eggmap3_2b[dim(eggmap3_2b)[1], , ] <- 0

    # no-flux boundary condition for j = 1
    eggmap3_2b[, 2, ] <- eggmap3_2b[, 2, ] + eggmap3_2b[, 1, ]
    eggmap3_2b[, 1, ] <- 0

    # no-flux boundary condition for j = trans_dim
    eggmap3_2b[, dim(eggmap3_2b)[2] - 1, ] <- eggmap3_2b[, dim(eggmap3_2b)[2] - 1, ] + eggmap3_2b[, dim(eggmap3_2b)[2], ]
    eggmap3_2b[, dim(eggmap3_2b)[2], ] <- 0

    # no-flux boundary condition for k = 1
    eggmap3_2b[, , 2] <- eggmap3_2b[, , 2] + eggmap3_2b[, , 1]
    eggmap3_2b[, , 1] <- 0

    # no-flux boundary condition for k = z_dim
    eggmap3_2b[, , dim(eggmap3_2b)[3] - 1] <- eggmap3_2b[, , dim(eggmap3_2b)[3] - 1] + eggmap3_2b[, , dim(eggmap3_2b)[3]]
    eggmap3_2b[, , dim(eggmap3_2b)[3]] <- 0


    if (store_output) {
      eggmap_temp3_2[[time]] <- eggmap3_2b
    }


    ## ascent eggs only####
    ##ascent water column egg bundles at various time points (1_3)
    eggmap3_3 <- blank_eggmap
    if (time %in% ascent_array_egg) { # Move up at ascent array time points
      # Moving everything up except for k boundaries
      for (i in 1:long_dim) {
        for (j in 1:trans_dim) {
          for (k in 2:(z_dim - 1)) {
            eggmap3_3[i, j, (k - 1)] <- eggmap3_2b[i, j, k]
          }
        }
      }
      # no-flux boundary condition at k=2, by moving k = 1 to 2
      eggmap3_3[, , 2] <- eggmap3_3[, , 1] + eggmap3_3[, , 2]
      # Reset boundaries to 0
      eggmap3_3[, , 1] <- 0
      eggmap3_3[, , z_dim] <- 0
    } else { # If not in ascent_array time, just copy the values
      eggmap3_3[1:long_dim, 1:trans_dim, 2:(z_dim - 1)] <- eggmap3_2b[1:long_dim, 1:trans_dim, 2:(z_dim - 1)]
    }


    if (store_output) {
      spemap_temp3_3[[time]] <- spemap3_3
    }

    if (store_output) {
      eggmap_temp3_3[[time]] <- eggmap3_3 # advected eggs
    }

    ## 3_4 fert####
    ## fertilisation section
    eggmap3_4 <- eggmap3_3
    if (time >= polarbody) {
      for (i in 2:(long_dim - 1)) # all rows
      {
        for (j in 2:(trans_dim - 1)) # all columns
        {
          for (k in 2:(z_dim - 1)) # all columns
          {
            if ((spemap3_3[i, j, k] > 0) == TRUE && (eggmap3_3[i, j, k] > 0) == TRUE) # if the sperm cell value is > 0 and the corresponding egg cell is > 0 . Unsure if '&&' is best here
            {
              S0_m3 <- spemap3_3[i, j, k] # sperm/m^3
              S0 <- S0_m3 / cell_vol_uL # (per uL)
              Bo <- pi * (rnorm(1, egg_rad, egg_rad_sd)^2) * rnorm(1, spe_speed, spe_speed_sd) # beta0 = egg-sperm collision rate constant (mm^3/s; Styan 1998 p291); merulina ampliata egg diameter ~0.339 mm so
              time_k <- dx # time (s) is the time the eggs are exposed to sperm
              xx <- fe * S0 / E0 * (1 - exp(-Bo * E0 * time_k)) # average number of potential fertilizers per egg
              bb <- fe * S0 / E0 * (1 - exp(-Bo * E0 * tb)) # mean number of extra fertilizing sperm that contact an egg in a time period tb in the population

              # fertilisation
              phi_poly <- (1 - exp(-xx) - xx * exp(-xx)) * (1 - exp(-bb)) # prop poly
              phi_mono <- 1 - exp(-xx) - phi_poly # (prop) fert.
              blank_phi_mono[i, j, k] <- phi_mono # store prop fert over space
              phi_mono_store[[time]] <- blank_phi_mono
              fert_no <- eggmap3_3[i, j, k] * phi_mono # instant numbers of fert. eggs per cell (counts but not integers).
              fertmap[i, j, k] <- fert_no
              fert_store[[time]] <- fertmap
              embmap[i, j, k] <- embmap[i, j, k] + fert_no # no of fert eggs/embryos (counts - not integers). Cumulative. Point of fert
              embmap_store[[time]] <- embmap

              # polyspermy
              poly <- eggmap3_3[i, j, k] * phi_poly ## instant poly eggs per cell (counts but not integers).
              polymap[i, j, k] <- polymap[i, j, k] + poly # add poly to previous poly  eggs (counts but not integers)
              polymap_store[[time]] <- polymap

              # remaining eggs
              eggmap3_4[i, j, k] <- eggmap3_3[i, j, k] - fert_no - poly # remove embryos from egg pool
              # eggmap3_4[,,2]  <- eggmap3_3[,,2] - embmap[,,2]  - polymap[,,2] #remove embryos from egg pool
              eggmap_temp3_4[[time]] <- eggmap3_4 # egg
            }
          }
        }
      }
    }
    ##output
    fertcounter[time] <- fertcounter[time] + (sum(embmap[, , ], na.rm = T) / tot_eggs) # % fert per time point. #the inside loop will keep on adding each cell until time changes
    spemap3_0 <- spemap3_3 # this to loop back to
    eggmap3_0 <- eggmap3_4
    print(c("time_s", time))
    print(c("fert_success", fertcounter[time])) # cumulative fert%
    if (debug_output) {
      print(c("poly", sum(polymap[, , ], na.rm = T)))
      print(c("embryos", sum(embmap[, , ], na.rm = T)))
      # print(c('bund_eggs1_1', sum(eggmap1_1[,,]*eggperb, na.rm = T))) #inital
      # print(c('bund_eggs1_2', sum(eggmap1_2[,,]*eggperb, na.rm = T))) #diserpse bundle
      # print(c('bund_eggs1_3', sum(eggmap1_3[,,]*eggperb, na.rm = T)))
      # print(c('bund_eggs2_1', sum(eggmap2_1[,,]*eggperb, na.rm = T)))
      # print(c('bund_eggs2_2', sum(eggmap2_2[,,]*eggperb, na.rm = T)))
      # print(c('bund_probs', bund_probs[time]))
      # print(c('bund_eggs3', sum(eggmap3[,,] * eggperb, na.rm = T)))
      # print(c('eggs3_1', sum(eggmap3_1[,,], na.rm = T)))
      # print(c('eggs3_2', sum(eggmap3_2[,,], na.rm = T)))
      # print(c('eggs3_4', sum(eggmap3_4[,,], na.rm = T)))
      # print(c('cons_mass_egg', sum(eggmap2_2[,,]*eggperb, na.rm = T) + sum(eggmap3_4[,,], na.rm = T) + sum(embmap[,,], na.rm = T)+ sum(polymap[,,], na.rm = T)))
      # if (time %% 100 == 0) {
      #   # Save the list to a file
      #   saveRDS(spemap_temp1_1, file = paste0('spemap_temp1_1_', "output_", time, ".rds"))
      #   saveRDS(eggmap_temp1_1, file = paste0('eggmap_temp1_1_',"output_", time, ".rds"))
      #   saveRDS(spemap_temp3_0, file = paste0('spemap_temp3_0_',"output_", time, ".rds"))
      #   saveRDS(eggmap_temp3_0, file = paste0('eggmap_temp3_0_',"output_", time, ".rds"))
      # }

      print(c("bund_spe1_1", format(sum(spemap1_1[, , ] * speperb, na.rm = T), scientific = TRUE)))
      print(c("bund_spe1_2", format(sum(spemap1_2[, , ] * speperb, na.rm = T), scientific = TRUE)))
      print(c("bund_spe1_3", format(sum(spemap1_3[, , ] * speperb, na.rm = T), scientific = TRUE)))
      print(c("bund_spe2_1", format(sum(spemap2_1[, , ] * speperb, na.rm = T), scientific = TRUE)))
      print(c("bund_spe2_2", format(sum(spemap2_2[, , ] * speperb, na.rm = T), scientific = TRUE)))
      print(c("bund_probs", format(bund_probs[time], scientific = TRUE)))
      print(c("bund_spe3", format(sum(spemap3[, , ] * speperb, na.rm = T), scientific = TRUE)))
      print(c("spe3_1", format(sum(spemap3_1[, , ], na.rm = T), scientific = TRUE)))
      print(c("spe3_2", format(sum(spemap3_2[, , ], na.rm = T), scientific = TRUE)))
      print(c("cons_mass_spe", format(sum(spemap2_2[, , ] * speperb, na.rm = T) + sum(spemap3_3[, , ], na.rm = T) + sum(embmap[, , ], na.rm = T) + sum(polymap[, , ], na.rm = T), scientific = TRUE)))
    }

    time <- time + 1
    # print(eggmap5)
    flush.console()
  }

  #
  end_time_opt <- Sys.time()
  time_diff <- end_time_opt - strt
  time_diff

}
# 7 Model end - run above ---------------------------------------------------------------
# save.image(file = file.path("./Rdata/500t.Rdata"))
#load("./Rdata/500t.Rdata")



