## Creating a new folder inside the 'sweeps' directory in the current working directory
dir_name <- "2024_01_23_ah_int1_size9_16_25"
dir.create(file.path(getwd(),  dir_name), showWarnings = T)
new_folder_path <- file.path(getwd(), dir_name)
library(foreach)
library(doParallel)

## Set up the parallel backend using doParallel and the number of cores
num_cores <- 3
registerDoParallel(cores = num_cores)
vec_x <- c(3,4,5)
foreach(i = 1:length(vec_x), .combine = 'c') %dopar% {

  #####PASTE FROM HERE#####################

  # load packages -----------------------------------------------------------
  library(abind)
  library(dplyr)
  library(lattice)
  library(gridExtra)
  library(ggplot2)
  library(RColorBrewer)
  library(truncnorm)
  library(profvis)
  library(magick)
  library(gganimate)
  library(Rcpp)
  source("./R/boundary_0_func.R")
  source("./R/config func.R")
  source("./R/plot_fun1.R")
  source("./R/plot_fun2.R")
  source("./R/dispersal_func5_diff.R")

  # Input variables -------------------------------------------------------

  ## space/time parameter
  grid_size <- 1
  cell_vol <- grid_size^3
  cell_vol_uL <- cell_vol * 10^9
  down_scale <- 1 / grid_size
  sub_steps <- 1
  dt <- 1 / sub_steps
  (no_timestep <- 500 * sub_steps)
  time <- 1
  polarbody <- 0 * sub_steps
  (bundlebreak <- 400 * sub_steps)
  bundlebreak_sd <- 100 * sub_steps
  bund_probs <- dtruncnorm(1:no_timestep, a = 0, mean = bundlebreak, sd = bundlebreak_sd)

  ## patch parameters
  patch_type <- "patch_bot"
  (int_col_spac <- vec_x[i] * down_scale)
  egg_viab <- 0.92
  colony_diam <- 27
  polyp_den <- 80.6
  fecun_zone <- 0.7
  bund_asc <- 0.008 * down_scale / sub_steps
  bund_asc_egg <- 0.0027 / sub_steps
  patch_density <- 4
  no_height <- sqrt(patch_density)
  no_width <- sqrt(patch_density)

  ## physical parameters
  Ks <- 0.27
  flow_rate <- 0.10 * down_scale
  van_K_c <- 0.41
  long_const <- 0.20
  trans_const <- 0.135
  vert_const <- 0.067
  depth_m <- 4
  (trans_dim <- 10 * down_scale)
  long_dim <- as.integer((flow_rate * no_timestep + 60) * down_scale)

  ## fertilisation kinetics parameters
  tb <- 0.26
  fe <- 0.0037
  E0 <- 0.0364 * cell_vol
  spe_speed <- 0.35
  spe_speed_sd <- 0.049
  egg_rad <- 0.517 / 2
  egg_rad_sd <- 0.005 / 2
  eggperb <- 6
  speperb <- 2.1 * 10^6

  # create grids ------------------------------------------------------------
  dx <- 1 * grid_size
  dy <- 1 * grid_size
  dz <- 1 * grid_size
  z_dim <- (depth_m * down_scale) + 2
  surface <- 2
  benthos <- z_dim - 1
  blank_spemap <- array(rep(0, long_dim * trans_dim * z_dim), c(long_dim, trans_dim, z_dim))
  blank_spemap <- boundary_0(x = blank_spemap, z_dim = z_dim)
  eggmap <- spemap <- blank_bo <- blank_unfert <- blank_phi_mono <- embmap <- blank_eggmap <- blank_spemap

  ## Create a patch and embed in the grid
  col_area <- pi * (colony_diam / 2)^2
  fecund <- col_area * polyp_den * fecun_zone

  # configuration -----------------------------------------------------------
  config_lst <- config(spemap = spemap, eggmap = eggmap, long_dim = long_dim, trans_dim = trans_dim, z_dim = z_dim, int_col_spac = int_col_spac, config1 = patch_type, benthos = benthos, no.height = no_height, no.width = no_width, patch_density = patch_density)
  spemap <- config_lst$spemap00
  eggmap <- config_lst$eggmap00
  smallmat <- config_lst$smallmat
  sum(smallmat)
  spe_config <- config_lst$plot1
  spe_config
  egg_config <- config_lst$plot2
  plot_fun1(spemap, benthos)
  plot_fun2(eggmap, benthos)
  coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
  levelplot(t(apply(smallmat[], 2, rev)),
            col.regions = coul, xlab = "Transverse",
            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
  )

  ##adults to bundles###
  spemap0 <- spemap * fecund
  eggmap0 <- eggmap * fecund
  p0 <- plot_fun1(spemap0, benthos)

  ## Advection
  h_s <- rev(0:(z_dim - 3) + 0.5) / down_scale
  drag_coef <- data.frame(depth_m = 1:10, dc = c(0.033, 0.019, 0.015, 0.012, 0.011, 0.010, 0.009, 0.009, 0.008, 0.008))
  drag_coef$shear <- sqrt(drag_coef$dc * (flow_rate / down_scale)^2)
  shear_surf <- drag_coef$shear[depth_m]
  z0 <- 0.05
  flow_rate_z <- ((shear_surf / van_K_c) * log(h_s / z0) * down_scale)
  advect_arrays <- lapply(flow_rate_z, function(fr) {
    as.integer(seq(1 / fr, no_timestep, 1 / fr))
  })

  ## Sperm sinking
  sperm_sink_cm_h <- 60
  sperm_sink_m_h <- sperm_sink_cm_h / 10
  sperm_sink_m_s <- sperm_sink_m_h / 3600

  ## Diffusion coefs
  D_L <- long_const * depth_m * shear_surf
  D_T <- trans_const * depth_m * shear_surf
  D_Z <- vert_const * depth_m * shear_surf

  ## Ascent array bund
  ascent_array <- as.integer(seq(1 / bund_asc, no_timestep, 1 / bund_asc)[1:(z_dim - 3)])
  ascent_array_egg <- as.integer(seq(1 / bund_asc_egg, no_timestep, 1 / bund_asc_egg)[1:(z_dim - 3)])[!is.na(as.integer(seq(1 / bund_asc_egg, no_timestep, 1 / bund_asc_egg)[1:(z_dim - 3)]))]

  # model prep-----------------------------------------------------------------

  ## blank maps###
  spemap3_3 <- spemap3_2 <- spemap3_1 <- spemap3_0 <- spemap3 <- spemap1_3 <- spemap2 <- spemap2_1 <- spemap1_2 <- blank_spemap
  fertmap <- polymap <- embmap <- eggmap3_4 <- eggmap3_3 <- eggmap3_2 <- eggmap3_1 <- eggmap3_0 <- eggmap3 <- eggmap2 <- eggmap1_3 <- eggmap2_1 <- eggmap1_2 <- blank_eggmap

  ## Create 10 empty lists for sperm and egg temp store
  labels_temp <- c("1", "1_0", "1_1", "1_2", "1_3", "2", "2_0", "2_1", "2_2", "3", "3_0", "3_1", "3_2", "3_3", "3_4")
  for (i in 1:length(labels_temp)) {
    list_name <- paste0("spemap_temp", labels_temp[i])
    list_name1 <- paste0("eggmap_temp", labels_temp[i])
    assign(list_name, list())
    assign(list_name1, list())
  }
  embmap_store <- list()
  phi_mono_store <- list()
  polymap_store <- list()
  fert_store <- list()
  spe_count_mat_store <- list()
  egg_count_mat_store <- list()
  fertcounter <- rep(0, (no_timestep) / 1)
  polycounter <- rep(0, (no_timestep) / 1)
  spemap1 <- spemap0
  eggmap1 <- eggmap0 * egg_viab
  tot_spe_bund <- sum(spemap1, na.rm = T)
  tot_egg_bund <- sum(eggmap1, na.rm = T)
  tot_sperm <- tot_spe_bund * speperb
  tot_eggs <- tot_egg_bund * eggperb
  bunds_sep_spe <- tot_spe_bund * bund_probs
  bunds_sep_egg <- tot_egg_bund * bund_probs
  bund_store <- numeric()

  # Start model -------------------------------------------------------------
  strt <- Sys.time()
  store_output <- T
  while (time <= no_timestep)
  {
    spemap1_1 <- spemap1
    eggmap1_1 <- eggmap1

    ## Advection loop
    spemap1_2a <- blank_spemap
    for (i in 2:(long_dim - 1)) {
      for (j in 2:(trans_dim - 1)) {
        for (k in 2:(z_dim - 1)) {
          if (!is.na(spemap1_1[i, j, k])) {
            if (time %in% advect_arrays[[k - 1]]) {
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
    spemap1_2b <- blank_spemap
    for (i in 1:(long_dim))
    {
      for (j in 1:(trans_dim))
      {
        for (k in 1:(z_dim))
        {
          if (is.na(spemap1_2a[i, j, k]) == FALSE) {
            dispdata <- dispersal5_diff(x = spemap1_2a, i = i, j = j, k = k, dx = dx, dy = dy, dz = dz, dt = dt, D_L = D_L, D_T = D_T, D_Z = D_Z, shear_surf = shear_surf, long_dim = long_dim, trans_dim = trans_dim, z_dim = z_dim, buoyant = TRUE)
            spemap1_2b[i, j, k] <- spemap1_2a[i, j, k] + dispdata
          }
        }
      }
    }

    ## no-flux boundary condition for i = 1
    spemap1_2b[2, , ] <- spemap1_2b[2, , ] + spemap1_2b[1, , ]
    spemap1_2b[1, , ] <- 0
    spemap1_2b[dim(spemap1_2b)[1] - 1, , ] <- spemap1_2b[dim(spemap1_2b)[1] - 1, , ] + spemap1_2b[dim(spemap1_2b)[1], , ]
    spemap1_2b[dim(spemap1_2b)[1], , ] <- 0
    spemap1_2b[, 2, ] <- spemap1_2b[, 2, ] + spemap1_2b[, 1, ]
    spemap1_2b[, 1, ] <- 0
    spemap1_2b[, dim(spemap1_2b)[2] - 1, ] <- spemap1_2b[, dim(spemap1_2b)[2] - 1, ] + spemap1_2b[, dim(spemap1_2b)[2], ]
    spemap1_2b[, dim(spemap1_2b)[2], ] <- 0
    spemap1_2b[, , 2] <- spemap1_2b[, , 2] + spemap1_2b[, , 1]
    spemap1_2b[, , 1] <- 0
    spemap1_2b[, , dim(spemap1_2b)[3] - 1] <- spemap1_2b[, , dim(spemap1_2b)[3] - 1] + spemap1_2b[, , dim(spemap1_2b)[3]]
    spemap1_2b[, , dim(spemap1_2b)[3]] <- 0

    ## Conservation of mass correction
    spemap1_2b <- spemap1_2b * (sum(spemap1_1, na.rm = TRUE) / sum(spemap1_2b, na.rm = TRUE))
    if (store_output) {
      spemap_temp1_1[[time]] <- spemap1_1
      spemap_temp1_2[[time]] <- spemap1_2b
    }
    eggmap1_2a <- blank_eggmap
    for (i in 2:(long_dim - 1)) {
      for (j in 2:(trans_dim - 1)) {
        for (k in 2:(z_dim - 1)) {
          if (!is.na(eggmap1_1[i, j, k])) {
            if (time %in% advect_arrays[[k - 1]]) {
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
    eggmap1_2b <- blank_spemap
    for (i in 1:(long_dim))
    {
      for (j in 1:(trans_dim))
      {
        for (k in 1:(z_dim))
        {
          if (is.na(eggmap1_2a[i, j, k]) == FALSE) {
            dispdata <- dispersal5_diff(x = eggmap1_2a, i = i, j = j, k = k, dx = dx, dy = dy, dz = dz, dt = dt, D_L = D_L, D_T = D_T, D_Z = D_Z, shear_surf = shear_surf, long_dim = long_dim, trans_dim = trans_dim, z_dim = z_dim, buoyant = TRUE)
            eggmap1_2b[i, j, k] <- eggmap1_2a[i, j, k] + dispdata
          }
        }
      }
    }

    ## no-flux boundary condition for i = 1
    eggmap1_2b[2, , ] <- eggmap1_2b[2, , ] + eggmap1_2b[1, , ]
    eggmap1_2b[1, , ] <- 0
    eggmap1_2b[dim(eggmap1_2b)[1] - 1, , ] <- eggmap1_2b[dim(eggmap1_2b)[1] - 1, , ] + eggmap1_2b[dim(eggmap1_2b)[1], , ]
    eggmap1_2b[dim(eggmap1_2b)[1], , ] <- 0
    eggmap1_2b[, 2, ] <- eggmap1_2b[, 2, ] + eggmap1_2b[, 1, ]
    eggmap1_2b[, 1, ] <- 0
    eggmap1_2b[, dim(eggmap1_2b)[2] - 1, ] <- eggmap1_2b[, dim(eggmap1_2b)[2] - 1, ] + eggmap1_2b[, dim(eggmap1_2b)[2], ]
    eggmap1_2b[, dim(eggmap1_2b)[2], ] <- 0
    eggmap1_2b[, , 2] <- eggmap1_2b[, , 2] + eggmap1_2b[, , 1]
    eggmap1_2b[, , 1] <- 0
    eggmap1_2b[, , dim(eggmap1_2b)[3] - 1] <- eggmap1_2b[, , dim(eggmap1_2b)[3] - 1] + eggmap1_2b[, , dim(eggmap1_2b)[3]]
    eggmap1_2b[, , dim(eggmap1_2b)[3]] <- 0

    ## Conservation of mass correction
    eggmap1_2b <- eggmap1_2b * (sum(spemap1_1, na.rm = TRUE) / sum(eggmap1_2b, na.rm = TRUE))
    if (store_output) {
      eggmap_temp1_1[[time]] <- eggmap1_1
      eggmap_temp1_2[[time]] <- eggmap1_2b
    }

    ## 1_3 ascent bund####

    ## ascent water column sperm bundles at various time points (1_3)
    spemap1_3 <- blank_spemap
    if (time %in% ascent_array) {
      for (i in 1:long_dim) {
        for (j in 1:trans_dim) {
          for (k in 2:(z_dim - 1)) {
            spemap1_3[i, j, (k - 1)] <- spemap1_2b[i, j, k]
          }
        }
      }
      spemap1_3[, , 2] <- spemap1_3[, , 1] + spemap1_3[, , 2]
      spemap1_3[, , 1] <- 0
      spemap1_3[, , z_dim] <- 0
    } else {
      spemap1_3[1:long_dim, 1:trans_dim, 2:(z_dim - 1)] <- spemap1_2b[1:long_dim, 1:trans_dim, 2:(z_dim - 1)]
    }
    if (store_output) {
      spemap_temp1_3[[time]] <- spemap1_3
    }

    ## ascent water column egg bundles at various time points (1_3)
    eggmap1_3 <- blank_eggmap
    if (time %in% ascent_array) {
      for (i in 1:long_dim) {
        for (j in 1:trans_dim) {
          for (k in 2:(z_dim - 1)) {
            eggmap1_3[i, j, (k - 1)] <- eggmap1_2b[i, j, k]
          }
        }
      }
      eggmap1_3[, , 2] <- eggmap1_3[, , 1] + eggmap1_3[, , 2]
      eggmap1_3[, , 1] <- 0
      eggmap1_3[, , z_dim] <- 0
    } else {
      eggmap1_3[1:long_dim, 1:trans_dim, 2:(z_dim - 1)] <- eggmap1_2b[1:long_dim, 1:trans_dim, 2:(z_dim - 1)]
    }
    if (store_output) {
      eggmap_temp1_3[[time]] <- eggmap1_3
    }
    spemap2_1 <- spemap1_3
    eggmap2_1 <- eggmap1_3
    if (store_output) {
      spemap_temp2_1[[time]] <- spemap2_1
    }
    if (store_output) {
      eggmap_temp2_1[[time]] <- eggmap2_1
    }

    ##dissociated separated bundles####

    ## Prepare maps for bundle dissociation. Separate popped bundles to map 3s
    spemap3 <- blank_spemap
    tot_bund_spemap2_1 <- sum(spemap2_1, na.rm = T)
    bundle_frac_spe <- bunds_sep_spe[time] / tot_spe_bund
    remove_spe <- spemap2_1 * bundle_frac_spe
    bund_store[time] <- sum(remove_spe, na.rm = T) * speperb
    spemap2_1 <- spemap2_1 - remove_spe
    spemap3 <- spemap3 + remove_spe
    spemap1 <- spemap2_1
    spemap2_2 <- spemap2_1
    if (store_output) {
      spemap_temp2_2[[time]] <- spemap2_2
      spemap_temp3[[time]] <- spemap3
    }

    ## Prepare maps for bundle dissociation. Move popped bundles to map 3s
    eggmap3 <- blank_spemap
    tot_bund_eggmap2_1 <- sum(eggmap2_1, na.rm = T)
    bundle_frac_egg <- bunds_sep_egg[time] / tot_bund_eggmap2_1
    remove_egg <- eggmap2_1 * bundle_frac_egg
    eggmap2_1 <- eggmap2_1 - remove_egg
    eggmap3 <- eggmap3 + remove_egg
    eggmap1 <- eggmap2_1
    eggmap2_2 <- eggmap2_1
    if (store_output) {
      eggmap_temp2_2[[time]] <- eggmap2_2
      eggmap_temp3[[time]] <- eggmap3
    }

    ##popped + existing gametes####

    ## dissociation sperm bundles
    spemap3_1 <- spemap3_0
    spe_count_mat <- spemap3 * speperb
    spemap3_1 <- spemap3_0 + spe_count_mat
    if (store_output) {
      spemap_temp3_0[[time]] <- spemap3_0
      spe_count_mat_store[[time]] <- spe_count_mat
      spemap_temp3_1[[time]] <- spemap3_1
    }
    eggmap3_1 <- eggmap3_0
    egg_count_mat <- eggmap3 * eggperb
    eggmap3_1 <- eggmap3_0 + egg_count_mat
    if (store_output) {
      eggmap_temp3_0[[time]] <- eggmap3_0
      eggmap_temp3_1[[time]] <- eggmap3_1
    }

    ## advect
    spemap3_2a <- blank_spemap
    for (i in 2:(long_dim - 1)) {
      for (j in 2:(trans_dim - 1)) {
        for (k in 2:(z_dim - 1)) {
          if (!is.na(spemap3_1[i, j, k])) {
            if (time %in% advect_arrays[[k - 1]]) {
              spemap3_2a[(i - 1), j, k] <- spemap3_1[i, j, k]
            } else {
              spemap3_2a[i, j, k] <- spemap3_1[i, j, k]
            }
          }
        }
      }
    }
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
        for (k in 1:(z_dim))
        {
          if (is.na(spemap3_2a[i, j, k]) == FALSE) {
            dispdata <- dispersal5_diff(x = spemap3_2a, i = i, j = j, k = k, dx = dx, dy = dy, dz = dz, dt = dt, D_L = D_L, D_T = D_T, D_Z = D_Z, shear_surf = shear_surf, long_dim = long_dim, trans_dim = trans_dim, z_dim = z_dim, buoyant = FALSE)
            spemap3_2b[i, j, k] <- spemap3_2a[i, j, k] + dispdata
          }
        }
      }
    }

    ## no-flux boundary condition for i = 1
    spemap3_2b[2, , ] <- spemap3_2b[2, , ] + spemap3_2b[1, , ]
    spemap3_2b[1, , ] <- 0
    spemap3_2b[dim(spemap3_2b)[1] - 1, , ] <- spemap3_2b[dim(spemap3_2b)[1] - 1, , ] + spemap3_2b[dim(spemap3_2b)[1], , ]
    spemap3_2b[dim(spemap3_2b)[1], , ] <- 0
    spemap3_2b[, 2, ] <- spemap3_2b[, 2, ] + spemap3_2b[, 1, ]
    spemap3_2b[, 1, ] <- 0
    spemap3_2b[, dim(spemap3_2b)[2] - 1, ] <- spemap3_2b[, dim(spemap3_2b)[2] - 1, ] + spemap3_2b[, dim(spemap3_2b)[2], ]
    spemap3_2b[, dim(spemap3_2b)[2], ] <- 0
    spemap3_2b[, , 2] <- spemap3_2b[, , 2] + spemap3_2b[, , 1]
    spemap3_2b[, , 1] <- 0
    spemap3_2b[, , dim(spemap3_2b)[3] - 1] <- spemap3_2b[, , dim(spemap3_2b)[3] - 1] + spemap3_2b[, , dim(spemap3_2b)[3]]
    spemap3_2b[, , dim(spemap3_2b)[3]] <- 0
    if (store_output) {
      spemap_temp3_2[[time]] <- spemap3_2b
    }
    spemap3_3 <- spemap3_2b

    ####
    moving_sperm <- sperm_sink_m_s * spemap3_2b[, , 2]
    spemap3_3[, , 3] <- spemap3_2b[, , 3] + moving_sperm
    spemap3_3[, , 2] <- spemap3_2b[, , 2] - moving_sperm

    ####

    ## advection  eggs
    eggmap3_2a <- blank_eggmap
    for (i in 2:(long_dim - 1)) {
      for (j in 2:(trans_dim - 1)) {
        for (k in 2:(z_dim - 1)) {
          if (!is.na(eggmap3_1[i, j, k])) {
            if (time %in% advect_arrays[[k - 1]]) {
              eggmap3_2a[(i - 1), j, k] <- eggmap3_1[i, j, k]
            } else {
              eggmap3_2a[i, j, k] <- eggmap3_1[i, j, k]
            }
          }
        }
      }
    }
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
        for (k in 1:(z_dim))
        {
          if (is.na(eggmap3_2a[i, j, k]) == FALSE) {
            dispdata <- dispersal5_diff(x = eggmap3_2a, i = i, j = j, k = k, dx = dx, dy = dy, dz = dz, dt = dt, D_L = D_L, D_T = D_T, D_Z = D_Z, shear_surf = shear_surf, long_dim = long_dim, trans_dim = trans_dim, z_dim = z_dim, buoyant = TRUE)
            eggmap3_2b[i, j, k] <- eggmap3_2a[i, j, k] + dispdata
          }
        }
      }
    }

    ## no-flux boundary condition for i = 1
    eggmap3_2b[2, , ] <- eggmap3_2b[2, , ] + eggmap3_2b[1, , ]
    eggmap3_2b[1, , ] <- 0
    eggmap3_2b[dim(eggmap3_2b)[1] - 1, , ] <- eggmap3_2b[dim(eggmap3_2b)[1] - 1, , ] + eggmap3_2b[dim(eggmap3_2b)[1], , ]
    eggmap3_2b[dim(eggmap3_2b)[1], , ] <- 0
    eggmap3_2b[, 2, ] <- eggmap3_2b[, 2, ] + eggmap3_2b[, 1, ]
    eggmap3_2b[, 1, ] <- 0
    eggmap3_2b[, dim(eggmap3_2b)[2] - 1, ] <- eggmap3_2b[, dim(eggmap3_2b)[2] - 1, ] + eggmap3_2b[, dim(eggmap3_2b)[2], ]
    eggmap3_2b[, dim(eggmap3_2b)[2], ] <- 0
    eggmap3_2b[, , 2] <- eggmap3_2b[, , 2] + eggmap3_2b[, , 1]
    eggmap3_2b[, , 1] <- 0
    eggmap3_2b[, , dim(eggmap3_2b)[3] - 1] <- eggmap3_2b[, , dim(eggmap3_2b)[3] - 1] + eggmap3_2b[, , dim(eggmap3_2b)[3]]
    eggmap3_2b[, , dim(eggmap3_2b)[3]] <- 0
    if (store_output) {
      eggmap_temp3_2[[time]] <- eggmap3_2b
    }

    ## ascent eggs only####
    eggmap3_3 <- blank_eggmap
    if (time %in% ascent_array_egg) {
      for (i in 1:long_dim) {
        for (j in 1:trans_dim) {
          for (k in 2:(z_dim - 1)) {
            eggmap3_3[i, j, (k - 1)] <- eggmap3_2b[i, j, k]
          }
        }
      }
      eggmap3_3[, , 2] <- eggmap3_3[, , 1] + eggmap3_3[, , 2]
      eggmap3_3[, , 1] <- 0
      eggmap3_3[, , z_dim] <- 0
    } else {
      eggmap3_3[1:long_dim, 1:trans_dim, 2:(z_dim - 1)] <- eggmap3_2b[1:long_dim, 1:trans_dim, 2:(z_dim - 1)]
    }
    if (store_output) {
      spemap_temp3_3[[time]] <- spemap3_3
    }
    if (store_output) {
      eggmap_temp3_3[[time]] <- eggmap3_3
    }

    ## 3_4 fert####

    ## fertilisation section
    eggmap3_4 <- eggmap3_3
    if (time >= polarbody) {
      for (i in 2:(long_dim - 1))
      {
        for (j in 2:(trans_dim - 1))
        {
          for (k in 2:(z_dim - 1))
          {
            if ((spemap3_3[i, j, k] > 0) == TRUE && (eggmap3_3[i, j, k] > 0) == TRUE)
            {
              S0_m3 <- spemap3_3[i, j, k]
              S0 <- S0_m3 / cell_vol_uL
              Bo <- pi * (rnorm(1, egg_rad, egg_rad_sd)^2) * rnorm(1, spe_speed, spe_speed_sd)
              time_k <- dx
              xx <- fe * S0 / E0 * (1 - exp(-Bo * E0 * time_k))
              bb <- fe * S0 / E0 * (1 - exp(-Bo * E0 * tb))
              phi_poly <- (1 - exp(-xx) - xx * exp(-xx)) * (1 - exp(-bb))
              phi_mono <- 1 - exp(-xx) - phi_poly
              blank_phi_mono[i, j, k] <- phi_mono
              phi_mono_store[[time]] <- blank_phi_mono
              fert_no <- eggmap3_3[i, j, k] * phi_mono
              fertmap[i, j, k] <- fert_no
              fert_store[[time]] <- fertmap
              embmap[i, j, k] <- embmap[i, j, k] + fert_no
              embmap_store[[time]] <- embmap
              poly <- eggmap3_3[i, j, k] * phi_poly
              polymap[i, j, k] <- polymap[i, j, k] + poly
              polymap_store[[time]] <- polymap
              eggmap3_4[i, j, k] <- eggmap3_3[i, j, k] - fert_no - poly
              eggmap_temp3_4[[time]] <- eggmap3_4
            }
          }
        }
      }
    }
    fertcounter[time] <- fertcounter[time] + (sum(embmap[, , ], na.rm = T) / tot_eggs)
    spemap3_0 <- spemap3_3
    eggmap3_0 <- eggmap3_4
    print(c("time_s", time))
    print(c("fert_success", fertcounter[time]))
    print(c("poly", sum(polymap[, , ], na.rm = T)))
    print(c("embryos", sum(embmap[, , ], na.rm = T)))
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
    time <- time + 1
    flush.console()
  }
  end_time_opt <- Sys.time()
  time_diff <- end_time_opt - strt
  time_diff

  # 7 Model end - run above ---------------------------------------------------------------

  #####PASTE TO HERE########################################
  print(vec_x)
  print(c('int_col_spac', int_col_spac))
  print(c('colony_diam', colony_diam))
  print(c('patch_type', patch_type))
  print(c('trans_dim', trans_dim))
  print(c('folder', dir_name))
  print(c('bundlebreak_sd', bundlebreak_sd))
  print(c('spe_speed', spe_speed))
  save.image(file = file.path(path = new_folder_path, paste0("last_run_", 'int_col_spac', int_col_spac, ".Rdata")))
  p0 = plot_fun1(spemap3_3, surface)
  p1 = plot_fun2(eggmap3_4, surface)
  coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
  p2 = levelplot(log10(t(apply(spemap_temp3_2[[no_timestep]][,,2], 2, rev))/10^6), col.regions = coul, at = c(0, 1, 2, 3, 3.5, 4, 5, 6), labels = T, cuts  = 5,
                 contour = T, region = T,xlab = "Transverse", ylab = "Longitudinal", main = "Final sperm concentration (log10(sperm/mL))")
  lay <- rbind(c(0,0,1,1,2,2),
               c(0,0,1,1,2,2),
               c(0,0,1,1,2,2),
               c(0,0,1,1,2,2),
               c(0,0,1,1,2,2),
               c(0,0,1,1,2,2)
  )
  gs = list(p0, p1, p2)
  plots1 = grid.arrange(grobs = gs, layout_matrix = lay)

  ####fert counter####
  fert_df = data.frame(time = 1:(no_timestep), fert = fertcounter)
  tail(fert_df)
  max_fert = max(fert_df$fert, na.rm = T)
  filename <- file.path(new_folder_path, paste0("fertcounter_", int_col_spac, ".RData"))
  source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2")
  p4 = ggplot()+geom_point(fert_df, mapping = aes(x = time, y = fert),position = position_jitter(width = .00, height = .00),
                           alpha = 0.50, size = 3) + theme_sleek2()
  p4 <- p4 + annotate("text", x = max(fert_df$time, na.rm = T)*0.3, y = max(fert_df$fert, na.rm = T)*0.5,
                      label = paste("Max fert:", round(max_fert*100, 2), '%' ),
                      hjust = 1, vjust = -1)
  p4 <- p4 + annotate("text", x = max(fert_df$time, na.rm = T)*0.3, y = max(fert_df$fert, na.rm = T)*0.4,
                      label = paste("Total embryos:", round(sum(embmap_store[[no_timestep-1]][,,2], na.rm = T),0) ),
                      hjust = 1, vjust = -1)
  p4 = p4 + labs(x=expression(Time~(s)),
                 y=expression(Cumulative~fert.~success~(prop.)))
  ggsave(file = paste0("clouds_", 'int_col_spac', int_col_spac,".pdf"), plots1, width = 10, height = 8, path = new_folder_path)
  ggsave(p4, file = paste0("fert_counter_",'int_col_spac', int_col_spac,".pdf"), width = 10, height = 8, path = new_folder_path)
}
