config <- function(spemap, eggmap, long_dim, trans_dim, z_dim, int_col_spac, patch_density, fecun_zone, polyp_den, colony_diam, benthos,
                   config1 = "patch_bot", no_height = 5, no_width = 4, target_position = 1, matrix_input = NULL, den_cell = 1) {
  col_area <- pi * (colony_diam / 2)^2
  fecund <- col_area * polyp_den * fecun_zone
  if (config1 == "patch_bot") {
    no_height <- no_height
    no_width <- no_width
    height_dim <- no_height * int_col_spac
    width_dim <- no_width * int_col_spac
    mat <- matrix(0, nrow = height_dim, ncol = width_dim)
    for (i in 1:dim(mat)[1]) {
      for (j in 1:dim(mat)[2]) {
        if ((row(mat)[i, j] %% int_col_spac) == 0 & (col(mat)[i, j] %% int_col_spac) == 0) {
          mat[i, j] <- mat[i, j] + 1
        }
      }
    }
    mat
    row_col <- row(mat)[which(!mat == 0)]
    col_col <- col(mat)[which(!mat == 0)]
    mid_row <- floor(long_dim / 2) + 1
    mid_col <- floor(trans_dim / 2) + 1
    small_mat <- mat[min(row_col):max(row_col), min(col_col):max(col_col)]
    small_mat
    small_mat[small_mat == 1] <- small_mat[small_mat == 1] * den_cell
    small_mat_egg <- small_mat
    small_mat_egg[] <- 0
    center <- c(nrow(small_mat) %/% 2 + 1, ncol(small_mat) %/% 2 + 1)
    indices <- which(small_mat > 0, arr.ind = TRUE)
    distances <- sqrt((indices[, 1] - center[1])^2 + (indices[, 2] - center[2])^2)
    closest <- which.min(distances)
    closest_position <- indices[closest, ]
    small_mat[closest_position[1], closest_position[2]] <- small_mat[closest_position[1], closest_position[2]] - 1
    small_mat_egg[closest_position[1], closest_position[2]] <- 1
    col_number <- sum(small_mat != 0)
    col_dens <- (height_dim * width_dim) / col_number
    row_index <- floor(nrow(spemap) * 9 / 10) - floor(nrow(small_mat) / 2)
    col_index <- mid_col - floor(ncol(small_mat) / 2)
    spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    spe_config <- levelplot(t(apply(spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                            col.regions = coul, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat_egg
    egg_config <- levelplot(t(apply(eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                            col.regions = coul2, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    spemap0 = spemap * fecund
    eggmap0 = eggmap * fecund
    lst1 <- list("spemap00" = spemap0, "eggmap00" = eggmap0, smallmat = small_mat, plot1 = spe_config, plot2 = egg_config)
    return(lst1)
  }
  if (config1 == "patch_bot_fine") {
    if(int_col_spac < 1){
      patch_dim = round(int_col_spac * (sqrt(patch_density)-1),0)
      no_height <- patch_dim
      no_width <- patch_dim
      converter = patch_density/patch_dim^2
      height_dim <- no_height
      width_dim <- no_width
      mat <- matrix(0, nrow = height_dim, ncol = width_dim)
      mid_row <- floor(long_dim / 2) + 1
      mid_col <- floor(trans_dim / 2) + 1
      small_mat = mat + converter
      small_mat_egg <- small_mat
      small_mat_egg[] <- 0
      center <- c(nrow(small_mat) %/% 2 + 1, ncol(small_mat) %/% 2 + 1)
      indices <- which(small_mat > 0, arr.ind = TRUE)
      distances <- sqrt((indices[, 1] - center[1])^2 + (indices[, 2] - center[2])^2)
      closest <- which.min(distances)
      closest_position <- indices[closest, ]
      small_mat[closest_position[1], closest_position[2]] <- small_mat[closest_position[1], closest_position[2]] - 1
      small_mat_egg[closest_position[1], closest_position[2]] <- 1
      col_number <- sum(small_mat != 0)
      col_dens <- (height_dim * width_dim) / col_number
      row_index <- floor(nrow(spemap) * 9 / 10) - floor(nrow(small_mat) / 2)
      col_index <- mid_col - floor(ncol(small_mat) / 2)
      spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat
      coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
      coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
      spe_config <- levelplot(t(apply(spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                              col.regions = coul, xlab = "Transverse",
                              ylab = "Longitudinal", main = "Conc. (cells/m^3)"
      )
      eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat_egg
      egg_config <- levelplot(t(apply(eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                              col.regions = coul2, xlab = "Transverse",
                              ylab = "Longitudinal", main = "Conc. (cells/m^3)"
      )
      spemap0 = spemap * fecund
      eggmap0 = eggmap * fecund
      lst1 <- list("spemap00" = spemap0, "eggmap00" = eggmap0, smallmat = small_mat, plot1 = spe_config, plot2 = egg_config)
      return(lst1)
    }else{
      no_height <- no_height
      no_width <- no_width
      height_dim <- no_height * int_col_spac
      width_dim <- no_width * int_col_spac
      mat <- matrix(0, nrow = height_dim, ncol = width_dim)
      for (i in 1:dim(mat)[1]) {
        for (j in 1:dim(mat)[2]) {
          if ((row(mat)[i, j] %% int_col_spac) == 0 & (col(mat)[i, j] %% int_col_spac) == 0) {
            mat[i, j] <- mat[i, j] + 1
          }
        }
      }
      mat
      row_col <- row(mat)[which(!mat == 0)]
      col_col <- col(mat)[which(!mat == 0)]
      mid_row <- floor(long_dim / 2) + 1
      mid_col <- floor(trans_dim / 2) + 1
      small_mat <- mat[min(row_col):max(row_col), min(col_col):max(col_col)]
      small_mat[small_mat == 1] <- small_mat[small_mat == 1] * den_cell
      small_mat_egg <- small_mat
      small_mat_egg[] <- 0
      center <- c(nrow(small_mat) %/% 2 + 1, ncol(small_mat) %/% 2 + 1)
      indices <- which(small_mat > 0, arr.ind = TRUE)
      distances <- sqrt((indices[, 1] - center[1])^2 + (indices[, 2] - center[2])^2)
      closest <- which.min(distances)
      closest_position <- indices[closest, ]
      small_mat[closest_position[1], closest_position[2]] <- small_mat[closest_position[1], closest_position[2]] - 1
      small_mat_egg[closest_position[1], closest_position[2]] <- 1
      col_number <- sum(small_mat != 0)
      col_dens <- (height_dim * width_dim) / col_number
      row_index <- floor(nrow(spemap) * 9 / 10) - floor(nrow(small_mat) / 2)
      col_index <- mid_col - floor(ncol(small_mat) / 2)
      spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat
      coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
      coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
      spe_config <- levelplot(t(apply(spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                              col.regions = coul, xlab = "Transverse",
                              ylab = "Longitudinal", main = "Conc. (cells/m^3)"
      )
      eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat_egg
      egg_config <- levelplot(t(apply(eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                              col.regions = coul2, xlab = "Transverse",
                              ylab = "Longitudinal", main = "Conc. (cells/m^3)"
      )
      spemap0 = spemap * fecund
      eggmap0 = eggmap * fecund
      lst1 <- list("spemap00" = spemap0, "eggmap00" = eggmap0, smallmat = small_mat, plot1 = spe_config, plot2 = egg_config)
      return(lst1)
    }
  }
  if (config1 == "2023_palau") {
    load("./Rdata/2023palau_coarse_grid_centre.RData")
    centre_patch = grid_coarse
    grid_coarse = matrix_input
    target_ind = grid_coarse
    dim(centre_patch)
    dim(target_ind)
    dim(spemap)
    bot_row_index <- floor(nrow(spemap) * 10 / 10)
    top_row_index <- bot_row_index - nrow(centre_patch)  +1
    mid_col <- floor(trans_dim / 2) + 1
    first_col_index <- mid_col - floor(ncol(centre_patch) / 2)
    last_col_index <- first_col_index + ncol(centre_patch)  - 1
    dim(spemap)
    range_x = top_row_index:bot_row_index
    range_y = first_col_index:last_col_index
    dim(centre_patch)
    length(range_x)
    length(range_y)
    spemap[range_x, range_y, benthos] <- centre_patch
    eggmap[range_x, range_y, benthos] <- target_ind
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    spe_config <- levelplot(t(apply(spemap[, , benthos], 2, rev)),
                            col.regions = coul, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    egg_config <- levelplot(t(apply(eggmap[, , benthos], 2, rev)),
                            col.regions = coul2, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    spemap0 = spemap * fecund
    eggmap0 = eggmap * fecund
    lst1 <- list("spemap00" = spemap0, "eggmap00" = eggmap0, plot1 = spe_config, plot2 = egg_config)
    return(lst1)
  }
  if (config1 == "2022ahyac") {
    no_height <- 4
    no_width <- 5
    int_col_spac = 1
    height_dim <- no_height * int_col_spac
    width_dim <- no_width * int_col_spac
    mat <- matrix(0, nrow = height_dim, ncol = width_dim)
    for (i in 1:dim(mat)[1]) {
      for (j in 1:dim(mat)[2]) {
        if ((row(mat)[i, j] %% int_col_spac) == 0 & (col(mat)[i, j] %% int_col_spac) == 0) {
          mat[i, j] <- mat[i, j] + 1
        }
      }
    }
    mat
    row_col <- row(mat)[which(!mat == 0)]
    col_col <- col(mat)[which(!mat == 0)]
    mid_row <- floor(long_dim / 2) + 1
    mid_col <- floor(trans_dim / 2) + 1
    small_mat <- mat[min(row_col):max(row_col), min(col_col):max(col_col)]
    small_mat
    small_mat_egg <- small_mat
    small_mat_egg[] <- 0
    center <- c(nrow(small_mat) %/% 2 + 1, ncol(small_mat) %/% 2 + 1)
    indices <- which(small_mat > 0, arr.ind = TRUE)
    distances <- sqrt((indices[, 1] - center[1])^2 + (indices[, 2] - center[2])^2)
    closest <- which.min(distances)
    closest_position <- indices[closest, ]
    small_mat[closest_position[1], closest_position[2]] <- 0
    small_mat_egg[closest_position[1], closest_position[2]] <- 1
    col_number <- sum(small_mat != 0)
    col_dens <- (height_dim * width_dim) / col_number
    row_index <- floor(nrow(spemap) * 9 / 10) - floor(nrow(small_mat) / 2)
    col_index <- mid_col - floor(ncol(small_mat) / 2)
    spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    spe_config <- levelplot(t(apply(spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                            col.regions = coul, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat_egg
    egg_config <- levelplot(t(apply(eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                            col.regions = coul2, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    spemap0 = spemap * fecund
    eggmap0 = eggmap * fecund
    lst1 <- list("spemap00" = spemap0, "eggmap00" = eggmap0, smallmat = small_mat, plot1 = spe_config, plot2 = egg_config)
    return(lst1)
  }
  if (config1 == "2021tenuis") {
    no_height <- no_height
    no_width <- no_width
    load("./Rdata/2021_ten_adu_coords_rect.RData")
    small_mat <- data1
    small_mat_egg <- small_mat
    small_mat_egg[] <- 0
    center <- c(nrow(small_mat) %/% 2 + 1, ncol(small_mat) %/% 2 + 1)
    indices <- which(small_mat > 0, arr.ind = TRUE)
    distances <- sqrt((indices[, 1] - center[1])^2 + (indices[, 2] - center[2])^2)
    closest <- which.min(distances)
    closest_position <- indices[closest, ]
    small_mat[closest_position[1], closest_position[2]] <- 0
    small_mat_egg[closest_position[1], closest_position[2]] <- 1
    height_dim <- nrow(small_mat)
    width_dim <- ncol(small_mat)
    col_number <- sum(small_mat != 0)
    col_dens <- (height_dim * width_dim) / col_number
    col_index <- floor((trans_dim - ncol(small_mat)) / 2) + 1
    buffer <- 2
    row_index <- (long_dim - buffer) - nrow(small_mat) + 1
    spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat
    eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat_egg
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    spe_config <- levelplot(t(apply(spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                            col.regions = coul, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    egg_config <- levelplot(t(apply(eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                            col.regions = coul2, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    spemap0 = spemap * fecund
    eggmap0 = eggmap * fecund
    lst1 <- list("spemap00" = spemap0, "eggmap00" = eggmap0, smallmat = small_mat, plot1 = spe_config, plot2 = egg_config)
    return(lst1)
  }
  if (config1 == "2021digit") {
    no_height <- no_height
    no_width <- no_width
    load("./Rdata/2021_dig_adu_coords_rect.RData")
    small_mat <- data1
    small_mat_egg <- small_mat
    small_mat_egg[] <- 0
    center <- c(nrow(small_mat) %/% 2 + 1, ncol(small_mat) %/% 2 + 1)
    indices <- which(small_mat > 0, arr.ind = TRUE)
    distances <- sqrt((indices[, 1] - center[1])^2 + (indices[, 2] - center[2])^2)
    small_mat[4, 3] <- 1
    closest <- which.min(distances)
    closest_position <- indices[closest, ]
    small_mat[closest_position[1], closest_position[2]] <- 0
    small_mat_egg[closest_position[1], closest_position[2]] <- 1

    ##random selection
    height_dim <- nrow(small_mat)
    width_dim <- ncol(small_mat)
    col_number <- sum(small_mat != 0)
    col_dens <- (height_dim * width_dim) / col_number
    col_index <- floor((trans_dim - ncol(small_mat)) / 2) + 1
    buffer <- 2
    row_index <- (long_dim - buffer) - nrow(small_mat) + 1
    spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat
    eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat_egg
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    spe_config <- levelplot(t(apply(spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                            col.regions = coul, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    egg_config <- levelplot(t(apply(eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                            col.regions = coul2, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    spemap0 = spemap * fecund
    eggmap0 = eggmap * fecund
    lst1 <- list("spemap00" = spemap0, "eggmap00" = eggmap0, smallmat = small_mat, plot1 = spe_config, plot2 = egg_config)
    return(lst1)
  }
  if (config1 == "patch_outer") {
    no_height <- no_height
    no_width <- no_width
    height_dim <- no_height * int_col_spac
    width_dim <- no_width * int_col_spac
    mat <- matrix(0:0, nrow = height_dim, ncol = width_dim)
    for (i in 1:dim(mat)[1]) {
      for (j in 1:dim(mat)[2]) {
        if ((row(mat)[i, j] %% int_col_spac) == 0 & (col(mat)[i, j] %% int_col_spac) == 0) {
          mat[i, j] <- mat[i, j] + 1
        }
      }
    }
    mat
    row_col <- row(mat)[which(!mat == 0)]
    col_col <- col(mat)[which(!mat == 0)]
    mid_row <- floor(long_dim / 2) + 1
    mid_col <- floor(trans_dim / 2) + 1
    small_mat <- mat[min(row_col):max(row_col), min(col_col):max(col_col)]
    small_mat
    small_mat_egg <- small_mat
    small_mat_egg[] <- 0
    center <- c(nrow(small_mat) %/% 2 + 1, ncol(small_mat) %/% 2 + 1)
    indices <- which(small_mat > 0, arr.ind = TRUE)
    distances <- sqrt((indices[, 1] - center[1])^2 + (indices[, 2] - center[2])^2)
    col_number <- sum(small_mat != 0)
    col_dens <- (height_dim * width_dim) / col_number
    row_index <- floor(nrow(spemap) * 3 / 4) - floor(nrow(small_mat) / 2)
    col_index <- mid_col - floor(ncol(small_mat) / 2)
    spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat
    row_col_all <- row(spemap[, , 2])[which(!spemap[, , benthos] == 0)]
    col_col_all <- col(spemap[, , 2])[which(!spemap[, , benthos] == 0)]
    spemap[row_col_all[1], col_col_all[1], benthos] <- 0
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    spe_config <- levelplot(t(apply(spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                            col.regions = coul, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    eggmap[row_col_all[1], col_col_all[1], benthos] <- 1
    egg_config <- levelplot(t(apply(eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                            col.regions = coul2, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    spemap0 = spemap * fecund
    eggmap0 = eggmap * fecund
    lst1 <- list("spemap00" = spemap0, "eggmap00" = eggmap0, smallmat = small_mat, plot1 = spe_config, plot2 = egg_config)
    return(lst1)
  }
  if (config1 == "patch_center") {
    no_height <- no_height
    no_width <- no_width
    height_dim <- no_height * int_col_spac
    width_dim <- no_width * int_col_spac
    mat <- matrix(0:0, nrow = height_dim, ncol = width_dim)
    for (i in 1:dim(mat)[1]) {
      for (j in 1:dim(mat)[2]) {
        if ((row(mat)[i, j] %% int_col_spac) == 0 & (col(mat)[i, j] %% int_col_spac) == 0) {
          mat[i, j] <- mat[i, j] + 1
        }
      }
    }
    mat
    row_col <- row(mat)[which(!mat == 0)]
    col_col <- col(mat)[which(!mat == 0)]
    mid_row <- floor(long_dim / 2) + 1
    mid_col <- floor(trans_dim / 2) + 1
    small_mat <- mat[min(row_col):max(row_col), min(col_col):max(col_col)]
    small_mat
    small_mat_egg <- small_mat
    small_mat_egg[] <- 0
    center <- c(nrow(small_mat) %/% 2 + 1, ncol(small_mat) %/% 2 + 1)
    indices <- which(small_mat > 0, arr.ind = TRUE)
    distances <- sqrt((indices[, 1] - center[1])^2 + (indices[, 2] - center[2])^2)
    closest <- which.min(distances)
    closest_position <- indices[closest, ]
    small_mat[closest_position[1], closest_position[2]] <- 0
    small_mat_egg[closest_position[1], closest_position[2]] <- 1
    closest <- which.min(distances)
    closest_position <- indices[closest, ]
    small_mat[closest_position[1], closest_position[2]] <- 0
    small_mat_egg[closest_position[1], closest_position[2]] <- 1
    col_number <- sum(small_mat != 0)
    col_dens <- (height_dim * width_dim) / col_number
    row_index <- mid_row - floor(nrow(small_mat) / 2)
    col_index <- mid_col - floor(ncol(small_mat) / 2)
    spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    spe_config <- levelplot(t(apply(spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                            col.regions = coul, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat_egg
    egg_config <- levelplot(t(apply(eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                            col.regions = coul2, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    spemap0 = spemap * fecund
    eggmap0 = eggmap * fecund
    lst1 <- list("spemap00" = spemap0, "eggmap00" = eggmap0, smallmat = small_mat, plot1 = spe_config, plot2 = egg_config)
    return(lst1)
  }
  if (config1 == "petal") {
    mid.long <- long_dim - int_col_spac - 2
    mid.trans <- floor(trans_dim / 2 - 0.5)
    eggmap <- blank.spemap
    eggmap[mid.long, mid.trans, benthos] <- 1
    eggmap[, , 6]
    spemap[mid.long + int_col_spac, mid.trans, benthos] <- 1
    spemap[mid.long - int_col_spac, mid.trans, benthos] <- 1
    spemap[ceiling(mid.long + 0.5 * int_col_spac), ceiling(mid.trans + 17 / 20 * int_col_spac), benthos] <- 1
    spemap[ceiling(mid.long + 0.5 * int_col_spac), floor(mid.trans - 17 / 20 * int_col_spac), benthos] <- 1
    spemap[floor(mid.long - 0.5 * int_col_spac), ceiling(mid.trans + 17 / 20 * int_col_spac), benthos] <- 1
    spemap[floor(mid.long - 0.5 * int_col_spac), floor(mid.trans - 17 / 20 * int_col_spac), benthos] <- 1
    spemap0 = spemap * fecund
    eggmap0 = eggmap * fecund
    lst1 <- list("spemap00" = spemap0, "eggmap00" = eggmap0)
    return(lst1)
  }
  if (config1 == "poisson") {
    no_height <- no_height
    no_width <- no_width
    small_mat <- matrix(0, nrow = no_height, ncol = no_width)
    placed_colonies <- 0
    while (placed_colonies < patch_density) {
      x_pos <- sample(1:(no_height), 1)
      y_pos <- sample(1:(no_width), 1)
      small_mat[x_pos, y_pos] <- small_mat[x_pos, y_pos] + 1
      placed_colonies <- placed_colonies + 1
    }
    sum(small_mat)
    mid_row <- floor(long_dim / 2) + 1
    mid_col <- floor(trans_dim / 2) + 1
    row_index <- floor(nrow(spemap) * 9 / 10) - floor(nrow(small_mat) / 2)
    col_index <- mid_col - floor(ncol(small_mat) / 2)
    spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat
    ones_positions <- which(small_mat == 1, arr.ind = TRUE)
    if (nrow(ones_positions) > 0) {
      selected_pos <- ones_positions[sample(nrow(ones_positions), 1), ]
      spemap[row_index + selected_pos[1] - 1, col_index + selected_pos[2] - 1, z_dim - 1] <- 0
      eggmap[row_index + selected_pos[1] - 1, col_index + selected_pos[2] - 1, z_dim - 1] <- 1
    }
    sum(small_mat)
    sum(spemap[, , z_dim - 1])
    sum(eggmap[, , z_dim - 1])
    spemap0 = spemap * fecund
    eggmap0 = eggmap * fecund
    lst1 <- list("spemap00" = spemap0, "eggmap00" = eggmap0, smallmat = small_mat)
    return(lst1)
  }
  if (config1 == "uniform") {
    no_height <- no_height
    no_width <- no_width
    small_mat <- matrix(0, nrow = no_height, ncol = no_width)
    total_colonies <- min(patch_density, no_height * no_width)
    grid_rows <- floor(sqrt(total_colonies))
    grid_cols <- ceiling(total_colonies / grid_rows)
    step_x <- floor(no_height / grid_rows)
    step_y <- floor(no_width / grid_cols)
    extra_rows <- no_height - (step_x * grid_rows)
    extra_cols <- no_width - (step_y * grid_cols)
    placed_colonies <- 0
    for (i in 1:grid_rows) {
      for (j in 1:grid_cols) {
        if (placed_colonies < total_colonies) {
          x_pos <- ((i - 1) * step_x + 1) + floor(step_x / 2)
          y_pos <- ((j - 1) * step_y + 1) + floor(step_y / 2)
          if (i == grid_rows && extra_rows > 0) {
            x_pos <- no_height - floor(extra_rows / 2)
          }
          if (j == grid_cols && extra_cols > 0) {
            y_pos <- no_width - floor(extra_cols / 2)
          }
          small_mat[x_pos, y_pos] <- 1
          placed_colonies <- placed_colonies + 1
        }
      }
    }
    sum(small_mat)
    sum(spemap[, , z_dim - 1])
    sum(eggmap[, , z_dim - 1])
    spemap0 = spemap * fecund
    eggmap0 = eggmap * fecund
    lst1 <- list("spemap00" = spemap0, "eggmap00" = eggmap0, smallmat = small_mat)
    return(lst1)
  }
  if (config1 == "0.5") {
    no_height <- 5
    mid_col <- floor(trans_dim / 2) + 1
    spemap[long_dim - 3, mid_col - 1, benthos] <- 1
    spemap[long_dim - 3, mid_col, benthos] <- 3
    spemap[long_dim - 3, mid_col + 1, benthos] <- 1
    spemap[long_dim - 4, mid_col - 1, benthos] <- 3
    spemap[long_dim - 4, mid_col, benthos] <- 8
    spemap[long_dim - 4, mid_col + 1, benthos] <- 3
    spemap[long_dim - 5, mid_col - 1, benthos] <- 1
    spemap[long_dim - 5, mid_col, benthos] <- 3
    spemap[long_dim - 5, mid_col + 1, benthos] <- 1
    small_mat <- spemap[(long_dim - 3):(long_dim - 5), (mid_col - 1):(mid_col + 1), benthos]
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    spe_config <- levelplot(t(apply(spemap[, , benthos], 2, rev)),
                            col.regions = coul, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    eggmap[long_dim - 4, mid_col, benthos] <- 1
    egg_config <- levelplot(t(apply(eggmap[, , z_dim - 1], 2, rev)),
                            col.regions = coul2, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    spemap0 = spemap * fecund
    eggmap0 = eggmap * fecund
    lst1 <- list("spemap00" = spemap0, "eggmap00" = eggmap0, smallmat = small_mat, plot1 = spe_config, plot2 = egg_config)
    return(lst1)
  }
  if (config1 == "1.5") {
    mid_col <- floor(trans_dim / 2) + 1
    spemap[long_dim - 3, mid_col - 3, benthos] <- 1
    spemap[long_dim - 3, mid_col - 1, benthos] <- 1
    spemap[long_dim - 3, mid_col, benthos] <- 1
    spemap[long_dim - 3, mid_col + 1, benthos] <- 1
    spemap[long_dim - 3, mid_col + 3, benthos] <- 1
    spemap[long_dim - 5, mid_col - 3, benthos] <- 1
    spemap[long_dim - 5, mid_col - 1, benthos] <- 1
    spemap[long_dim - 5, mid_col, benthos] <- 1
    spemap[long_dim - 5, mid_col + 1, benthos] <- 1
    spemap[long_dim - 5, mid_col + 3, benthos] <- 1
    spemap[long_dim - 6, mid_col - 3, benthos] <- 1
    spemap[long_dim - 6, mid_col - 1, benthos] <- 1
    spemap[long_dim - 6, mid_col, benthos] <- 0
    spemap[long_dim - 6, mid_col + 1, benthos] <- 1
    spemap[long_dim - 6, mid_col + 3, benthos] <- 1
    spemap[long_dim - 7, mid_col - 3, benthos] <- 1
    spemap[long_dim - 7, mid_col - 1, benthos] <- 1
    spemap[long_dim - 7, mid_col, benthos] <- 1
    spemap[long_dim - 7, mid_col + 1, benthos] <- 1
    spemap[long_dim - 7, mid_col + 3, benthos] <- 1
    spemap[long_dim - 9, mid_col - 3, benthos] <- 1
    spemap[long_dim - 9, mid_col - 1, benthos] <- 1
    spemap[long_dim - 9, mid_col, benthos] <- 1
    spemap[long_dim - 9, mid_col + 1, benthos] <- 1
    spemap[long_dim - 9, mid_col + 3, benthos] <- 1
    small_mat <- spemap[(long_dim - 3):(long_dim - 9), (mid_col - 3):(mid_col + 3), benthos]
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    spe_config <- levelplot(t(apply(spemap[, , benthos], 2, rev)),
                            col.regions = coul, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    eggmap[long_dim - 4, mid_col, benthos] <- 1
    egg_config <- levelplot(t(apply(eggmap[, , z_dim - 1], 2, rev)),
                            col.regions = coul2, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    spemap0 = spemap * fecund
    eggmap0 = eggmap * fecund
    lst1 <- list("spemap00" = spemap0, "eggmap00" = eggmap0, smallmat = small_mat, plot1 = spe_config, plot2 = egg_config)
    return(lst1)
  }
  if (config1 == "patch_bot_min1") {
    no_height <- no_height
    no_width <- no_width
    height_dim <- no_height * int_col_spac
    width_dim <- no_width * int_col_spac
    mat <- matrix(0, nrow = height_dim, ncol = width_dim)
    for (i in 1:dim(mat)[1]) {
      for (j in 1:dim(mat)[2]) {
        if ((row(mat)[i, j] %% int_col_spac) == 0 & (col(mat)[i, j] %% int_col_spac) == 0) {
          mat[i, j] <- mat[i, j] + 1
        }
      }
    }
    mat
    row_col <- row(mat)[which(!mat == 0)]
    col_col <- col(mat)[which(!mat == 0)]
    mid_row <- floor(long_dim / 2) + 1
    mid_col <- floor(trans_dim / 2) + 1
    small_mat <- mat[min(row_col):max(row_col), min(col_col):max(col_col)]
    small_mat <- small_mat[-nrow(small_mat), ]
    small_mat_egg <- small_mat
    small_mat_egg[] <- 0
    center <- c(nrow(small_mat) %/% 2 + 1, ncol(small_mat) %/% 2 + 1)
    indices <- which(small_mat > 0, arr.ind = TRUE)
    distances <- sqrt((indices[, 1] - center[1])^2 + (indices[, 2] - center[2])^2)
    closest <- which.min(distances)
    closest_position <- indices[closest, ]
    small_mat[closest_position[1], closest_position[2]] <- small_mat[closest_position[1], closest_position[2]] - 1
    small_mat_egg[closest_position[1], closest_position[2]] <- 1
    col_number <- sum(small_mat != 0)
    col_dens <- (height_dim * width_dim) / col_number
    row_index <- floor(nrow(spemap) * 9 / 10) - floor(nrow(small_mat) / 2)
    col_index <- mid_col - floor(ncol(small_mat) / 2)
    spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    spe_config <- levelplot(t(apply(spemap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                            col.regions = coul, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1] <- small_mat_egg
    egg_config <- levelplot(t(apply(eggmap[row_index:(row_index + nrow(small_mat) - 1), col_index:(col_index + ncol(small_mat) - 1), z_dim - 1], 2, rev)),
                            col.regions = coul2, xlab = "Transverse",
                            ylab = "Longitudinal", main = "Conc. (cells/m^3)"
    )
    spemap0 = spemap * fecund
    eggmap0 = eggmap * fecund
    lst1 <- list("spemap00" = spemap0, "eggmap00" = eggmap0, smallmat = small_mat, plot1 = spe_config, plot2 = egg_config)
    return(lst1)
  }

  ## 5 x 5b 1 int reg grid. (superceeded?)

  ### poisson approach###

  ####
}
