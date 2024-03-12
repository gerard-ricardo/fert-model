config = function(spemap, eggmap, long_dim, trans_dim, z_dim, int_col_spac, patch_density, benthos, config1 = 'patch_bot', no.height =  4, no.width = 5){
  
  if(config1 == 'patch_bot') {
    no.height =  no.height
    no.width =   no.width
    height.dim = no.height * int_col_spac #+ 1  #creates enough space
    width.dim = no.width * int_col_spac #+ 1
    mat <- matrix(0:0, nrow = height.dim, ncol = width.dim)
    for (i in 1:dim(mat)[1]) {
      for (j in 1:dim(mat)[2]) {
        if((row(mat)[i,j] %% int_col_spac) == 0   & (col(mat)[i,j] %% int_col_spac) == 0) {   #%% returns the remainder of x/y i.e x/2. If true, then even number
          mat[i, j] <- mat[i, j] + 1
        }}}
    mat  #if divisible by intercolonial spacing , then return fecund
    row.col = row(mat)[which(!mat == 0)] #gets the index values for non-zeros
    col.col = col(mat)[which(!mat == 0)]
    # Compute the midpoints
    mid_row = floor(long_dim / 2) + 1
    mid_col = floor(trans_dim / 2) + 1
    #small.mat = mat[int_col_spac:(height.dim + int_col_spac-1), int_col_spac: (width.dim + int_col_spac-1) ]  #this shaves off the edges
    small.mat = mat[min(row.col): max(row.col) , min(col.col): max(col.col)]  #this shaves off the edges
    small.mat
    small.mat_egg = small.mat
    small.mat_egg[] <- 0  #make all values zero
    center <- c(nrow(small.mat) %/% 2 + 1, ncol(small.mat) %/% 2 + 1)
    indices <- which(small.mat > 0, arr.ind = TRUE)
    distances <- sqrt((indices[,1] - center[1])^2 + (indices[,2] - center[2])^2)
    closest <- which.min(distances)
    closest_position <- indices[closest, ]
    small.mat[closest_position[1],closest_position[2]] <- 0
    small.mat_egg[closest_position[1],closest_position[2]] <- 1
    col.number = sum(small.mat != 0)  #colony number
    col.dens = (height.dim * width.dim) / col.number  #colony density
    row_index = floor(nrow(spemap) * 9/10) - floor(nrow(small.mat) / 2)
    col_index = mid_col - floor(ncol(small.mat) / 2)
    spemap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1] <- small.mat
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    spe_config= levelplot(t(apply(spemap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1], 2, rev)), col.regions = coul, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    eggmap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1] <- small.mat_egg
    egg_config= levelplot(t(apply(eggmap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1], 2, rev)), col.regions = coul2, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    lst1 <- list('spemap00'= spemap, 'eggmap00'=eggmap, smallmat = small.mat, plot1 = spe_config, plot2 = egg_config) #store output in list
    return(lst1)
  }
  
  if(config1 == '2022ahyac') {
    no.height =  4
    no.width =   5
    height.dim = no.height * int_col_spac #+ 1  #creates enough space
    width.dim = no.width * int_col_spac #+ 1
    mat <- matrix(0:0, nrow = height.dim, ncol = width.dim)
    for (i in 1:dim(mat)[1]) {
      for (j in 1:dim(mat)[2]) {
        if((row(mat)[i,j] %% int_col_spac) == 0   & (col(mat)[i,j] %% int_col_spac) == 0) {   #%% returns the remainder of x/y i.e x/2. If true, then even number
          mat[i, j] <- mat[i, j] + 1
        }}}
    mat  #if divisible by intercolonial spacing , then return fecund
    row.col = row(mat)[which(!mat == 0)] #gets the index values for non-zeros
    col.col = col(mat)[which(!mat == 0)]
    # Compute the midpoints
    mid_row = floor(long_dim / 2) + 1
    mid_col = floor(trans_dim / 2) + 1
    #small.mat = mat[int_col_spac:(height.dim + int_col_spac-1), int_col_spac: (width.dim + int_col_spac-1) ]  #this shaves off the edges
    small.mat = mat[min(row.col): max(row.col) , min(col.col): max(col.col)]  #this shaves off the edges
    small.mat
    small.mat_egg = small.mat
    small.mat_egg[] <- 0  #make all values zero
    center <- c(nrow(small.mat) %/% 2 + 1, ncol(small.mat) %/% 2 + 1)
    indices <- which(small.mat > 0, arr.ind = TRUE)
    distances <- sqrt((indices[,1] - center[1])^2 + (indices[,2] - center[2])^2)
    closest <- which.min(distances)
    closest_position <- indices[closest, ]
    small.mat[closest_position[1],closest_position[2]] <- 0
    small.mat_egg[closest_position[1],closest_position[2]] <- 1
    col.number = sum(small.mat != 0)  #colony number
    col.dens = (height.dim * width.dim) / col.number  #colony density
    row_index = floor(nrow(spemap) * 9/10) - floor(nrow(small.mat) / 2)
    col_index = mid_col - floor(ncol(small.mat) / 2)
    spemap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1] <- small.mat
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    spe_config= levelplot(t(apply(spemap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1], 2, rev)), col.regions = coul, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    eggmap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1] <- small.mat_egg
    egg_config= levelplot(t(apply(eggmap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1], 2, rev)), col.regions = coul2, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    lst1 <- list('spemap00'= spemap, 'eggmap00'=eggmap, smallmat = small.mat, plot1 = spe_config, plot2 = egg_config) #store output in list
    return(lst1)
  }
  
  
  if(config1 == 'patch_outer') {
    no.height =  no.height
    no.width =   no.width
    height.dim = no.height * int_col_spac #+ 1  #creates enough space
    width.dim = no.width * int_col_spac #+ 1
    mat <- matrix(0:0, nrow = height.dim, ncol = width.dim)
    for (i in 1:dim(mat)[1]) {
      for (j in 1:dim(mat)[2]) {
        if((row(mat)[i,j] %% int_col_spac) == 0   & (col(mat)[i,j] %% int_col_spac) == 0) {   #%% returns the remainder of x/y i.e x/2. If true, then even number
          mat[i, j] <- mat[i, j] + 1
        }}}
    mat  #if divisible by intercolonial spacing , then return fecund
    row.col = row(mat)[which(!mat == 0)] #gets the index values for non-zeros
    col.col = col(mat)[which(!mat == 0)]
    # Compute the midpoints
    mid_row = floor(long_dim / 2) + 1
    mid_col = floor(trans_dim / 2) + 1
    #small.mat = mat[int_col_spac:(height.dim + int_col_spac-1), int_col_spac: (width.dim + int_col_spac-1) ]  #this shaves off the edges
    small.mat = mat[min(row.col): max(row.col) , min(col.col): max(col.col)]  #this shaves off the edges
    small.mat
    small.mat_egg = small.mat
    small.mat_egg[] <- 0  #make all values zero
    center <- c(nrow(small.mat) %/% 2 + 1, ncol(small.mat) %/% 2 + 1)
    indices <- which(small.mat > 0, arr.ind = TRUE)
    distances <- sqrt((indices[,1] - center[1])^2 + (indices[,2] - center[2])^2)
    # closest <- which.min(distances)
    # closest_position <- indices[closest, ]
    # small.mat[closest_position[1],closest_position[2]] <- 0
    # small.mat_egg[closest_position[1],closest_position[2]] <- 1
    col.number = sum(small.mat != 0)  #colony number
    col.dens = (height.dim * width.dim) / col.number  #colony density
    # Compute the starting indices, subtract half of small.mat's dimensions from the midpoints
    # row_index = mid_row - floor(nrow(small.mat) / 2)
    # col_index = mid_col - floor(ncol(small.mat) / 2)
    row_index = floor(nrow(spemap) * 3/4) - floor(nrow(small.mat) / 2)
    col_index = mid_col - floor(ncol(small.mat) / 2)
    spemap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1] <- small.mat
    
    #remove a spemap and replace into eggmap 
    row.col.all = row(spemap[,,2])[which(!spemap[,, benthos] == 0)] #finds all row vals
    col.col.all = col(spemap[,,2])[which(!spemap[,, benthos] == 0)]
    spemap[row.col.all[1], col.col.all[1], benthos] <- 0
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    
    spe_config= levelplot(t(apply(spemap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1], 2, rev)), col.regions = coul, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    #eggmap = blank.spemap
    eggmap[row.col.all[1], col.col.all[1], benthos] <- 1
    egg_config= levelplot(t(apply(eggmap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1], 2, rev)), col.regions = coul2, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    lst1 <- list('spemap00'=spemap, 'eggmap00'=eggmap, smallmat = small.mat, plot1 = spe_config, plot2 = egg_config) #store output in list
    return(lst1)
  }
  
  if(config1 == 'patch_center') {   #middle of the patch and more flexible
    #flexible patch and intercolonial density
    no.height =  no.height
    no.width =   no.width
    height.dim = no.height * int_col_spac #+ 1  #creates enough space
    width.dim = no.width * int_col_spac #+ 1
    mat <- matrix(0:0, nrow = height.dim, ncol = width.dim)
    for (i in 1:dim(mat)[1]) {
      for (j in 1:dim(mat)[2]) {
        if((row(mat)[i,j] %% int_col_spac) == 0   & (col(mat)[i,j] %% int_col_spac) == 0) {   #%% returns the remainder of x/y i.e x/2. If true, then even number
          mat[i, j] <- mat[i, j] + 1
        }}}
    mat  #if divisible by intercolonial spacing , then return fecund
    row.col = row(mat)[which(!mat == 0)] #gets the index values for non-zeros
    col.col = col(mat)[which(!mat == 0)]
    # Compute the midpoints
    mid_row = floor(long_dim / 2) + 1
    mid_col = floor(trans_dim / 2) + 1
    #small.mat = mat[int_col_spac:(height.dim + int_col_spac-1), int_col_spac: (width.dim + int_col_spac-1) ]  #this shaves off the edges
    small.mat = mat[min(row.col): max(row.col) , min(col.col): max(col.col)]  #this shaves off the edges
    small.mat
    small.mat_egg = small.mat
    small.mat_egg[] <- 0  #make all values zero

    center <- c(nrow(small.mat) %/% 2 + 1, ncol(small.mat) %/% 2 + 1)
    indices <- which(small.mat > 0, arr.ind = TRUE)
    distances <- sqrt((indices[,1] - center[1])^2 + (indices[,2] - center[2])^2)
    closest <- which.min(distances)
    closest_position <- indices[closest, ]
    small.mat[closest_position[1],closest_position[2]] <- 0
    small.mat_egg[closest_position[1],closest_position[2]] <- 1
    col.number = sum(small.mat != 0)  #colony number
    col.dens = (height.dim * width.dim) / col.number  #colony density
    # Compute the starting indices, subtract half of small.mat's dimensions from the midpoints
    row_index = mid_row - floor(nrow(small.mat) / 2)
    col_index = mid_col - floor(ncol(small.mat) / 2)
    spemap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1] <- small.mat
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    spe_config= levelplot(t(apply(spemap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1], 2, rev)), col.regions = coul, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    eggmap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1] <- small.mat_egg
    egg_config= levelplot(t(apply(eggmap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1], 2, rev)), col.regions = coul2, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    
    lst1 <- list('spemap00'=spemap, 'eggmap00'=eggmap, smallmat = small.mat, plot1 = spe_config, plot2 = egg_config) #store output in list
    return(lst1)
  }
  
  if(config1 == '2021tenuis') {
    no.height =  no.height
    no.width =   no.width
    load('./Rdata/2021_ten_adu_coords_rect.RData') #tenuis
    small.mat = data1
    small.mat_egg = small.mat
    small.mat_egg[] <- 0  #make all values zero
    center <- c(nrow(small.mat) %/% 2 + 1, ncol(small.mat) %/% 2 + 1)
    indices <- which(small.mat > 0, arr.ind = TRUE)
    distances <- sqrt((indices[,1] - center[1])^2 + (indices[,2] - center[2])^2)
    closest <- which.min(distances)
    closest_position <- indices[closest, ]
    small.mat[closest_position[1],closest_position[2]] <- 0
    small.mat_egg[closest_position[1],closest_position[2]] <- 1
    height.dim = nrow(small.mat)
    width.dim = ncol(small.mat)
    col.number = sum(small.mat != 0)  #colony number
    col.dens = (height.dim * width.dim) / col.number  #colony density
    #inset small.at in to mat
    col_index <- floor((trans_dim - ncol(small.mat)) / 2) + 1  # left position
    buffer = 2 #(m)
    row_index <- (long_dim - buffer) - nrow(small.mat)+1  #top position
    spemap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1] <- small.mat
    eggmap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1] <- small.mat_egg
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    spe_config= levelplot(t(apply(spemap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1], 2, rev)), col.regions = coul, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    egg_config= levelplot(t(apply(eggmap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1], 2, rev)), col.regions = coul2, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    lst1 <- list('spemap00'=spemap, 'eggmap00'=eggmap, smallmat = small.mat, plot1 = spe_config, plot2 = egg_config) #store output in list
    return(lst1)
  }
  
  if(config1 == '2021digit') {
    no.height =  no.height
    no.width =   no.width
    load('./Rdata/2021_dig_adu_coords_rect.RData') #tenuis
    small.mat = data1
    small.mat_egg = small.mat
    small.mat_egg[] <- 0  #make all values zero
    center <- c(nrow(small.mat) %/% 2 + 1, ncol(small.mat) %/% 2 + 1)
    indices <- which(small.mat > 0, arr.ind = TRUE)
    distances <- sqrt((indices[,1] - center[1])^2 + (indices[,2] - center[2])^2)
    closest <- which.min(distances)
    closest_position <- indices[closest, ]
    small.mat[4,3] <- 1   #added back lost colony from conversion
    small.mat[closest_position[1],closest_position[2]] <- 0
    small.mat_egg[closest_position[1],closest_position[2]] <- 1
    height.dim = nrow(small.mat)
    width.dim = ncol(small.mat)
    col.number = sum(small.mat != 0)  #colony number
    col.dens = (height.dim * width.dim) / col.number  #colony density
    #inset small.at in to mat
    col_index <- floor((trans_dim - ncol(small.mat)) / 2) + 1  # left position
    buffer = 2 #(m)
    row_index <- (long_dim - buffer) - nrow(small.mat)+1  #top position
    spemap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1] <- small.mat
    eggmap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1] <- small.mat_egg
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    spe_config= levelplot(t(apply(spemap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1], 2, rev)), col.regions = coul, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    egg_config= levelplot(t(apply(eggmap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1], 2, rev)), col.regions = coul2, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    lst1 <- list('spemap00'=spemap, 'eggmap00'=eggmap, smallmat = small.mat, plot1 = spe_config, plot2 = egg_config) #store output in list
    return(lst1)
  }
    
  
  if(config1 == 'petal') {
    #petal
    #int_col_spac = vec.x[loop]
    #Middle
    mid.long = long_dim - int_col_spac - 2
    mid.trans = floor(trans_dim / 2 - 0.5)
    eggmap = blank.spemap
    eggmap[mid.long, mid.trans, benthos] <- 1  #76,76
    eggmap[,,6]
    
    
    #Sperm
    spemap[mid.long + int_col_spac, mid.trans, benthos] <-1 #N   #96, 76
    spemap[mid.long - int_col_spac, mid.trans, benthos] <-1 #S  #56, 76
    spemap[ceiling(mid.long + 0.5 * int_col_spac), ceiling(mid.trans + 17/20 * int_col_spac),  benthos] <- 1  #86, 93     #10 diff, 17
    spemap[ceiling(mid.long + 0.5 * int_col_spac), floor(mid.trans - 17/20 * int_col_spac),  benthos] <- 1
    spemap[floor(mid.long - 0.5 * int_col_spac), ceiling(mid.trans + 17/20 * int_col_spac),  benthos] <- 1
    spemap[floor(mid.long - 0.5 * int_col_spac), floor(mid.trans - 17/20 * int_col_spac),  benthos] <- 1
    lst1 <- list('spemap00'=spemap, 'eggmap00'=eggmap) #store output in list
    return(lst1)
  }
  
  if(config1 == 'poisson') {
    #set.seed(123) # For reproducibility
    
    no.height =  no.height
    no.width =   no.width
    small.mat = matrix(0, nrow = no.height, ncol = no.width) # Initialize matrix with zeros
    #mean_colonies <- 10 # Mean number of colonies, adjust based on your criteria
    #patch_density <- 10 #rpois(1, lambda = mean_colonies) # Number of colonies to place
    
    placed_colonies <- 0
    
    while (placed_colonies < patch_density) {
      x_pos <- sample(1:(no.height), 1) # Random x position, avoiding edges
      y_pos <- sample(1:(no.width), 1) # Random y position, avoiding edges
      
      #Check if the cell is already occupied
      # if (small.mat[x_pos, y_pos] == 0) {
      #   small.mat[x_pos, y_pos] <- 1 # Mark colony presence (for sperm map or egg map as needed)
      #   placed_colonies <- placed_colonies + 1
      # }
      small.mat[x_pos, y_pos] <- small.mat[x_pos, y_pos] + 1
      placed_colonies <- placed_colonies + 1
    
    }
    
    sum(small.mat)
    mid_row = floor(long_dim / 2) + 1
    mid_col = floor(trans_dim / 2) + 1
    row_index = floor(nrow(spemap) * 9/10) - floor(nrow(small.mat) / 2)
    col_index = mid_col - floor(ncol(small.mat) / 2)
    spemap[row_index:(row_index + nrow(small.mat) - 1), col_index:(col_index + ncol(small.mat) - 1), z_dim-1] <- small.mat    
    ones_positions <- which(small.mat == 1, arr.ind = TRUE)
    if(nrow(ones_positions) > 0) { # Ensure there is at least one '1' to remove
      # Select one position randomly
      selected_pos <- ones_positions[sample(nrow(ones_positions), 1), ]
      
      # 2. Update spemap: Set the corresponding position to '0'
      spemap[row_index + selected_pos[1] - 1, col_index + selected_pos[2] - 1, z_dim - 1] <- 0
      
      # 3. Update eggmap: Set the same position to '1'
      eggmap[row_index + selected_pos[1] - 1, col_index + selected_pos[2] - 1, z_dim - 1] <- 1
    }
    
    # Display the sum of '1's in small.mat, spemap, and eggmap for verification
    sum(small.mat)
    sum(spemap[,,z_dim - 1])
    sum(eggmap[,,z_dim - 1])
    
    lst1 <- list('spemap00'=spemap, 'eggmap00'=eggmap,smallmat = small.mat) #store output in list
    return(lst1)
  }
  
  if(config1 == 'uniform') {   #not working great
    no.height = no.height
    no.width = no.width
    small.mat = matrix(0, nrow = no.height, ncol = no.width) # Initialize matrix with zeros
    
    total_colonies = min(patch_density, no.height * no.width)
    
    # Calculate rows and columns for grid sections
    grid_rows = floor(sqrt(total_colonies))
    grid_cols = ceiling(total_colonies / grid_rows)
    
    # Calculate step size for each section
    step_x = floor(no.height / grid_rows)
    step_y = floor(no.width / grid_cols)
    
    # Adjust for any rounding down that occurred during step size calculation
    extra_rows = no.height - (step_x * grid_rows)
    extra_cols = no.width - (step_y * grid_cols)
    
    # Initialize counters
    placed_colonies = 0
    for(i in 1:grid_rows) {
      for(j in 1:grid_cols) {
        if (placed_colonies < total_colonies) {
          # Determine the cell within the section
          x_pos = ((i-1) * step_x + 1) + floor(step_x / 2)
          y_pos = ((j-1) * step_y + 1) + floor(step_y / 2)
          
          # Adjust if we're in the last row/column that compensates for rounding
          if (i == grid_rows && extra_rows > 0) {
            x_pos = no.height - floor(extra_rows / 2)
          }
          if (j == grid_cols && extra_cols > 0) {
            y_pos = no.width - floor(extra_cols / 2)
          }
          
          small.mat[x_pos, y_pos] <- 1
          placed_colonies <- placed_colonies + 1
        }
      }
    }
    
    # Display the sum of '1's in small.mat, spemap, and eggmap for verification
    sum(small.mat)
    sum(spemap[,,z_dim - 1])
    sum(eggmap[,,z_dim - 1])
    
    lst1 <- list('spemap00'=spemap, 'eggmap00'=eggmap,smallmat = small.mat) #store output in list
    return(lst1)
  } 
  
  
  if(config1 == '0.5') {
    #this assumes the 4 x 5 colonies is spread over 3x 3 cells. Some trans cells are uneven
    no.height = 5
    #mid_row = floor(long_dim / 2) + 1
    mid_col = floor(trans_dim / 2) + 1
    #spemap[long_dim-3,mid_col-2,benthos]<- 1
    spemap[long_dim-3,mid_col-1,benthos]<- 1
    spemap[long_dim-3,mid_col,benthos]<- 3
    spemap[long_dim-3,mid_col+1,benthos]<- 1
    #spemap[long_dim-3,mid_col+2,benthos]<- 1
    #spemap[long_dim-4,mid_col-2,benthos]<- 1
    spemap[long_dim-4,mid_col-1,benthos]<- 3
    spemap[long_dim-4,mid_col, benthos]<-  8  #take 1 x 0.5^2 out of 4 x 0.5^2
    spemap[long_dim-4,mid_col+1,benthos]<- 3
    #spemap[long_dim-4,mid_col+2,benthos]<- 1
    #spemap[long_dim-5,mid_col-2,benthos]<- 1
    spemap[long_dim-5,mid_col-1,benthos]<- 1
    spemap[long_dim-5,mid_col,benthos]<- 3
    spemap[long_dim-5,mid_col+1,benthos]<- 1
    #spemap[long_dim-5,mid_col+2,benthos]<- 1
    # spemap[long_dim-6,mid_col-2,benthos]<- 1
    # spemap[long_dim-6,mid_col-1,benthos]<- 1
    # spemap[long_dim-6,mid_col,benthos]<- 1
    # spemap[long_dim-6,mid_col+1,benthos]<- 1
    #spemap[long_dim-6,mid_col+2,benthos]<- 1
    # spemap[long_dim-7,mid_col-1,benthos]<- 1
    # spemap[long_dim-7,mid_col,benthos]<- 1
    # spemap[long_dim-7,mid_col+1,benthos]<- 1
    small.mat = spemap[(long_dim-3): (long_dim-5), (mid_col - 1) : (mid_col + 1), benthos]
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    
    spe_config= levelplot(t(apply(spemap[, , benthos], 2, rev)), col.regions = coul, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    eggmap[long_dim-4, mid_col, benthos]<- 1
    egg_config= levelplot(t(apply(eggmap[, , z_dim-1], 2, rev)), col.regions = coul2, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    lst1 <- list('spemap00'=spemap, 'eggmap00'=eggmap, smallmat = small.mat, plot1 = spe_config, plot2 = egg_config) #store output in list
    return(lst1)
  }
  
  if(config1 == '1.5') {
    #body
    #no.height = 5
    #mid_row = floor(long_dim / 2) + 1
    mid_col = floor(trans_dim / 2) + 1
    spemap[long_dim-3,mid_col-3,benthos]<- 1
    #spemap[long_dim-3,mid_col-1,benthos]<- 1
    spemap[long_dim-3,mid_col-1, benthos]<- 1
    spemap[long_dim-3,mid_col,benthos]<- 1
    spemap[long_dim-3,mid_col+1,benthos]<- 1
    spemap[long_dim-3,mid_col+3,benthos]<- 1
    spemap[long_dim-5,mid_col-3,benthos]<- 1
    #spemap[long_dim-4,mid_col-1,benthos]<- 1
    #spemap[long_dim-4,mid_col, benthos]<-  3/4  #take 1 x 0.5^2 out of 4 x 0.5^2
    spemap[long_dim-5,mid_col-1,benthos]<- 1
    spemap[long_dim-5,mid_col,benthos]<- 1
    spemap[long_dim-5,mid_col+1,benthos]<- 1
    spemap[long_dim-5,mid_col+3,benthos]<- 1
    spemap[long_dim-6,mid_col-3,benthos]<- 1
    spemap[long_dim-6,mid_col-1,benthos]<- 1
    spemap[long_dim-6,mid_col,benthos]<- 0  #removed for egg
    spemap[long_dim-6,mid_col+1,benthos]<- 1
    spemap[long_dim-6,mid_col+3,benthos]<- 1
    spemap[long_dim-7,mid_col-3,benthos]<- 1
    spemap[long_dim-7,mid_col-1,benthos]<- 1
    spemap[long_dim-7,mid_col, benthos]<- 1
    spemap[long_dim-7,mid_col+1,benthos]<- 1
    spemap[long_dim-7,mid_col+3,benthos]<- 1
    spemap[long_dim-9,mid_col-3,benthos]<- 1
    spemap[long_dim-9,mid_col-1,benthos]<- 1
    spemap[long_dim-9,mid_col, benthos]<- 1
    spemap[long_dim-9,mid_col+1,benthos]<- 1
    spemap[long_dim-9,mid_col+3,benthos]<- 1
    small.mat = spemap[(long_dim-3): (long_dim-9), (mid_col - 3) : (mid_col + 3), benthos]
    coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
    coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
    
    spe_config= levelplot(t(apply(spemap[, , benthos], 2, rev)), col.regions = coul, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    eggmap[long_dim-4, mid_col, benthos]<- 1
    egg_config= levelplot(t(apply(eggmap[, , z_dim-1], 2, rev)), col.regions = coul2, xlab = "Transverse",
                          ylab = "Longitudinal", main = "Conc. (cells/m^3)")
    lst1 <- list('spemap00'=spemap, 'eggmap00'=eggmap, smallmat = small.mat, plot1 = spe_config, plot2 = egg_config) #store output in list
    return(lst1)
  }
  
  
  # plot_fun1(eggmap, 6)
  # plot_fun1(spemap, 6)
  
  # #original teo spacing
  # eggmap[127,127,2]<-5.7*10^4  #original code creates a circle around the egg map.
  # spemap[107,127,2]<-5.7*10^4  #20m
  # spemap[147,127,2]<-5.7*10^4 #20m
  # spemap[117,144,2]<-5.7*10^4  #19.72m
  # spemap[137,144,2]<-5.7*10^4  #19.72m
  # spemap[117,110,2]<-5.7*10^4  #19.72m
  # spemap[137,110,2]<-5.7*10^4  #19.72m
  
  # #staggered grid
  # pattern = rep(c(0,fecund), 25)  #creates a pattern
  # pattern.mat  <- array(data = pattern, dim = c(21, 20))  #creates a pattern matrix. Use odd and even to stagger.
  # height.dim = 4
  # width.dim = 3
  # small.mat = pattern.mat[1:height.dim, 1:width.dim]  #subset for small patch
  # spemap[(trans_dim / 2):((trans_dim / 2)+height.dim - 1), (long_dim / 2):((long_dim/2) + width.dim - 1),2] = small.mat  #embed into grid
  # spemap[,,2]
  
  ##5 x 5b 1 int reg grid. (superceeded?)
  # pattern = rep(fecund, 25)  #creates a pattern
  # pattern.mat  <- array(data = pattern, dim = c(21, 20))  #creates a pattern matrix. Use odd and even to stagger.
  # height.dim = 5
  # width.dim = 5
  # small.mat = pattern.mat[1:height.dim, 1:width.dim]
  # left = floor((trans_dim / 2) - ((width.dim-1) / 2) )
  # right = floor((trans_dim / 2) + ((width.dim-1) / 2))
  # top = floor((long_dim / 2) - ((height.dim-1) / 2))
  # bot = floor((long_dim / 2) + ((height.dim-1) / 2))
  # spemap[top:bot, left:right, 2] = small.mat  #embed into grid
  # spemap[(trans_dim / 2), (long_dim / 2) , 2] <-0
  # eggmap[(trans_dim / 2), (long_dim / 2) ,2] <-fecund 
  # spemap[,,2]
  # eggmap[,,2]
  
  
  #old patch code - can delete if code above works
  # if(config1 == 'patch') {
  #   #flexible patch and intercolonial density
  #   no.height =  no.height
  #   no.width =   no.width
  #   height.dim = no.height * int_col_spac #+ 1  #creates enough space
  #   width.dim = no.width * int_col_spac #+ 1
  #   mat <- matrix(0:0, nrow = height.dim, ncol = width.dim)
  #   for (i in 1:dim(mat)[1]) {
  #     for (j in 1:dim(mat)[2]) {
  #       if((row(mat)[i,j] %% int_col_spac) == 0   & (col(mat)[i,j] %% int_col_spac) == 0) {   #%% returns the remainder of x/y i.e x/2. If true, then even number
  #         mat[i, j] <- mat[i, j] + 1
  #       }}}
  #   mat  #if divisible by intercolonial spacing , then return fecund
  #   row.col = row(mat)[which(!mat == 0)] #gets the index values for non-zeros
  #   col.col = col(mat)[which(!mat == 0)]
  #   #small.mat = mat[int_col_spac:(height.dim + int_col_spac-1), int_col_spac: (width.dim + int_col_spac-1) ]  #this shaves off the edges
  #   small.mat = mat[min(row.col): max(row.col) , min(col.col): max(col.col)]  #this shaves off the edges
  #   small.mat
  #   col.number = sum(small.mat != 0)  #colony number
  #   col.dens = (height.dim * width.dim) / col.number  #colony density
  #   #inset small.at in to mat
  #   left = floor((trans_dim / 2) - ((nrow(small.mat)-1) / 2) )
  #   right = floor((trans_dim / 2) + ((nrow(small.mat)-1) / 2))
  #   # top = floor((long_dim / 2) - ((ncol(small.mat)-1) / 2)) centred
  #   #bot = floor((long_dim / 2) + ((ncol(small.mat)-1) / 2))  #centred
  #   buffer = 2 #(m)
  #   top = long_dim -ncol(small.mat)- buffer+1
  #   bot = long_dim - buffer
  #   #spemap = blank.spemap
  #   spemap[top:bot, left:right, z_dim-1] = small.mat  #embed into grid
  #   #remove a spemap and replace into eggmap 
  #   row.col.all = row(spemap[,,2])[which(!spemap[,,benthos] == 0)] #finds all row vals
  #   col.col.all = col(spemap[,,2])[which(!spemap[,,benthos] == 0)]
  #   #row.col.end = row(spemap[,,2])[which(!spemap[,,benthos] == 0)][length(row.col)]  
  #   #col.col.end = col(spemap[,,2])[which(!spemap[,,benthos] == 0)][length(col.col)] 
  #   spemap[median(row.col.all),median(col.col.all),benthos] <- 0
  #   #eggmap = blank.spemap
  #   eggmap[median(row.col.all),median(col.col.all),benthos] <- 1
  #   lst1 <- list('spemap'=spemap, 'eggmap'=eggmap) #store output in list
  #   return(lst1)
  # }
  
  ###poisson approach###
  # coin <- c(0, fecund)
  # rand.col = as.numeric(sample(coin, size = long_dim*trans_dim, replace = T, prob = c(0.9, 0.1)))  #can also add clustering here
  # spemap[, , 2] <-rand.col  #but not pois
  # 
  # library(secr)
  # temppop <- sim.popn (D = 10, expand.grid(x = c(0,30), y = c(0,30)), buffer = 30) #D = desmity/hectare
  # plot(temppop, pch = 1, cex= 1.5)
  ####
}
