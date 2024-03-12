#dispersal5 test bundle dispersal

#this function allows dispersal3 but treats  boundaries serperately

dispersal5_diff  <- function(x, i,j,k, dx, dy, dz, dt, D_L, D_T, D_Z, shear_surf, long_dim, trans_dim, z_dim,  buoyant=TRUE){  #x will be spemap. buoyant=T is default
  #difference <- NA
  difference <- 0
  
  #h_s = ((z_dim - 1) - k) + 0.5
  #h_s <- ifelse(h_s < 0, 0, h_s)
  # height above the benthos at the current depth level. e.g 6 dim -1 - 5  = 0 (k=5 is benthos). 6 dim -1 - 2
  #drag_coeff <- -0.04031 * k + 0.40806 #assuming linear relationship between depth and drag coef
  
  
  #shear_surf <- flow_rate * van_K_c / (log((h_s * 30) / Ks))
  if (k == 1 || k == z_dim) {
    shear_surf <- 0  # Set shear_surf to 0 for these special cases
  } else {
    #shear_surf <- flow_rate * van_K_c / (log((h_s * 30) / Ks))
    #shear_surf <- sqrt(drag_coeff * flow_rate^2)
    shear_surf <- shear_surf
  }
  
  #using law of wall for U(z) flow rate at each height and using z0/yo from Lentz instead of Baird
  #flow_rate_z = (shear_surf / van_K_c) * log(h_s / 0.05)
  
  
  # Computing the advection term
  #if(!is.na(x[i+1,j,k]) && !is.na(x[i,j,k] && !is.na(x[i-1,j,k]))) {
  # if (i == 1 || i == long_dim) {  #special condition for water surface
  #   dC_dx <- 0 
  # } else {
  #   #dC_dx <- (x[i+1, j, k] - x[i-1, j, k]) / (2 * dx)
  #   dC_dx <- (x[i, j, k] - x[i-1, j, k]) / dx
  # }
  
  if (!buoyant) {    ##buoyant = F
    
    #i
    if (i == 1) {  #special condition for water surface
      d2C_dx2 <- (x[i+1,j,k] - x[i,j,k]) / dx   #sets forward difference rather than central for surface boundary
    } else if (i == long_dim) {  # Special condition for bottom boundary
      d2C_dx2 <- (x[i,j,k] - x[i-1,j,k]) / dx  # Backward difference for bottom boundary
    } else {  #Rest of grid
      d2C_dx2 <- (x[i+1,j,k] - 2*x[i,j,k] + x[i-1,j,k]) / dx^2  #FTCS
    }
    
    #j
    if (j == 1) {  #special condition for water surface
      d2C_dy2 <- (x[i,j+1,k] - x[i,j,k]) / dy   #sets forward difference rather than central for surface boundary
    } else if (j == trans_dim) {  # Special condition for bottom boundary
      d2C_dy2 <- (x[i,j,k] - x[i,j-1,k]) / dy  # Backward difference for bottom boundary
    } else {  #Rest of grid
      d2C_dy2 <- (x[i,j+1,k] - 2*x[i,j,k] + x[i,j-1,k]) / dy^2
    }
    
    #k
    if (k == 1) {  #special condition for water surface
      d2C_dz2 <- (x[i, j, k+1] - x[i, j, k]) / dz   #sets forward difference rather than central for surface boundary
    } else if (k == z_dim) {  # Special condition for bottom boundary
      d2C_dz2 <- (x[i,j,k] - x[i, j, k-1]) / dz  # Backward difference for bottom boundary
    } else { #Rest of grid
      d2C_dz2 <- (x[i,j,k+1] - 2*x[i,j,k] + x[i,j,k-1]) / dz^2
    }
    
    #difference =  dt * (flow_rate_z * dC_dx + D_L * d2C_dx2 + D_T * d2C_dy2 + D_Z * d2C_dz2)  # for non-buoyant case
    difference = dt * (D_L * d2C_dx2 + D_T * d2C_dy2 + D_Z * d2C_dz2)
    
    
  } else {   #buoyant = T
    
    #i
    if (i == 1) {  #special condition for forward layer
      d2C_dx2 <- (x[i+1, j, k] - x[i, j, k]) / dx   #sets forward difference rather than central for surface boundary
    } else if (i == long_dim) {  # Special condition for backwards boundary
      d2C_dx2 <- (x[i, j, k] - x[i-1, j, k]) / dx  # Backward difference for bottom boundary
      #(spemap1_1[160, 1, 1] - spemap1_1[160-1, 1, 1]) / 1
    } else {  #Rest of grid
      d2C_dx2 <- (x[i+1, j,k] - 2*x[i, j,k] + x[i-1, j,k]) / dx^2
    }
    
    #j
    if (j == 1) {  #special condition for left boundary
      d2C_dy2 <- (x[i,j+1,k] - x[i,j,k]) / dy   #sets forward difference rather than central for surface boundary
    } else if (j == trans_dim) {  # Special condition for right boundary
      d2C_dy2 <- (x[i, j, k] - x[i, j-1, k]) / dy  # Backward difference for bottom boundary
    } else { #Rest of grid
      d2C_dy2 <- (x[i,j+1,k] - 2*x[i,j,k] + x[i,j-1,k]) / dy^2
    }
    
    #difference =  dt * (flow_rate_z * dC_dx + D_L * d2C_dx2 + D_T * d2C_dy2)  # for non-buoyant case
    difference = dt * (D_L * d2C_dx2 + D_T * d2C_dy2)
    
  }
  
  #
  if (is.na(difference)) {  #if na, make 0
    difference <- 0 # or some other value
  } else if (abs(difference) < 10^-3) {  #if (abs(difference) < 10^-3) difference <- 0
    difference <- 0
  }
  return(difference)   #spits out the cell new value
}
