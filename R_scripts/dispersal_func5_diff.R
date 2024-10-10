dispersal5_diff  <- function(x, i,j,k, dx, dy, dz, dt, D_L, D_T, D_Z, shear_surf, long_dim, trans_dim, z_dim,  buoyant=TRUE){
  difference <- 0
  if (!buoyant) {
    if (i == 1) {
      d2C_dx2 <- (x[i+1,j,k] - x[i,j,k]) / dx
    } else if (i == long_dim) {
      d2C_dx2 <- (x[i,j,k] - x[i-1,j,k]) / dx
    } else {
      d2C_dx2 <- (x[i+1,j,k] - 2*x[i,j,k] + x[i-1,j,k]) / dx^2
    }
    if (j == 1) {
      d2C_dy2 <- (x[i,j+1,k] - x[i,j,k]) / dy
    } else if (j == trans_dim) {
      d2C_dy2 <- (x[i,j,k] - x[i,j-1,k]) / dy
    } else {
      d2C_dy2 <- (x[i,j+1,k] - 2*x[i,j,k] + x[i,j-1,k]) / dy^2
    }
    if (k == 1) {
      d2C_dz2 <- (x[i, j, k+1] - x[i, j, k]) / dz
    } else if (k == z_dim) {
      d2C_dz2 <- (x[i,j,k] - x[i, j, k-1]) / dz
    } else {
      d2C_dz2 <- (x[i,j,k+1] - 2*x[i,j,k] + x[i,j,k-1]) / dz^2
    }
    difference = dt * (D_L * d2C_dx2 + D_T * d2C_dy2 + D_Z * d2C_dz2)
  } else {
    if (i == 1) {
      d2C_dx2 <- (x[i+1, j, k] - x[i, j, k]) / dx
    } else if (i == long_dim) {
      d2C_dx2 <- (x[i, j, k] - x[i-1, j, k]) / dx
    } else {
      d2C_dx2 <- (x[i+1, j,k] - 2*x[i, j,k] + x[i-1, j,k]) / dx^2
    }
    if (j == 1) {
      d2C_dy2 <- (x[i,j+1,k] - x[i,j,k]) / dy
    } else if (j == trans_dim) {
      d2C_dy2 <- (x[i, j, k] - x[i, j-1, k]) / dy
    } else {
      d2C_dy2 <- (x[i,j+1,k] - 2*x[i,j,k] + x[i,j-1,k]) / dy^2
    }
    difference = dt * (D_L * d2C_dx2 + D_T * d2C_dy2)
  }
  if (is.na(difference)) {
    difference <- 0
  } else if (abs(difference) < 10^-3) {
    difference <- 0
  }
  return(difference)
}
