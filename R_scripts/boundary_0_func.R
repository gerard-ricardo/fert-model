boundary_0 <- function(x, z_dim){
  # Dimensions of the array
  x_dim <- dim(x)[1]
  y_dim <- dim(x)[2]
  
  # Top and bottom z boundaries
  x[ , , 1] <- 0
  x[ , , z_dim] <- 0
  
  # Front and back x boundaries
  x[1, , ] <- 0
  x[x_dim, , ] <- 0
  
  # Left and right y boundaries
  x[ , 1, ] <- 0
  x[ , y_dim, ] <- 0
  
  return(x)
}
