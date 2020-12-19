
#Jacobian multiplication for each fibers in a voxel (Geenralized)
Get_jacob_val = function(params2, params3) {
  # params2=c(pi/8,pi/2)
  # params3=c(pi/9,pi/2)
  temp_j = 1
  for (fib in 1:(num.fib)) {
    temp_j = temp_j * abs(sin(params2[fib]))
  }
  return(temp_j)
}

