
#From all fibers within a voxel, this function Selects most likely fiber aligned to a certain direction
##Input: 
  #(a)fibers:array of theta and phi representing array of fibers.
  #(b)direction: The direction towards we seek the alignment

##Output:
  #(a):(theta,phi) array of one theta and one phi, representing one fiber which is aligned to the "direction"

Get_select_fiber = function(fibers, direction) {
  

  fibers=matrix(unlist(fibers),nrow = num.fib,byrow = T)
  val = c()
  for (i in 1:nrow(fibers)) {
    
    val = c(val, Get_polar_2_axial((fibers[i,1]), (fibers[i,2])) %*%
              direction)
    
  }
  j <- sample(which(abs(val) == max(abs(val))), 1)
  return (fibers[j,])
}

