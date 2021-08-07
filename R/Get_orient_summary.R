
#orient_summary: finds principal eigen vector for a given theta and phi
Get_orient_summary = function(dir_mat) {
  ev <- eigen(dir_mat)
  k <- sample(which(abs(ev$values) == max(abs(ev$values))), 1)
  ev$vectors[, k] / sqrt(sum(ev$vectors[, k] ^ 2))
  
}
