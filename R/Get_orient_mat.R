
#orient_summary: finds principal eigen vector for a given theta and phi
Get_orient_mat = function(th, ph) {
  mat = matrix(0, nrow = 3, ncol = 3)
  
  stored_x = sin(th) * cos(ph)
  stored_y = sin(th) * sin(ph)
  stored_z = cos(th)
  mat[1, 1] = sum(stored_x * stored_x)
  mat[1, 2] = sum(stored_x * stored_y)
  mat[1, 3] = sum(stored_x * stored_z)
  
  mat[2, 1] = sum(stored_y * stored_x)
  mat[2, 2] = sum(stored_y * stored_y)
  mat[2, 3] = sum(stored_y * stored_z)
  
  mat[3, 1] = sum(stored_z * stored_x)
  mat[3, 2] = sum(stored_z * stored_y)
  mat[3, 3] = sum(stored_z * stored_z)
  return(mat)
}

