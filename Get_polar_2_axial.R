##### TRANSFORM POLAR COORDINATE TO AXIAL COORDINATE:

Get_polar_2_axial = function(params2, params3) {
  vector_x = sin(params2) * cos(params3)
  vector_y = sin(params2) * sin(params3)
  vector_z = cos(params2)
  r = sqrt(vector_x ^ 2 + vector_y ^ 2 + vector_z ^ 2)
  
  
  cbind(sign(vector_z) * vector_x / r,
    sign(vector_z) * vector_y / r,
    abs(vector_z / r))
}

