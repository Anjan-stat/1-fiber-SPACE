
Get_axial_2_polar = function(x, y, z) {
  r = sqrt(x ^ 2 + y ^ 2 + z ^ 2)
  x = sign(z) * x / r
  y = sign(z) * y / r
  z = abs(z / r)
  
  th = acos(z)
  ph = atan(y / x)
  if (is.nan(ph))
  {
    ph = runif(1, min = .001, max = 2 * pi)
  } else{
    if ((x > 0) & (y > 0)) {
      ph = acos(x / sin(th))
      
    } else if ((x < 0) & (y > 0)) {
      ph = acos(x / sin(th))
    } else if ((x < 0) & (y < 0)) {
      ph = atan(y / x) + pi
    } else{
      ph = asin(y / sin(th)) %% (2 * pi)
    }
  }
  return (c(th, ph))
}
