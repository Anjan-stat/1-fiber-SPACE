##### FIXING RANGE OF THETA AND PHI by taking the mod of the angles:
Get_chck_th_ph = function(th, ph) {
  x = sin(th) * cos(ph)
  y = sin(th) * sin(ph)
  z = cos(th)
  
  x = sign(z) * x
  y = sign(z) * y
  th_out = acos(abs(z))
  
  if ((x > 0) & (y > 0)) {
    ph_out = acos(x / sin(th_out))
    
  } else if ((x < 0) & (y > 0)) {
    ph_out = acos(x / sin(th_out))
  } else if ((x < 0) & (y < 0)) {
    ph_out = atan(y / x) + pi
  } else{
    ph_out = asin(y / sin(th_out)) %% (2 * pi)
  }
  return (c(th_out, ph_out))
}

