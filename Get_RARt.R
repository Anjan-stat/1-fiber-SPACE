############# DEFINE BSM FUNCTIONS ################# (OK)

Get_RARt = function(pars)
{
  th = pars[1]
  ph = pars[2]
  matrix(
    c(
      sin(th) ^ 2 * cos(ph) ^ 2,
      sin(th) ^ 2 * cos(ph) * sin(ph),
      sin(th) * cos(th) * cos(ph),
      sin(th) ^ 2 * cos(ph) * sin(ph),
      sin(th) ^ 2 * sin(ph) ^ 2,
      sin(th) * cos(th) * sin(ph),
      cos(th) * sin(th) * cos(ph),
      sin(th) * cos(th) * sin(ph),
      cos(th) ^ 2
    ),
    ncol = 3,
    byrow = T
  )
}

