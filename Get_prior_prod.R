
#Multiplication of prior value and fibers direction (Generalized)
Get_prior_prod = function(kapa, params2, params3, prior_meu) {
  #
  ####Test value
  # params1=f_old[j,]
  # params2=th_old[j, ]
  # params3=ph_old[j, ]
  # params4= rem.pars[j,]
  # prior_meu= wtsn_meu[]
  # ro=j
  # kapa=0
  # sigm=.01

  
  temp_p = 1
  for (fib in 1:(num.fib)) {
    temp_p = temp_p * exp(kapa * ((
      ####For multifiber
      # prior_meu[fib, ] * Get_polar_2_axial(params2[fib], params3[fib])
      #####End
      as.vector(prior_meu) %*% as.vector(Get_polar_2_axial(params2[fib], params3[fib]))
    ) ^ 2))
  }
  temp_p
}
