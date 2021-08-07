
#target distribution is the posterior distribution you are interested in (ro is the voxel number from master_data)
# params1=f
# params2=th_old[j,]
# params3=ph_old[j,]
# params4= c(S0_propose,rem.pars[2])
# prior_meu=wtsn_meu[] **** it's a matrix for multifiber model.
# ro=j**** row from master data
# sigm=sigm_lik

####Test value
# params1=c(.5,.25,.25)
# params2=c(pi/4,pi/10)
# params3=c(pi/4,pi/2)
# params4= c(100,.7)
# prior_meu=matrix(c(.80,.14,.9,.2,.3,.3,.5,.7),nrow = 2, byrow = TRUE)
# ro=5
# kapa=5
# sigm=.01

# ####Test value
# params1=f_old[j,]
# params2=th_old[j, ]
# params3=ph_old[j, ]
# params4= c(100,rem.pars[j,2])
# prior_meu= wtsn_meu[]
# ro=j
# sigm=sigm_lik_old[j]
# 

Get_target = function(params1,
                      params2,
                      params3,
                      params4,
                      kapa,
                      prior_meu,
                      ro,
                      sigm)
  
{
  
  prior_dir = Get_prior_prod(kapa, params2, params3, prior_meu)
  prior_d = dgamma(params4[2],d_shape,d_rate,log = FALSE)
  # prior calculation: prior is in cartesian coordinates. We input the value of polar. So this is atransformation from Cartesian to polar.
  jacobian = Get_jacob_val(params2, params3) #jacobian for transformation from  cartesian to polar
  log_lik_fun = 0 # WE NEED TO ADD LOGLIKELIHOOD TO CALCULATE THE RATIO. BECAUSE, IF WE TAKE RATIO ITSELF, THE LIKELIHOOD RATIO BECOMES UNDEFINED DUE TO LOW LIKELIHOOD VALUE
  for (k in 1:(n.sig)) {
    log_lik_fun = log_lik_fun - ((
      master_data[ro, k + 5] - Get_signal_pred(params1, params2, params3, params4, bvls[k], bdir[k, ])
    ) ^ 2) / (2 * sigm ^ 2)
    
  }

  
  results = list("dir_prior" = prior_dir, "diff_prior" = prior_d, "jacob"=jacobian, "log_lik"=log_lik_fun,"sum"=(log(prior_dir) + log(prior_d) + log(jacobian) + log_lik_fun))
  return(results)
  # log(prior_dir) + log(prior_d) + log(jacobian) + log_lik_fun
}


# target(f,th_old[j,],ph_old[j,],c(S0_propose,rem.pars[2]),kapa,wtsn_meu[],j,sigm_lik)

# Check function with dummy values:  Get_target(c(.5,.25,.25),c(pi/4,pi/10),c(pi/4,pi/2),c(100,.7),5,matrix(c(.80,.14,.9,.2),nrow = 2, byrow = TRUE) ,2,.01) : OK
