##### FUNCTION TO CALCULATE INDIVIDUAL BSM MEAN SIGNALS (THAT IS FOR SINGLE VOXEL AND MULTIPLE FIBER.) Check:OK

# params1: f (array of num.fib+1)
# params2: theta
# params3: phi
# params4[1]: S0
# params4[2]: d
# b=bvls[k]
# r=bdir[k, ]

Get_signal_pred = function(params1,params2,params3,params4,b,r){
  ## compute component signals ------- num.fibers anisotropic components
  temp.s = rep(NA,num.fib+1)
  temp.s[1] = exp(-b*params4[2])
  #cat("temp1",temp.s[1],"\n")
  for(fib in 2:(num.fib+1)){
    
    temp.s[fib] = exp(-b*params4[2]*t(r)%*%Get_RARt(c(params2[fib-1],params3[fib-1]))%*%r)
    # cat("comp:",fib,"\n","temp:",temp.s[fib],"\n")
  }
  params4[1]*sum(params1*temp.s)
}

#test
# Get_signal_pred(c(.5,.25,.25),c(pi/4,pi/10),c(pi/4,pi/2),c(100,.7),3000,c(.80,.14,.59))

