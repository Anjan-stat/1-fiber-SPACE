###Function for adaptive MCMC
Get_Adapt_MCMC <- function(accept_rate, sigma_par, batch_no) {
  if (accept_rate < .44) {
    sigma_par = sigma_par - (1 / batch_no)
    sigma_par
    
  } else if (accept_rate >= .44) {
    sigma_par = sigma_par + (1 / batch_no)
    sigma_par
  }  else{
    sigma_par
  }
  
}

