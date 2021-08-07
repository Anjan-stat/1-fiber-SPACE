rm(list = ls())
library("Directional")
library("erer")
library(pkgmaker)
library(MCMCprecision) #To use rdirichlet function for sampling volume fraction

setwd("C:/Single_fiber")
source("Get_polar_2_axial.R")
source("Get_chck_th_ph.R")
source("Get_axial_2_polar.R")
source("Get_orient_mat.R")
source("Get_orient_summary.R")
source("Get_RARt.R")
source("Get_signal_pred.R")
source("Get_prior_prod.R")
source("Get_jacob_val.R")
source("Get_target.R")
source("Get_check_voxel.R")
source("Get_meu_fun.R")
source("Get_Adapt_MCMC.R")
source("Get_select_fiber.R")
source("Get_norm_vec.R")

setwd("C:/1_Non_system_files/UNLV/Semesters/10_Fall_2020/Neuroimage/Code/Single_fiber")
start_time <- Sys.time()
batch = 100
itrn = 100 ## number of iteration
burn_in = 20
track_voxl = 14
num.fib = 1 #Max Number of nerves possible
SNR=5

original_data = read.csv(file = "Ground_truth/Simulation_source_1fib_connected.csv")
original_data = original_data[order(original_data$Z_axis,
                                    original_data$Y_axis,
                                    original_data$X_axis), ]

#######################################################################################################################################
#######################~~~~~~~READ THE MASTER DATA AND UPDATE THE VARIABLES RELATED TO IT~~~~~~~~######################################
#######################################################################################################################################


master_data = read.csv(file = paste("Input/Input_file_mcmc_3X3_V2_SNR_",SNR,"_connected.csv", sep = ""),sep = ",")
master_data = master_data[order(master_data$Z_axis, master_data$Y_axis, master_data$X_axis), ]

bvls = scan("Input/hardi-scheme.bval.txt")
bdir = scan("Input/hardi-scheme.bvec.txt")
bdir = matrix(bdir, ncol = 3)
n.sig = nrow(bdir) # Number of signal available for each voxel: currently 65
voxl_cnt = nrow(master_data) #Number of voxels in master data: currently 27
# |||CREATE||| a Vox_cube is the cube of dimension (max_vox+1) X (max_vox+1) X (max_vox+1)

max_vox = max(master_data$X_axis,master_data$Y_axis,master_data$Z_axis) #We are going to create a (max_vox+1)^3 sized cube.
vox_cube_old <- array(3, dim = c(max_vox + 1, max_vox + 1, max_vox + 1, 3*num.fib+2))
vox_cube <- vox_cube_old 
dir_cnt <- matrix(0, voxl_cnt, 13)
# colnames(dir_cnt) <- c("1st", "2nd", "3th", "4th", "5th", "6th")


#######################################################################################################################################
# ##################################~~~~~~~MCMC PARAMETER INITIALIZATION~~~~~~~~########################################################
#######################################################################################################################################
##Fix the hyperparmaeter for the parameters:
# f_hyp = c(1,1,1) #Hyperparamter corresponding to Dirichlet distribution for fractional volume
d_shape = 1 # Shape paramter for gamma distribution for adc
d_rate=1 # Rate paramter for gamma distribution for adc
kapa = 200 #Prior for Watson distibution
signal_est = rep(NA,n.sig)
sig_shp=1
sig_rt=.001


#Create the parameters which will be updated.
th_old = matrix(NA, voxl_cnt, num.fib)
ph_old = matrix(NA, voxl_cnt, num.fib)
rem.pars = matrix(NA, voxl_cnt, 2)
f_old = matrix(NA, voxl_cnt, num.fib)
sigm_lik_old=matrix(NA, voxl_cnt,1)


#Initialize the parameters
# th_old[, 1] = runif(1,0,pi/2)
# ph_old[, 1] = runif(1,0,2*pi)
# rem.pars[,1] = mean(master_data$S1) #params4# (S0, d)
# rem.pars[,2] = rgamma(1, shape = d_shape, rate = d_rate) #params4# (S0, d)
f_old[] = .1 #params1#: depends on num.fib
# f_old = cbind(1-f_old[],f_old[])
# sigm_lik_old[] = rem.pars[1,1]/SNR # Where Signal follows N(Mean of ball and stick model. sigma_lik)

th_old[, 1] = original_data$theta_1
ph_old[, 1] = original_data$phi_1
rem.pars[,1] = original_data$S_0 #params4# (S0, d)
rem.pars[,2] = original_data$d #params4# (S0, d)
# f_old[] = original_data$f #params1#: depends on num.fib
f_old = cbind(1-f_old[],f_old[])
sigm_lik_old[] = rem.pars[1,1]/SNR # Where Signal follows N(Mean of ball and stick model. sigma_lik)

track_val = c(th_old[track_voxl,],ph_old[track_voxl,],rem.pars[track_voxl,],f_old[track_voxl,2],sigm_lik_old[track_voxl],Get_meu_fun(master_data$X_axis[track_voxl],master_data$Y_axis[track_voxl],master_data$Z_axis[track_voxl],vox_cube_old[,,,1:2])$out)


##Log Variance for sampling Proposal value
lsig_theta1 = -3
lsig_phi1 = -3
lsig_s0 = -.2
lsig_d = -1
lsig_f=-1
lsig_sigm_lik=-1


#Create updating placeholders.
th_propose_1 = 0
ph_propose_1 = 0
th_propose_2 = 0
ph_propose_2 = 0



### VARIABLES TO STORE PARAMETER SAMPLES OF ALL ITERATIONS WITHIN 1 BATCH FOR ALL VOXEL
strd_itrn_th = matrix(NA, nrow = voxl_cnt, ncol = itrn) 
strd_itrn_ph = matrix(NA, nrow = voxl_cnt, ncol = itrn)
strd_itrn_S0= matrix(NA, nrow = voxl_cnt, ncol = itrn)
strd_itrn_d=matrix(NA, nrow = voxl_cnt, ncol = itrn)
strd_itrn_f=matrix(NA, nrow = voxl_cnt, ncol = itrn)


##### strd_batch_ Keeps addding the values from each batch (will discard the burn-in values)----so that we can average them in the end ##########################
strd_batch_1 = array(0, dim = c(voxl_cnt, 3, 3))
strd_batch_2 = matrix(0, nrow = voxl_cnt, ncol = num.fib+2)
# strd_batch_3 = array(0,dim = c(voxl_cnt,1))
# strd_batch_4 = array(0,dim = c(voxl_cnt,1))
# strd_batch_5 = array(0, dim = c(voxl_cnt,num.fib))
#####Required to convert the array of the theta and phi into final angles##########################


###########################################################################################################################################
for (cnt in 1:voxl_cnt) {
  vox_cube[master_data$X_axis[cnt], master_data$Y_axis[cnt], master_data$Z_axis[cnt], 1] = th_old[cnt,]
  vox_cube[master_data$X_axis[cnt], master_data$Y_axis[cnt], master_data$Z_axis[cnt], 2] = ph_old[cnt,]
  vox_cube[master_data$X_axis[cnt], master_data$Y_axis[cnt], master_data$Z_axis[cnt], 3:4] = rem.pars[cnt,]
  vox_cube[master_data$X_axis[cnt], master_data$Y_axis[cnt], master_data$Z_axis[cnt], 5]= f_old[cnt,1]
}

vox_cube_old = vox_cube
########################### Loop starts ####################################################################################################


#Loop for h(batch)----> i(iteration)-------j(voxel)

for (h in 1:batch) {
  acc_rate_theta1 = 0 #INITIALIZING THE ACCEPTANCE RATE FOR THETA. OUTPUT
  acc_rate_phi1 = 0 #INITIALIZING THE ACCEPTANCE RATE FOR PHI. OUTPUT
  acc_rate_s0 = 0 #INITIALIZING THE ACCEPTANCE RATE FOR S0. OUTPUT
  acc_rate_d = 0 #INITIALIZING THE ACCEPTANCE RATE FOR d OUTPUT
  acc_rate_f = 0 #INITIALIZING THE ACCEPTANCE RATE FOR f OUTPUT
  acc_rate_sigm_lik =0 #INITIALIZING THE ACCEPTANCE RATE FOR likelihood sigma OUTPUT

  
  for (i in 1:itrn) {
    vox_cube_old = vox_cube ##This step is necessary as we want to update the prior based on whatever theta and phi we calculate from previous iteration.
    
    for (j in 1:voxl_cnt) {
      ### We have to Calculate Meu beforehand for a particular iteration. Why? Because previous voxel's updated value should not affect the next voxel's prior.
      ### CAUTION: Here we are indexing the voxels by master data coordinates. So Master data contains the initial values.
      ########but the theta and phi will be updated within the voxel cube.
      ### REMEDY: We are taking the theta and phi values from Vox_cube which has a fixed theta phi when iteration number is fixed.
      #####Vox_cube is getting updated within each iteration as we proceed with each voxel. but vox_cube_old is fixed for a iteration.
      
      wtsn_meu = Get_meu_fun(master_data$X_axis[j],
                             master_data$Y_axis[j],
                             master_data$Z_axis[j],
                             vox_cube_old[,,,1:2])$out
      dir_cnt[j, ] = dir_cnt[j, ] + Get_meu_fun(master_data$X_axis[j],
                                                master_data$Y_axis[j],
                                                master_data$Z_axis[j],
                                                vox_cube_old[,,,1:2])$cnt
      
      
      ### PARAMETER UPDATES:
      
      ########################################~~~~f~~~~########################################
      
      
      f_propose = rnorm(1,f_old[j,(num.fib+1)],exp(lsig_f)) ##******** what should be the value of the standard deviation? ***************##

      if((f_propose<0)||(f_propose>1)){
        accept=0
      } else {
      f_propose=c(1-f_propose,f_propose)
      
      ratio = min(1, exp(
      Get_target(f_propose, th_old[j, ],ph_old[j, ],rem.pars[j,], kapa,wtsn_meu[],j,sigm_lik_old[j])$sum - 
        Get_target(f_old[j,],th_old[j, ], ph_old[j, ],rem.pars[j,],kapa,wtsn_meu[],j,sigm_lik_old[j])$sum
                        )
                  )
      accept = (runif(1) < ratio)
      f_old[j,] = ifelse(c(accept,accept), f_propose, f_old[j,])
      }
      acc_rate_f = acc_rate_f + (j == track_voxl) * accept
      
      # if (j == track_voxl){
      #   print(c("f",i,ratio,accept))
      # }
      ########################################~~~~f~~~~########################################
 
      
      ########################################~~~~S0~~~~########################################
      # S0_propose = rnorm(1, rem.pars[j,1], exp(2 * lsig_s0)) ##******** what should be the value of the standard deviation? ***************##
      # 
      # if(S0_propose<0){
      #   accept=0
      # } else {
      # ratio = min(1, exp(Get_target(f_old[j,],th_old[j, ],ph_old[j, ],c(S0_propose, rem.pars[j,2]),kapa,wtsn_meu[],j,sigm_lik_old[j])$sum - 
      #                      Get_target(f_old[j,], th_old[j, ], ph_old[j, ], rem.pars[j,], kapa, wtsn_meu[], j, sigm_lik_old[j])$sum
      #                   )
      #             )
      # accept = (runif(1) < ratio)
      # rem.pars[j,1] = ifelse(accept, S0_propose, rem.pars[j,1])
      #         }
      # 
      # acc_rate_s0 = acc_rate_s0 + (j == track_voxl) * accept
      
      # if (j == track_voxl){
      #   print(c("S0",i,ratio,accept))
      # }
      
      ########################################~~~~S0~~~~########################################
      
      ########################################~~~~d~~~~########################################
      
      # d_propose = rnorm(1, rem.pars[j,2], exp(2 * lsig_d))
      # if(d_propose<0){
      #   accept=0
      # } else {
      # ratio = min(1, exp(
      #   Get_target(f_old[j,],th_old[j, ],ph_old[j, ],c(rem.pars[j,1], d_propose),kapa,wtsn_meu[],j,sigm_lik_old[j])$sum -
      #     Get_target(f_old[j,], th_old[j, ], ph_old[j, ], rem.pars[j,], kapa, wtsn_meu[], j, sigm_lik_old[j])$sum
      #                   )
      #             )
      # accept = (runif(1) < ratio)
      # rem.pars[j,2] = ifelse(accept, d_propose, rem.pars[j,2])
      # }
      # 
      # acc_rate_d = acc_rate_d++(j == track_voxl) * accept
      
      # if (j == track_voxl){
      #   print(c("d",i,ratio,accept))
      # }
      
      
      ########################################~~~~d~~~~########################################
      
      
      ########################################~~~~theta~~~~########################################
      
      # th_propose = rnorm(1,th_old[j, ],exp(lsig_theta1))
      # 
      # if((th_propose<0)||(th_propose>pi/2)){
      #   accept=0
      # } else {
      # 
      # ratio = min(1, exp(
      #   Get_target(f_old[j,],
      #     th_propose,
      #     ph_old[j, ],
      #     rem.pars[j,],
      #     kapa,
      #     wtsn_meu[],
      #     j,
      #     sigm_lik_old[j]
      #   )$sum - Get_target(f_old[j,], th_old[j, ], ph_old[j, ], rem.pars[j,], kapa, wtsn_meu[], j, sigm_lik_old[j])$sum
      # ))
      # accept = (runif(1) < ratio)
      # th_old[j, ][1] = ifelse(accept, th_propose[1], th_old[j, ][1])
      # }
      # acc_rate_theta1 = acc_rate_theta1++(j == track_voxl) * accept
      
      # if (j == track_voxl){
      #   print(c("theta",i,ratio,accept))
      # }
      ########################################~~~~theta~~~~########################################
      
      
      
      ########################################~~~~phi~~~~########################################
      
      # ph_propose = rnorm(1,ph_old[j, ],exp(lsig_phi1))
      # 
      # if((ph_propose<0)||(ph_propose>pi*2)){
      #   accept=0
      # } else {
      #   
      # 
      # ratio = min(1, exp(
      #   Get_target(
      #     f_old[j,],
      #     th_old[j, ],
      #     ph_propose,
      #     rem.pars[j,],
      #     kapa,
      #     wtsn_meu[],
      #     j,
      #     sigm_lik_old[j]
      #   )$sum - Get_target(
      #     f_old[j,],
      #     th_propose,
      #     ph_old[j, ],
      #     rem.pars[j,],
      #     kapa,
      #     wtsn_meu[],
      #     j,
      #     sigm_lik_old[j]
      #   )$sum
      # ))
      # accept = (runif(1) < ratio)
      # ph_old[j, ][1] = ifelse(accept, ph_propose[1], ph_old[j, ][1])
      # 
      # }
      # acc_rate_phi1 = acc_rate_phi1++(j == track_voxl) * accept
      
      # if (j == track_voxl){
      #   print(c("phi",i,ratio,accept))
      # }
      # 
      ########################################~~~~phi~~~~########################################
      
      
       
      ########################################~~~~Sigma of Likelihood~~~~########################################
      
      # for (k in 1:n.sig){
      #   signal_est[k] = Get_signal_pred(params1= f_old[j,],params2=th_old[j, ],params3=ph_old[j, ],params4=rem.pars[j,],bvls[k],bdir[k,])
      # }
      # 
      # sigm_lik_old[j] = sqrt(1/rgamma(1,sig_shp+n.sig/2,sig_rt+1/2*sum((master_data[j,6:ncol(master_data)]-signal_est)^2))) 
      # 
      ########################################~~~~Sigma of Likelihood~~~~########################################
      
      
      #######################################################################################
      # vox_cube[master_data$X_axis[j], master_data$Y_axis[j], master_data$Z_axis[j], 1] = th_old[j, ][1]
      # vox_cube[master_data$X_axis[j], master_data$Y_axis[j], master_data$Z_axis[j], 2] = ph_old[j, ][1]
      # 
      
      # #######################################################################################
      # #######################################################################################
  
      # Each direction is Required for Watson direction calculation and restare for saving the values for future updation of code based on volume fraction
      vox_cube[master_data$X_axis[j], master_data$Y_axis[j], master_data$Z_axis[j], 1] = th_old[j,][1]
      vox_cube[master_data$X_axis[j], master_data$Y_axis[j], master_data$Z_axis[j], 2] = ph_old[j,][1]
      vox_cube[master_data$X_axis[j], master_data$Y_axis[j], master_data$Z_axis[j], 3:4] = rem.pars[j,]
      vox_cube[master_data$X_axis[j], master_data$Y_axis[j], master_data$Z_axis[j], 5]= f_old[j,2]
      #######################################################################################
      
      
      
      strd_itrn_th[j, i] = th_old[j, ][1]
      strd_itrn_ph[j, i] = ph_old[j, ][1]
      strd_itrn_S0[j, i] = rem.pars[j,1] 
      strd_itrn_d[j, i] = rem.pars[j,2]
      strd_itrn_f[j, i] = f_old[j,2]
      
      
      
      

      # cat("batch:",h,"\n","iteration:",i,"\n","voxels:",j)
    }#Loop ends for j---------> voxel
    
    track_mu=Get_meu_fun(master_data$X_axis[track_voxl], master_data$Y_axis[track_voxl],master_data$Z_axis[track_voxl],vox_cube_old[,,,1:2])$out
    track_val = rbind(track_val, c(th_old[track_voxl, ], ph_old[track_voxl, ], rem.pars[track_voxl,],f_old[track_voxl,2],sigm_lik_old[track_voxl],track_mu))
    
  }#Loop ends for i---------> iteration
  
  ####### AFTER EACH BATCH:  
  ################### 1. update the MCMC sigma based on acceptance rate 
  
  lsig_theta1 = Get_Adapt_MCMC(acc_rate_theta1 / itrn, lsig_theta1, h)

  lsig_phi1 = Get_Adapt_MCMC(acc_rate_phi1 / itrn, lsig_phi1, h)

  lsig_s0 = Get_Adapt_MCMC(acc_rate_s0 / itrn, lsig_s0, h)
  lsig_d = Get_Adapt_MCMC(acc_rate_d / itrn, lsig_d, h)
  lsig_f=Get_Adapt_MCMC(acc_rate_f/itrn, lsig_f, h)
  # lsig_sigm_lik= Get_Adapt_MCMC(acc_rate_sigm_lik/itrn, lsig_sigm_lik, h)
  
  ##################### 2. update the values of the parameters based on the iterated samples.
  for (p in 1:voxl_cnt) {
    strd_batch_1[p, , ] = strd_batch_1[p, , ] + Get_orient_mat(strd_itrn_th[p, ], strd_itrn_ph[p, ])
    # strd_batch_2[p, ] = strd_batch_2[p, ] +   cbind(rowMeans(strd_itrn_S0), rowMeans(strd_itrn_d), rowMeans(strd_itrn_f) )[p,]
  }
  strd_batch_2 = strd_batch_2 +   cbind(rowMeans(strd_itrn_S0), rowMeans(strd_itrn_d), rowMeans(strd_itrn_f) )
  dir_cnt = dir_cnt * (h > burn_in)
  strd_batch_1 = strd_batch_1 * (h > burn_in)
  strd_batch_2 = strd_batch_2 * (h > burn_in)

  }

#Loop ends for h---------> batch



for (p in 1:voxl_cnt) {
  fiber_1 = Get_orient_summary(strd_batch_1[p, , ])
  th_old[p, ][1] = Get_axial_2_polar(fiber_1[1], fiber_1[2], fiber_1[3])[1]
  ph_old[p, ][1] = Get_axial_2_polar(fiber_1[1], fiber_1[2], fiber_1[3])[2]
  strd_batch_2[p,]=strd_batch_2[p,]/(batch-burn_in)

}


final_tbl = cbind(original_data[, 4:5],
                  th_old[, 1],
                  ph_old[, 1])

colnames(final_tbl)=c("original_theta1", "original_phi1", "estimated_theta1", "estimated_phi1")
# (if(h >= batch) break("complete"))
head(final_tbl)





end_time <- Sys.time()
tim = end_time - start_time
sigm_info = cbind(
  c(lsig_theta1, acc_rate_theta1 / itrn),
  c(lsig_phi1, acc_rate_phi1 / itrn),
  c(lsig_s0, acc_rate_s0 / itrn),
  c(lsig_d, acc_rate_d / itrn),
  c(lsig_f, acc_rate_f / itrn),
  c(lsig_sigm_lik,acc_rate_sigm_lik/itrn)
)
colnames(sigm_info) = c("theta1", "phi1", "S0", "d","f", "lik_sigma")
rownames(sigm_info) = c("s.d.", "Acceptance rate")
# final_tbl=cbind(final_tbl,original_mat, estimated_mat, dir_cnt,log_lik)
final_tbl = cbind(final_tbl, dir_cnt)
final_output = list(final_tbl, sigm_info, tim)

write.list(
  final_output,
  file = paste(
    "Test_1202/",
    "Fib1_all_random_batch",
    batch,
    "Itrn_",
    i,
    "_SNR",
    SNR,
    "_kapa_",
    kapa,
    "Output_connected.csv",
    sep = ""
  ),
  t.name = c("Fiber_direction", "MCMC info", "Time"),
  row.names = T
)

colnames(track_val)=c("theta","phi","S0","d",	"f","sigma","meu_x","meu_y","meu_z")
path <- paste0("Test_1202/","Fib1_all_random_batch",batch,"Itrn_",i,"_SNR_",SNR,"_kapa_",kapa,"Chain_2_connected.csv", sep = "")


write.table(track_val,file = path,row.names = F, sep = ",")


#########################################################################################################################################

##############################################################Code for Analysis##########################################################

#########################################################################################################################################

test_data = read.csv(file = path ,sep = ",")
head(test_data)

Dir_Prior=matrix(NA,nrow(test_data),1)
Diffusivity_Prior=matrix(NA,nrow(test_data),1)
Jacobian_Val=matrix(NA,nrow(test_data),1)
Log_Likelihood=matrix(NA,nrow(test_data),1)
Total_Val=matrix(NA,nrow(test_data),1)


for (m in 1:nrow(test_data)){
  wt_mu=c(test_data$meu_x[m],test_data$meu_y[m],test_data$meu_z[m])
  Values= Get_target(params1 =test_data$f[m],params2=test_data$theta[m],params3=test_data$phi[m],params4=c(test_data$S0[m],test_data$d[m]),kapa = 10,prior_meu =wt_mu,ro=14,sigm =test_data$sigma[m])
  
  Dir_Prior[m,]= Values$dir_prior
  Diffusivity_Prior[m,]= Values$diff_prior
  Jacobian_Val[m,]= Values$jacob
  Log_Likelihood[m,]= Values$log_lik
  Total_Val[m,]= Values$sum
  
}
test_data=cbind(test_data,Dir_Prior,Diffusivity_Prior,Jacobian_Val,Log_Likelihood,Total_Val)
head(test_data)

write.table(test_data,file = paste("Test_1202/Test_data_SNR_",SNR,"_kapa_", kapa, ".csv",sep = ""),row.names = F,sep = ",")














