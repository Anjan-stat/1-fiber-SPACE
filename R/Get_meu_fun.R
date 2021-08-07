
#Calculating the prior from surrounding information

Get_meu_fun = function(crd_i, crd_j, crd_k, vox_array) {
  
  
  # crd_i= master_data$X_axis[j]
  # crd_j= master_data$Y_axis[j]
  # crd_k=master_data$Z_axis[j]
  # vox_array=vox_cube_old[,,,1:2]
  # 
  
  ####From Upper Layer(UL) Left Corner(LC) Front Side(FS) --------- Lower Layer(LL), Right Corner(RC), Back Side (BS) 
  v_1 = Get_check_voxel((crd_i - 1), (crd_j - 1), (crd_k + 1), vox_array)
  a_1 = Get_check_voxel((crd_i + 1), (crd_j + 1), (crd_k - 1), vox_array)
  p_1 =  Get_polar_2_axial(pi / 4, 5 * pi / 4)
  
  ####From Upper Layer(UL) Left Corner(LC) Front Side(FS) --------- Lower Layer(LL), Right Corner(RC), Back Side (BS)
  v_2 = Get_check_voxel((crd_i), (crd_j - 1), (crd_k + 1), vox_array)
  a_2 = Get_check_voxel((crd_i), (crd_j + 1), (crd_k - 1), vox_array)
  p_2 =  Get_polar_2_axial(pi / 4, 6 * pi / 4)
  
  v_3 = Get_check_voxel((crd_i + 1), (crd_j - 1), (crd_k + 1), vox_array)
  a_3 = Get_check_voxel((crd_i - 1), (crd_j + 1), (crd_k - 1), vox_array)
  p_3 =  Get_polar_2_axial(pi / 4, 7 * pi / 4)
  
  v_4 = Get_check_voxel((crd_i + 1), (crd_j), (crd_k + 1), vox_array)
  a_4 = Get_check_voxel((crd_i - 1), (crd_j), (crd_k - 1), vox_array)
  p_4 =  Get_polar_2_axial(pi / 4, 0)
  
  
  v_5 = Get_check_voxel((crd_i + 1), (crd_j + 1), (crd_k + 1), vox_array)
  a_5 = Get_check_voxel((crd_i - 1), (crd_j - 1), (crd_k - 1), vox_array)
  p_5 =  Get_polar_2_axial(pi / 4, pi / 4)
  
  v_6 = Get_check_voxel((crd_i), (crd_j + 1), (crd_k + 1), vox_array)
  a_6 = Get_check_voxel((crd_i), (crd_j - 1), (crd_k - 1), vox_array)
  p_6 =  Get_polar_2_axial(pi / 4, 2 * pi / 4)
  
  v_7 = Get_check_voxel((crd_i - 1), (crd_j + 1), (crd_k + 1), vox_array)
  a_7 = Get_check_voxel((crd_i + 1), (crd_j - 1), (crd_k - 1), vox_array)
  p_7 =  Get_polar_2_axial(pi / 4, 3 * pi / 4)
  
  v_8 = Get_check_voxel((crd_i - 1), (crd_j), (crd_k + 1), vox_array)
  a_8 = Get_check_voxel((crd_i + 1), (crd_j), (crd_k - 1), vox_array)
  p_8 =  Get_polar_2_axial(pi / 4, 4 * pi / 4)
  
  v_9 = Get_check_voxel((crd_i), (crd_j), (crd_k + 1), vox_array)
  a_9 = Get_check_voxel((crd_i), (crd_j), (crd_k - 1), vox_array)
  p_9 =  c(0, 0, 1)
  
  v_10 = Get_check_voxel((crd_i - 1), (crd_j - 1), (crd_k), vox_array)
  a_10 = Get_check_voxel((crd_i + 1), (crd_j + 1), (crd_k), vox_array)
  p_10 =  Get_polar_2_axial(pi / 2, 5 * pi / 4)
  
  v_11 = Get_check_voxel((crd_i), (crd_j - 1), (crd_k), vox_array)
  a_11 = Get_check_voxel((crd_i), (crd_j + 1), (crd_k), vox_array)
  p_11 =  c(0, 1, 0)
  
  v_12 = Get_check_voxel((crd_i + 1), (crd_j - 1), (crd_k), vox_array)
  a_12 = Get_check_voxel((crd_i - 1), (crd_j + 1), (crd_k), vox_array)
  p_12 =  Get_polar_2_axial(pi / 2, 7 * pi / 4)
  
  v_13 = Get_check_voxel((crd_i + 1), (crd_j), (crd_k), vox_array)
  a_13 = Get_check_voxel((crd_i - 1), (crd_j), (crd_k), vox_array)
  p_13 =  c(1, 0, 0)
  
  
  prior = matrix(NA, 13, 3)
  
  vx = list(v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8, v_9, v_10, v_11, v_12, v_13)
  ax = list(a_1, a_2, a_3, a_4, a_5, a_6, a_7, a_8, a_9, a_10, a_11, a_12, a_13)
  px = list(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10, p_11, p_12, p_13)
  
  wt = c()
  ## In the next loop select weight and prior from all 13 directions by
  #(a) choosing candidate fiber from vx and ax (if non-NA)
  #(b) then find the weight and select the prior
  for (l in 1:length(vx)) {
    if (is.na(vx[l]) & (is.na(ax[l]))) {
      prior[l, ] = c(0.5, 0.5, 0.7071068)
      wt <- c(wt, 0)
      
    } else if (is.na(vx[l])) {
      fib2 = Get_select_fiber(ax[l], unlist(px[l]))
      th2 = unlist(fib2)[1]
      ph2 = unlist(fib2)[2]
      prior[l, ] = Get_polar_2_axial(th2, ph2)
      wt <- c(wt, prior[l, ] %*% unlist(px[l]))
      
      
    } else if (is.na(ax[l])) {
      fib1 = Get_select_fiber(vx[l], unlist(px[l]))
      th1 = unlist(fib1)[1]
      ph1 = unlist(fib1)[2]
      prior[l, ] = Get_polar_2_axial(th1, ph1)
      wt <- c(wt, prior[l, ] %*% unlist(px[l]))
      
    } else {
      fib1 = Get_select_fiber(vx[l], unlist(px[l]))
      fib2 = Get_select_fiber(ax[l], unlist(px[l]))
      th1 = unlist(fib1)[1]
      ph1 = unlist(fib1)[2]
      th2 = unlist(fib2)[1]
      ph2 = unlist(fib2)[2]
      # x=Get_orient_summary(Get_orient_mat(c(th1,th2),c(ph1,ph2)))
      # Get_axial_2_polar(x[1],x[2],x[3])
      prior[l, ] = Get_orient_summary(Get_orient_mat(c(th1, th2), c(ph1, ph2)))
      wt <- c(wt, prior[l, ] %*% unlist(px[l]))
    }
    
  }# for loop ends
  
  wt[is.na(wt)] <- 0
  
  test = sample(seq(1:13), num.fib, replace = F, abs(wt))
  out = prior[test, ]
  cnt = rep(0, 13)
  cnt[test] = 1
  results = list("out" = out, "cnt" = cnt)
  return(results)
  
} #meu calculation function ends

