
########Define the Matching function to calculate Watson meu given any X,Y,Z coordinate of the voxel and the voxel cube

Get_check_voxel = function(crd_i, crd_j, crd_k, vox_array) {
  vx_mx = 3
  vx_min = 1
  mx = max(crd_i, crd_j, crd_k)
  mn = min(crd_i, crd_j, crd_k)
  if (mx > vx_mx | mn < vx_min) {
    return(NA)
  } else {
    return(vox_array[crd_i, crd_j, crd_k, ])
  }
}


