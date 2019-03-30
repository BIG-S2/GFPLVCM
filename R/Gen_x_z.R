Gen_x_z<-function(time_number, time_point, u_time, len_k){

  ## generate x and z for one subject
  phi=function(t,k){
    if(k==1)  return(1)  else  return(sqrt(2)*cos((k-1)*pi*t) );
  }


  ###################################################################
  ### generate the correlation relationship within the subject ######
  ###################################################################
 # sigma_z_temp=diag( time_number )

 #  for(j in 1: time_number ){
  #  for(k in j: time_number ){
   #   sigma_z_temp[j, k] = exp(- abs ( time_point[j] - time_point[k] ) )
    #  sigma_z_temp[k, j]=sigma_z_temp[j, k]
  #  }
 #  }

  sigma_z_temp = exp( -abs( outer(time_point, time_point, "-") ) )

  z_temp_value = MASS::mvrnorm( n=1,rep(0,time_number), sigma_z_temp )


  ######################################
  ### generate x for a subject #########
  ######################################

  kk=matrix(rep(c( 1:len_k ), u_time ),  len_k, u_time );
  u=matrix(rep(seq(0,1,length=u_time),len_k ), len_k, u_time, byrow=T);

  phi_u=matrix(mapply(function(t,k)phi(t,k),u,kk), len_k, u_time);


  M = MASS::mvrnorm( n=len_k,rep(0,time_number), sigma_z_temp )  ## len_k*length of S

  x_temp_value=0

  for(k in 1: len_k){
    x_temp_value=x_temp_value +  as.matrix( M[k, ] )%*%phi_u[k, ]*(-1)^(k+1)/k
  }


  #################################
  ### generate beta ###############
  #################################
  kk_s=matrix(rep(c( 1:len_k ), time_number ),  len_k, time_number );
  s=matrix(rep(time_point ,len_k ), len_k, time_number, byrow=T);

  phi_s=matrix(mapply(function(t,k)phi(t,k),s,kk_s), len_k, time_number);

  beta_temp_value=0

  for(k in 1: len_k){
    beta_temp_value=beta_temp_value +  as.matrix( phi_s[k, ] )%*%phi_u[k, ]*(-1)^(k+1)/(k^2)   ## s*u
  }


  list(z_temp_value=z_temp_value, x_temp_value=x_temp_value, beta_temp_value=beta_temp_value)


}
