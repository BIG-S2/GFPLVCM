beta_calculate <- function(s_time_point, u_time, len_k, B){
  phi=function(t,k){
    if(k==1)  return(1)  else  return(sqrt(2)*cos((k-1)*pi*t) );
  }

  time_number <- length(s_time_point)


  kk=matrix(rep(c( 1:len_k ), u_time ),  len_k, u_time );
  u=matrix(rep(seq(0,1,length=u_time),len_k ), len_k, u_time, byrow=T);

  phi_u=matrix(mapply(function(t,k)phi(t,k),u,kk), len_k, u_time);


  kk_s=matrix(rep(c( 1:len_k ), time_number ),  len_k, time_number );
  s=matrix(rep(s_time_point ,len_k ), len_k, time_number, byrow=T);

  phi_s=matrix(mapply(function(t,k)phi(t,k),s,kk_s), len_k, time_number);

  beta_temp_value=0

  for(k in 1: len_k){
    beta_temp_value=beta_temp_value +  as.matrix( phi_s[k, ] )%*%phi_u[k, ]*(-1)^(k+1)/(k^2)   ## s*u
  }

  beta_value=beta_temp_value*B

  return(beta_value)

}
