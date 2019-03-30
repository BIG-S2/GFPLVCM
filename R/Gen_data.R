Gen_data<-function(n, pois_par, u_time, gamma, len_k, B){
  ## n: number of observations
  ## pois_par: number of observation times of x and y follow a possion distribution with this parameter
  ## u_time: number of observation times for the dimension u
  ## B: signal of beta
  phi=function(t,k){
    if(k==1)  return(1)  else  return(sqrt(2)*cos((k-1)*pi*t) );
  }


  y_time_number= stats::rpois(n, pois_par)
  y_time_number[y_time_number==0]=pois_par
  z_time_number= stats::rpois(n, pois_par)
  z_time_number[z_time_number==0]=pois_par
  x_time_number= z_time_number    ### t dimension of x shares the same time point with z, but different from y

  #y_time_number=rep(pois_par, n) # same number
  #z_time_number=rep(pois_par, n) #
  #x_time_number= z_time_number #

  y_ID=NULL
  z_ID=NULL
  y_time_point=NULL
  z_time_point=NULL
  y_value=NULL
  z_value=NULL
  x_value=NULL

  #y_temp_time_point = sort( runif( pois_par, 0 , 1) )  # same time point
  #z_temp_time_point = sort( runif( pois_par, 0 , 1) )  #

  for(i in 1:n){
    y_ID=c(y_ID, rep(i, y_time_number[i]) )
    z_ID=c(z_ID, rep(i, z_time_number[i]) )

    y_temp_time_point = sort( stats::runif( y_time_number[i], 0 , 1) )
    z_temp_time_point = sort( stats::runif( z_time_number[i], 0 , 1) )

    y_time_point=c(y_time_point, y_temp_time_point)
    z_time_point=c(z_time_point, z_temp_time_point)

    range_time_point <- c( y_temp_time_point , z_temp_time_point)
    len_y <- y_time_number[i]
    len_z <- z_time_number[i]

    ###########################################
    ### generate x and z ######################
    ###########################################

    gen_xz=Gen_x_z( len_y + len_z, range_time_point, u_time, len_k )

    z_value=c(z_value, gen_xz$z_temp_value[-c(1:len_y)])

    x_value=rbind(x_value, gen_xz$x_temp_value[-c(1:len_y), ])


    ###################################
    ### generate y ####################
    ###################################
    z_temp_value_for_y <- gen_xz$z_temp_value[c(1:len_y)]
    x_temp_value_for_y <- matrix( gen_xz$x_temp_value[c(1:len_y), ], len_y, u_time)
    beta_temp_value_for_y <- matrix( gen_xz$beta_temp_value[c(1:len_y), ], len_y, u_time )*B

    sigma_epsilon <- 2^( -abs( outer(y_temp_time_point, y_temp_time_point, "-") ) )

    epsilon_temp <- MASS::mvrnorm( n=1,rep(0,y_time_number[i]), sigma_epsilon )


    y_temp_value= z_temp_value_for_y*gamma + rowSums( x_temp_value_for_y*beta_temp_value_for_y )/u_time + epsilon_temp

    y_value=c(y_value, y_temp_value)

  }

  y=list("y_ID"=y_ID, "y_time_point"=y_time_point, "y_value"=y_value)
  z=list("z_ID"=z_ID, "z_time_point"=z_time_point, "z_value"=as.matrix(z_value) )
  x=list("x_ID"=z_ID, "x_time_point"=z_time_point, "x_value"=x_value)


  list(y=y, x=x, z=z)


}
