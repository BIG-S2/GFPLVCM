parameter_estimate<-function(y, x, z, lambda_s, lambda_u, bd, order_1, order_2, breaks,  pois_par, grid){
  ### lambda_s: tuning parameter on the direction s of beta(s,u)
  ### lambda_u: tuning parameter on the direction u of beta(s,u)
  ### bd: bandwidth of the kernel function

  ############################
  eta <- fda::create.bspline.basis( norder=order_1, breaks=breaks )  ## order_1 + 4 basis
  theta <- fda::create.bspline.basis( norder=order_2, breaks=breaks )


  #########################################
  ### create the penalty matrix ###########
  #########################################
  R <- fda::bsplinepen( eta )
  S <- fda::bsplinepen( theta )

  J_eta <- fda::bsplinepen( eta, Lfdobj= 0 )
  J_theta <- fda::bsplinepen( theta, Lfdobj= 0 )

  Pen_matrix_s <- kronecker( J_eta, R  )   ### symmetric
  Pen_matrix_u <- kronecker( S, J_theta  )


  k_h <- function(t, h){
  return( 0.75*max( 0, (1 - (t/h)^2   ) )/h )
  }

  p <- dim(as.matrix(z$z_value) )[2]

  unique_ID <- unique(y$y_ID)
  n <- length( unique_ID )

  time_range <- c( y$y_time_point, x$x_time_point )    ### time range of x, z and y

  u_time <- dim( x$x_value )[2]

  if(length(grid)==1){
    grid <- seq(0,1,length=u_time)
  }

  ###### create tilde x ############
  x_ID <- x$x_ID
  nr_x <- length(x_ID)
  x_time_point <- x$x_time_point

  tilde_x <- list(ID=x_ID, time_point=x_time_point, value = NULL)

  theta_value <- fda::bsplineS( grid, breaks=breaks, norder=order_2, nderiv=0, returnMatrix=FALSE )

  for(i in 1:nr_x){
     temp_time_point <- x_time_point[ i ]
     eta_value <- fda::bsplineS( temp_time_point, breaks=breaks, norder=order_1, nderiv=0, returnMatrix=FALSE ) ## eta_k change with the subject

     tilde_x_value <-  t(eta_value) %*% x$x_value[i, ] %*% theta_value/u_time   ### K_1 *K_2 matrix, col: 1, ... K_1; row: 1, ..., K_2
     tilde_x$value <- rbind(tilde_x$value, as.vector(tilde_x_value ) )

  }

  tilde_z <- list(ID=x_ID, time_point = x_time_point, value = cbind(z$z_value, tilde_x$value) )





  ###################################
  ## create the penalty matrix ######
  ###################################
  K_1 <- order_1 + length(breaks) - 2
  K_2 <- order_2 + length(breaks) - 2

  Pen_b <- lambda_s*Pen_matrix_s + lambda_u*Pen_matrix_u
  Pen_all <- diag(rep(0, K_1 * K_2+p) )
  Pen_all[ (p+1): (K_1 * K_2+p),   (p+1): (K_1 * K_2+p)]=Pen_b


  ###################################
  ### calculate the summations ######
  ###################################

  sum_y_z <- 0

  sum_z_z <- 0

  for( i in 1:n){

    index_z <- NULL
    index_y <- NULL
    index_z <- which( tilde_z$ID == unique_ID[i] )
    index_y <- which( y$y_ID == unique_ID[i] )

    len_index_z <- length(index_z)
    len_lndex_y <- length(index_y)

    time_z <- tilde_z$time_point[ index_z ]
    time_y <- y$y_time_point[ index_y ]

    temp_z_value <- matrix( tilde_z$value[ index_z, ], len_index_z, p+K_1*K_2 )
    temp_y_value <- y$y_value[ index_y ]

    gap_time_point <- outer(time_y, time_z, "-")

    kernel_value <-  apply( gap_time_point, 1:2, function(t) k_h(t, bd) )  ## K_h (t -s ) within the same subject

    #Kernel_with_y <- kernel_value * matrix( rep( temp_y_value, len_index_z ), len_index_y , len_index_z ) ### K*y

    sum_y_z <- sum_y_z + temp_y_value %*% kernel_value%*%temp_z_value

    temp_kernel_s <- matrix( rep( colSums(kernel_value),  dim(temp_z_value)[2] ), len_index_z,  dim(temp_z_value)[2])

    sum_z_z <- sum_z_z + t( temp_kernel_s*temp_z_value )%*% temp_z_value
  }


  parameter_est <- solve( sum_z_z + n*Pen_all )%*%t( sum_y_z )

  gamma_est <- parameter_est[ 1: p ]

  b_est <- parameter_est[ -c(1:p )]



  ######################################
  ### calculate mse of beta ############
  ######################################
  if(length(pois_par)==1 ){
    s_time_point <- stats::runif(2*pois_par, 0, 1)
  }else{
    s_time_point <- pois_par
  }


 # beta_true <- beta_calculate(s_time_point, u_time, len_k, B)

  B_est <- matrix(b_est, K_1, K_2)

  eta_value_for_beta <- fda::bsplineS( s_time_point, breaks=breaks, norder=order_1, nderiv=0, returnMatrix=FALSE ) ## eta_k change with the subject
  theta_value_for_beta <- fda::bsplineS( grid, breaks=breaks, norder=order_2, nderiv=0, returnMatrix=FALSE )

  beta_est <- eta_value_for_beta%*%B_est%*%t(theta_value_for_beta)

 # mse_beta <- mean( ( beta_est - beta_true )^2 )

 # mse_gamma <- sum((gamma_est - gamma_real)^2)
  gamma_stat <- sum(gamma_est^2)
  b_stat <- sum(b_est^2)


  list(gamma_est = gamma_est, b_est = b_est,  beta_est=beta_est, gamma_stat=gamma_stat, b_stat=b_stat)


}
