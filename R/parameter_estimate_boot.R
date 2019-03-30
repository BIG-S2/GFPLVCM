parameter_estimate_boot <-function(y, x, z, lambda_s, lambda_u, bd, order_1, order_2, breaks,boot_R, gamma_est_or, b_est_or){
  ### lambda_s: tuning parameter on the direction s of beta(s,u)
  ### lambda_u: tuning parameter on the direction u of beta(s,u)
  ### bd: bandwidth of the kernel function
  ### signal: for RMSE of beta
  ### xi: weight

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

  ##############

  p <- dim(as.matrix(z$z_value) )[2]

  unique_ID <- unique(y$y_ID)
    ### total observations of y
  n <- length( unique_ID )

  time_range <- c( y$y_time_point, x$x_time_point )    ### time range of x, z and y

  u_time <- dim( x$x_value )[2]

  ###### create tilde x ############
  x_ID <- x$x_ID
  nr_x <- length(x_ID)
  x_time_point <- x$x_time_point

  tilde_x <- list(ID=x_ID, time_point=x_time_point, value = NULL)

  theta_value <- fda::bsplineS( seq(0,1,length=u_time), breaks=breaks, norder=order_2, nderiv=0, returnMatrix=FALSE )

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


  ##################################
  ## bootstrap the distribution ####
  ##################################
  boot_gamma <- NULL
  boot_b <- NULL

  for(k in 1: boot_R){

    #print(k)

    #xi <- rexp(n)

	xi <- as.numeric( stats::rmultinom(1, n, rep(1/n,n)) )

    sum_y_z <- 0

    sum_z_z <- 0

    for( i in 1:n){

      #####################################
      #### estimates ######################
      #####################################


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

      kernel_value <-  apply( gap_time_point, 1:2, function(t) k_h(t, bd) )*xi[i]  ## K_h (t -s ) within the same subject

      #Kernel_with_y <- kernel_value * matrix( rep( temp_y_value, len_index_z ), len_index_y , len_index_z ) ### K*y

      sum_y_z <- sum_y_z + temp_y_value %*% kernel_value%*%temp_z_value

      temp_kernel_s <- matrix( rep( colSums(kernel_value),  dim(temp_z_value)[2] ), len_index_z,  dim(temp_z_value)[2])

      sum_z_z <- sum_z_z + t( temp_kernel_s*temp_z_value )%*% temp_z_value
    }

    design_matrix <- solve( sum_z_z + sum(xi)*Pen_all )

    parameter_est <- design_matrix%*%t( sum_y_z )

    gamma_est <- parameter_est[ 1: p ]

    b_est <- parameter_est[ -c(1:p )]




    ##############################
    ### calculate statistics #####
    ##############################
    #gamma_stat <- t( gamma_est - final_beta$gamma_est )%*% solve( gamma_variance  ) %*% ( gamma_est - final_beta$gamma_est )
    #b_stat <- t(b_est - final_beta$b_est)%*%ginv(b_variance)%*%(b_est - final_beta$b_est)

	gamma_stat <- ( gamma_est - gamma_est_or )^2
	b_stat <-  sum( (b_est - b_est_or)^2 )

    boot_gamma <- rbind(boot_gamma, gamma_stat)
    boot_b <- c(boot_b, b_stat)
  }

  gamma_stat_or <- matrix(rep( gamma_est_or^2, boot_R), boot_R, p, byrow=T)
  b_stat_or <- sum(b_est_or^2)

  gamma_pvalue <- ( 1 + colSums(boot_gamma > gamma_stat_or )) /(1 + boot_R)

  b_pvalue <- ( 1 + sum(boot_b > b_stat_or ) ) /(1 + boot_R)

  gamma_judge <- gamma_pvalue < 0.05
  b_judge <- b_pvalue < 0.05


  list(boot_gamma = boot_gamma, boot_b =boot_b, gamma_judge=gamma_judge, b_judge=b_judge,
       gamma_pvalue = gamma_pvalue, b_pvalue=b_pvalue)


}
