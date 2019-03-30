tuning_parameter_selection <- function( lambda_range, h_range, y, x, z, order_1, order_2, breaks, pois_par, grid){

  K_1 <- order_1 + length(breaks) - 2
  K_2 <- order_2 + length(breaks) - 2
  ###############################
  ##### split sample ############
  ###############################

  ID_y <- unique( y$y_ID )
  n  <- length(ID_y)

  sample_1_size <- round(n/2 ,0)

  sample_1_y_id <- NULL
  sample_1_z_id <- NULL

  for(i in 1: sample_1_size){

    sample_1_y_id <- c(sample_1_y_id, which( y$y_ID ==  ID_y[i] )  )
    sample_1_z_id <- c(sample_1_z_id, which( z$z_ID ==  ID_y[i] )  )

  }
  sample_1_x_id <- sample_1_z_id

  sample_1_y <- list( y_ID = y$y_ID[sample_1_y_id], y_time_point = y$y_time_point[sample_1_y_id], y_value = y$y_value[sample_1_y_id]  )
  sample_2_y <- list( y_ID = y$y_ID[-sample_1_y_id], y_time_point = y$y_time_point[-sample_1_y_id], y_value = y$y_value[-sample_1_y_id]  )

  sample_1_x <- list( x_ID = x$x_ID[sample_1_x_id], x_time_point = x$x_time_point[sample_1_x_id], x_value = x$x_value[sample_1_x_id, ]  )
  sample_2_x <- list( x_ID = x$x_ID[-sample_1_x_id], x_time_point = x$x_time_point[-sample_1_x_id], x_value = x$x_value[-sample_1_x_id, ]  )

  sample_1_z <- list( z_ID = z$z_ID[sample_1_z_id], z_time_point = z$z_time_point[sample_1_z_id], z_value = z$z_value[sample_1_z_id, ]  )
  sample_2_z <- list( z_ID = z$z_ID[-sample_1_z_id], z_time_point = z$z_time_point[-sample_1_z_id], z_value = z$z_value[-sample_1_z_id, ]  )

  ####################################
  #### calculate criterion ###########
  ####################################

  length_lambda <- length( lambda_range )
  length_h <- length( h_range )
  criterion_beta <- matrix( rep(100, length_lambda*length_h ), length_lambda, length_h)


  for(i in 1:length_lambda){
    for(j in 1:length_h){

      lambda_temp <- lambda_range[i]
      bd_temp <- h_range[j]
      sample_1_est <- parameter_estimate(sample_1_y, sample_1_x, sample_1_z, lambda_temp, lambda_temp, bd_temp,
                                         order_1, order_2, breaks,  pois_par, grid)

      sample_1_b_est <- sample_1_est$b_est
      sample_1_gamma_est <- sample_1_est$gamma_est


      sample_2_est <- parameter_estimate(sample_2_y, sample_2_x, sample_2_z, lambda_temp, lambda_temp, bd_temp,
                                         order_1, order_2, breaks,  pois_par, grid)

      sample_2_b_est <- sample_2_est$b_est
      sample_2_gamma_est <- sample_2_est$gamma_est

      #####
      bias_square <-   lambda_temp + K_1*K_2*bd_temp^4

      criterion_temp <- mean( ( sample_1_b_est -  sample_2_b_est)^2  )/4 + sum( (sample_1_gamma_est- sample_1_gamma_est )^2)/4
      criterion_beta[i, j] <- criterion_temp + bias_square

    }
  }

  position_beta <- which( criterion_beta == min( criterion_beta  ), arr.ind=T  )

  lambda_beta <- lambda_range[position_beta[1]]
  bd_beta <- h_range[position_beta[2]]


  list( final_lambda=lambda_beta,  final_bd= bd_beta )

}
