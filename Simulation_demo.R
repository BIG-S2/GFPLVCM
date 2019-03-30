library(GFPLVCM)  ### dependencies: MASS, fda, stats

n=100              ## sample size
pois_par=10        ## parameter for the Poission distribution to generate the observation times

order_1=4; order_2=4   ## order of the B-splines in the s and u directions
breaks=c(0,.2, .4, .6, .8, 1)   ## knots of the B-splines

###########################################################
u_time=200;      ## observation times for the functional variable in the u direction
len_k=50 ;       ## number of basis functions to generate the functional parameter
gamma_real=0.3;  ## true value of gamma
B = 4          ## signal of the functional parameter
pre_n <- 200   #### sample size of prediction
boot_R <- 1000  ### times for bootstrap
grid <- u_time

lambda_range <- seq(n^(-2), n^(-0.4), length=10)  ## range of the tuning parameter lambda

h_range <- seq(n^(-1), n^(-0.2), length=10)       ## range of the tuning parameter bandwidth

 
data <- Gen_data(n, pois_par, u_time, gamma_real, len_k, B)   ### generate the data

para <- tuning_parameter_selection( lambda_range, h_range, data$y, data$x, data$z, order_1, order_2, breaks, pois_par, grid)

final_beta <- parameter_estimate_prediction(data$y, data$x, data$z, para$final_lambda, para$final_lambda, para$final_bd, order_1, order_2, breaks,  pois_par,len_k, gamma_real, pre_n, B)

boot <- parameter_estimate_boot(data$y, data$x, data$z, para$final_lambda, para$final_lambda, para$final_bd, order_1, order_2, breaks, boot_R, final_beta$gamma_est, final_beta$b_est)

## do estimation, prediction and testing
