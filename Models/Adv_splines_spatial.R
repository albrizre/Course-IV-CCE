Adv_splines_spatial <- nimbleCode({
  
  lambda0 ~ dgamma(0.1*lambda0_base,0.1)
  
  for (j in 1:N_splines){
    beta[j] ~ dnorm(0,sd=10)
  }
  
  epsilon[1] ~ dnorm(0,sd=sigma_epsilon)
  for (W in 2:53){
    epsilon[W] ~ dnorm(epsilon[W-1],sd=sigma_epsilon)
  }
  sigma_epsilon ~ dunif(0,10)
  
  for (i in 1:N_events) {
    
    lambda_event[i] <- lambda0*exp(beta[1]*X_splines[Id_cell[i],1]+beta[2]*X_splines[Id_cell[i],2]+
                                   beta[3]*X_splines[Id_cell[i],3]+beta[4]*X_splines[Id_cell[i],4]+
                                   beta[5]*X_splines[Id_cell[i],5]+beta[6]*X_splines[Id_cell[i],6]+
                                   beta[7]*X_splines[Id_cell[i],7]+beta[8]*X_splines[Id_cell[i],8]+
                                   beta[9]*X_splines[Id_cell[i],9])*
                               exp(epsilon[Week_events[i]])
    log_L[i] <- log(lambda_event[i])
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  for (i in (N_events+1):N) { 
    
    lambda_int[1:S,i-N_events] <- lambda0*exp(beta[1]*X_splines[1:S,1]+beta[2]*X_splines[1:S,2]+
                                              beta[3]*X_splines[1:S,3]+beta[4]*X_splines[1:S,4]+
                                              beta[5]*X_splines[1:S,5]+beta[6]*X_splines[1:S,6]+
                                              beta[7]*X_splines[1:S,7]+beta[8]*X_splines[1:S,8]+
                                              beta[9]*X_splines[1:S,9])*
                                          exp(epsilon[Week_time_points[i-N_events]])
    # Integral approximation 
    log_L[i] <- -sum(lambda_int[1:S,i-N_events])*area_cell*1  
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
})