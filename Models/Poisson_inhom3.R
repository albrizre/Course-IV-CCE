Poisson_inhom3 <- nimbleCode({
  
  lambda0 ~ dgamma(0.1*lambda0_base,0.1)
  
  for (j in 1:4){
    beta[j] ~ dnorm(0,sd=10)
  }
  
  epsilon[1] ~ dnorm(0,sd=sigma_epsilon)
  for (W in 2:53){
    epsilon[W] ~ dnorm(epsilon[W-1],sd=sigma_epsilon)
  }
  sigma_epsilon ~ dunif(0,10)
  
  for (i in 1:N_events) {
    
    lambda_event[i] <- lambda0*exp(beta[1]*X1[Id_cell[i]]+beta[2]*X2[Id_cell[i]]+beta[3]*X3[Id_cell[i]]+beta[4]*X4[Id_cell[i]])*
                               exp(epsilon[Week_events[i]])
    log_L[i] <- log(lambda_event[i])
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  for (i in (N_events+1):N) { 
    
    lambda_int[1:S,i-N_events] <- lambda0*exp(beta[1]*X1[1:S]+beta[2]*X2[1:S]+beta[3]*X3[1:S]+beta[4]*X4[1:S])*
                                          exp(epsilon[Week_time_points[i-N_events]])
    # Integral approximation 
    log_L[i] <- -sum(lambda_int[1:S,i-N_events])*area_cell*1  
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
})