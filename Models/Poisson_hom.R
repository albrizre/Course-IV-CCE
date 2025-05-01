Poisson_hom <- nimbleCode({
  
  lambda0 ~ dgamma(0.1*lambda0_base,0.1)
  
  for (i in 1:N_events) {
    
    lambda_event[i] <- lambda0
    log_L[i] <- log(lambda_event[i])
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  for (i in (N_events+1):N) { 
    
    lambda_int[1:S,i-N_events] <- rep(lambda0,S)
    # Integral approximation 
    log_L[i] <- -sum(lambda_int[1:S,i-N_events])*area_cell*1  
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
})