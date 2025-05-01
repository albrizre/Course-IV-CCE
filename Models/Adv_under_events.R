Adv_under_events <- nimbleCode({
  
  lambda0 ~ dgamma(0.1*lambda0_base,0.1)
  
  alpha ~ dnorm(0,sd=10)
  
  epsilon[1] ~ dnorm(0,sd=sigma_epsilon)
  for (W in 2:53){
    epsilon[W] ~ dnorm(epsilon[W-1],sd=sigma_epsilon)
    pi[W] <- exp(alpha+epsilon[W])/(1+exp(alpha+epsilon[W]))
  }
  sigma_epsilon ~ dunif(0,10)
  
  for (i in 1:N_events) {
    
    lambda_event[i] <- lambda0*pi[Week_events[i]]
    log_L[i] <- log(lambda_event[i])
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  for (i in (N_events+1):N) { 
    
    lambda_int[1:S,i-N_events] <- rep(lambda0*pi[Week_time_points[i-N_events]],S)
    # Integral approximation 
    log_L[i] <- -sum(lambda_int[1:S,i-N_events])*area_cell*1  
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
})