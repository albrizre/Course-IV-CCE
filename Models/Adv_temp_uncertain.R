Adv_temp_uncertain <- nimbleCode({
  
  lambda0 ~ dgamma(0.1*lambda0_base,0.1)
  
  for (j in 1:4){
    beta[j] ~ dnorm(0,sd=10)
  }
  
  omega ~ dunif(0,1000)
  
  for (i in 1:N_events) {
    
    time_events[i] ~ dunif(time_events_start[i],time_events_end[i])
    
  }
  
  for (i in 1:N_events) {
    
    lambda_event[i] <- lambda0*exp(beta[1]*X1[Id_cell[i]]+beta[2]*X2[Id_cell[i]]+beta[3]*X3[Id_cell[i]]+beta[4]*X4[Id_cell[i]])*
      exp(-abs(time_events[i]-200)/omega)
    log_L[i] <- log(lambda_event[i])
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
  for (i in (N_events+1):N) { 
    
    lambda_int[1:S,i-N_events] <- lambda0*exp(beta[1]*X1[1:S]+beta[2]*X2[1:S]+beta[3]*X3[1:S]+beta[4]*X4[1:S])*
      exp(-abs(time_points[i-N_events]+0.5-200)/omega)
    # Integral approximation 
    log_L[i] <- -sum(lambda_int[1:S,i-N_events])*area_cell*1  
    z[i] <- -log_L[i]+C
    zeros[i] ~ dpois(z[i])
    
  }
  
})