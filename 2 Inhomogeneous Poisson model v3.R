library(nimble)
library(sp)
library(sf)
library(coda)
library(RColorBrewer)
library(ggplot2)
library(ggspatial)

# setwd(...) in a folder containing all the folders of the repository
source("Models/Poisson_inhom3.R")

# Data loading
load("Data/covariates.rda")
load("Data/grid_val.rda")
load("Data/events_study_area.rda")
grid_val_sf=as(grid_val,"sf")

# Auxiliary variables
time_events=events_study_area$Time
area_cell=max(st_area(grid_val_sf))
Tmax=365
time_points=seq(0,Tmax,1)

# Constants for model setting
constants <- list(
  
  N = length(time_events)+length(time_points),
  N_events = length(time_events),
  N_integral = length(time_points),
  time_events = time_events,
  time_points = time_points,
  lambda0_base = length(time_events)/(sum(st_area(grid_val_sf))*Tmax),
  time_all = c(time_events,time_points),
  C = 100000000,
  S = length(grid_val),
  Id_cell = match(events_study_area$IdS,grid_val$IdS),
  Tmax = Tmax,
  area_cell = area_cell,
  X1 = covariates$poblacion_65_mas,
  X2 = covariates$poblacion_15_29,
  X3 = covariates$poblacion_extranjera,
  X4 = covariates$renta_hogar,
  Week_events = floor(time_events/7)+1,
  Week_time_points = floor(time_points/7)+1
  
)

# Fit model
set.seed(12345)
data <- list(zeros = rep(0,constants$N))
inits <- function() list(lambda0 = as.numeric(constants$lambda0_base),
                         beta = rep(0,4),
                         epsilon = rep(0,53),
                         sigma_epsilon = 5)

monitors_pars <- c("lambda0","beta","sigma_epsilon","epsilon")
Sys.time()
mcmc.output <- nimbleMCMC(Poisson_inhom3, data = data, inits = inits, constants = constants,
                          monitors = monitors_pars, 
                          niter = 20000, nburnin = 10000, nchains = 2, thin = 10,
                          summary = TRUE, WAIC = TRUE)
Sys.time()
# save(mcmc.output,file="Outputs/Poisson_inhom3.rda")

# Check results
head(mcmc.output$summary$all.chains)
tail(mcmc.output$summary$all.chains)

# Temporal effect

df_plot=c()
epsilons=mcmc.output$samples$chain1[,grep("epsilon",colnames(mcmc.output$samples$chain1))][,-1] # remove sigma_epsilon column
for (W in seq(1,53,1)){
  temp_effect=exp(epsilons[,W])
  df_plot=rbind(df_plot,data.frame(W,
                                   Lo=quantile(temp_effect,0.025),
                                   Mean=mean(temp_effect),
                                   Up=quantile(temp_effect,0.975)))
}
ggplot(df_plot,aes(x=W,y=Mean))+
  geom_ribbon(aes(ymin=Lo,ymax=Up),fill="gray60")+
  geom_line(linewidth=2)+
  theme_bw()+
  ylab("Week effect on the intensity")+
  geom_hline(yintercept = 1,lty="dashed")+
  theme(text = element_text(size = 24),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        title = element_text(colour = "black",face = "bold"),
        legend.position = "right")
