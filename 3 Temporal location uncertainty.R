library(nimble)
library(sp)
library(sf)
library(coda)
library(RColorBrewer)
library(ggplot2)
library(ggspatial)

# setwd(...) in a folder containing all the folders of the repository
source("Models/Adv_temp_uncertain.R")

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

# Inserting NAs in some time_events and adding temporal interval
set.seed(123)
unc_events=sort(sample(1:length(time_events),size=0.4*length(time_events)))
time_events_start=time_events
time_events_end=time_events
for (j in unc_events){
  time_events_start[j]=time_events[j]-rexp(1,rate=1/3)
  time_events_end[j]=time_events[j]+rexp(1,rate=1/3)
  if (time_events_start[j]<0){time_events_start[j]=0}
  if (time_events_end[j]>365){time_events_end[j]=364.5}
}
head(cbind(time_events,time_events_start,time_events_end))

# Constants for model setting
constants <- list(
  
  N = length(time_events)+length(time_points),
  N_events = length(time_events),
  N_integral = length(time_points),
  time_events_start = time_events_start,
  time_events_end = time_events_end,
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
  X4 = covariates$renta_hogar
  
)

# Fit model
set.seed(12345)
data <- list(zeros = rep(0,constants$N))
inits <- function() list(lambda0 = as.numeric(constants$lambda0_base),
                         beta = rep(0,4),
                         omega = 1,
                         time_events = time_events)

monitors_pars <- c("lambda0","beta","lambda_int","omega","time_events")
Sys.time()
mcmc.output <- nimbleMCMC(Adv_temp_uncertain, data = data, inits = inits, constants = constants,
                          monitors = monitors_pars, 
                          niter = 20000, nburnin = 10000, nchains = 1, thin = 10,
                          summary = TRUE, WAIC = TRUE)
Sys.time()

# Check results
head(mcmc.output$summary$all.chains)
tail(mcmc.output$summary$all.chains)
View(mcmc.output$summary)
mcmc.output$summary[grep("time_events",rownames(mcmc.output$summary))[1:6],]
