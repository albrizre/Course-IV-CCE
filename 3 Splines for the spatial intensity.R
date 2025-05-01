library(nimble)
library(sp)
library(sf)
library(coda)
library(RColorBrewer)
library(ggplot2)
library(ggspatial)
library(npreg)

# setwd(...) in a folder containing all the folders of the repository
source("Models/Adv_splines_spatial.R")

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

# Spline basis over the spatial window
centroids_grid=sf::st_centroid(grid_val_sf)
aux=as.numeric(unlist(centroids_grid$geometry))
x=aux[seq(1,length(aux),2)]
y=aux[seq(2,length(aux),2)]
input=cbind(x,y)
k=4
knots=expand.grid(quantile(x,1:(k-1)/k),quantile(y,1:(k-1)/k))
X_splines=basis.tps(input, knots, m = 2, rk = TRUE, intercept = FALSE, ridge = FALSE)

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
  Week_time_points = floor(time_points/7)+1,
  X_splines = scale(X_splines[,3:ncol(X_splines)])
  
)
constants$N_splines=ncol(constants$X_splines)

# Fit model
set.seed(12345)
data <- list(zeros = rep(0,constants$N))
inits <- function() list(lambda0 = as.numeric(constants$lambda0_base),
                         beta = rep(0,constants$N_splines),
                         epsilon = rep(0,53),
                         sigma_epsilon = 5)

monitors_pars <- c("lambda0","beta","lambda_int")
Sys.time()
mcmc.output <- nimbleMCMC(Adv_splines_spatial, data = data, inits = inits, constants = constants,
                          monitors = monitors_pars, 
                          niter = 20000, nburnin = 10000, nchains = 1, thin = 10,
                          summary = TRUE, WAIC = TRUE)
Sys.time()

# Check results
head(mcmc.output$summary$all.chains)
tail(mcmc.output$summary$all.chains)

# Plot splines
aux=SpatialPolygonsDataFrame(grid_val,data.frame(constants$X_splines),match.ID = F)
proj4string(aux)="+proj=utm +zone=30 ellps=WGS84"
aux_gg=as(aux,"sf")

ggplot() +
  annotation_spatial(aux_gg) +
  layer_spatial(aux_gg,aes(fill=knot.2))+
  scale_fill_gradientn(colors=brewer.pal(8,"YlOrRd"),name=expression(X[italic(s)]))+
  theme_bw()+
  ggtitle(paste0("Spline function ",2))+
  theme(text = element_text(size = 24),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        title = element_text(colour = "black",face = "bold"),
        legend.position = "right")

# Intensity map
select=grep("lambda_int",colnames(mcmc.output$samples))
day=200
estimates=as.numeric(mcmc.output$summary[select[((day-1)*length(grid_val)+1):(day*length(grid_val))],1])
aux=SpatialPolygonsDataFrame(grid_val,data.frame(estimates),match.ID = F)
proj4string(aux)="+proj=utm +zone=30 ellps=WGS84"
aux_gg=as(aux,"sf")

ggplot() +
  annotation_spatial(aux_gg) +
  layer_spatial(aux_gg,aes(fill=estimates))+
  scale_fill_gradientn(colors=brewer.pal(8,"YlOrRd"),name=expression(hat(lambda)[italic(st)]))+
  theme_bw()+
  ggtitle(paste0("t = ",day))+
  theme(text = element_text(size = 24),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        title = element_text(colour = "black",face = "bold"),
        legend.position = "right")
