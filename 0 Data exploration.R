library(ggplot2)
library(ggspatial)
library(RColorBrewer)

# setwd(...) in a folder containing all the folders of the repository

load("Data/covariates.rda")
load("Data/grid_val.rda")
load("Data/events_study_area.rda")

plot(covariates)
plot(grid_val)
plot(events_study_area)
plot(grid_val)
plot(events_study_area,add=T,col="red")

# ggplot

grid_val$poblacion_65_mas=covariates$poblacion_65_mas
grid_val$poblacion_15_29=covariates$poblacion_15_29
grid_val$poblacion_extranjera=covariates$poblacion_extranjera
grid_val$renta_hogar=covariates$renta_hogar
aux_gg=as(grid_val,"sf")

ggplot() +
  annotation_spatial(aux_gg) +
  layer_spatial(aux_gg)+
  theme_bw()+
  theme(text = element_text(size = 24),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        title = element_text(colour = "black",face = "bold"),
        legend.position = "right")

ggplot() +
  annotation_spatial(aux_gg) +
  layer_spatial(aux_gg,aes(fill=renta_hogar))+
  theme_bw()+
  theme(text = element_text(size = 24),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        title = element_text(colour = "black",face = "bold"),
        legend.position = "right")

proj4string(events_study_area)="+proj=utm +zone=30 ellps=WGS84"
aux_events_gg=as(events_study_area,"sf")
time_events=events_study_area$Time
ggplot() +
  layer_spatial(aux_events_gg,aes(col=time_events),size=3)+
  scale_color_gradientn(colors=brewer.pal(6,"BuPu")[2:5],name=expression(italic(t)[italic(i)]),
                        limits=c(0,365),breaks=c(0,100,200,300,365))+
  theme_bw()+
  theme(text = element_text(size = 24),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        title = element_text(colour = "black",face = "bold"),
        legend.position = "right")
