library(nimble)
library(sp)
library(sf)
library(coda)
library(RColorBrewer)
library(spatstat)
library(sp)
library(ggplot2)
library(ggspatial)
library(RColorBrewer)

# setwd(...) in a folder containing all the folders of the repository

# Data loading
load("Data/covariates.rda")
load("Data/grid_val.rda")
load("Data/events_study_area.rda")
grid_val_sf=as(grid_val,"sf")

# Residual analysis (spatial) for a given model

# Model loading
load("Outputs/Poisson_hom.rda")

# First we obtain a non-parametric estimate of the intensity
pattern=ppp(x = events_study_area@coords[,1],y = events_study_area@coords[,2],
            window = as.owin(sf::st_union(grid_val_sf)))
bw_select=spatstat.explore::bw.scott(pattern,isotropic = T)
funcion_density=densityfun(pattern,sigma = bw_select)
centroids=sf::st_coordinates(sf::st_centroid(grid_val_sf,byid = T))
density_values=funcion_density(centroids[,1],centroids[,2])

# Second, we compute the smoothed values of the fitted intensity
X=ppp(x=centroids[,1],y=centroids[,2],
      window = as.owin(sf::st_union(grid_val_sf)))
lambdas=as.data.frame(mcmc.output$summary$chain1[grep("lambda_int",
                      rownames(mcmc.output$summary$chain1)),])
aux=rep(1:647,366)
lambdas$Cell=aux
lambdas_avg_cell=sqldf::sqldf("SELECT Cell, SUM(Mean) as Model from lambdas
                               GROUP BY Cell")
marks(X)=lambdas_avg_cell
smooth_X=Smooth(X, bw_select,at = "points")[,2]

# All together and plots. Dif corresponds to the residual
aux=SpatialPolygonsDataFrame(grid_val,data.frame(Model=lambdas_avg_cell$Model,
                                                 Model_smoothed=smooth_X,
                                                 KDE=density_values,
                                                 Dif=density_values-smooth_X),match.ID = F)
proj4string(aux)="+proj=utm +zone=30 ellps=WGS84"
aux_gg=as(aux,"sf")

ggplot() +
  annotation_spatial(aux_gg) +
  layer_spatial(aux_gg,aes(fill=Model_smoothed),col="transparent")+ #,col="transparent"
  scale_fill_gradientn(colors=brewer.pal(8,"YlOrRd"),name=expression(hat(lambda)[italic(c)]),limits=c(0,1.5*10^(-4)))+
  theme_bw()+
  ggtitle("Average estimated intensity")+
  labs(subtitle = "Model-based")+
  theme(text = element_text(size = 24),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        title = element_text(colour = "black",face = "bold"),
        legend.position = "right")

ggplot() +
  annotation_spatial(aux_gg) +
  layer_spatial(aux_gg,aes(fill=KDE),col="transparent")+ #,col="transparent"
  scale_fill_gradientn(colors=brewer.pal(8,"YlOrRd"),name=expression(hat(lambda)[italic(c)]),limits=c(0,1.5*10^(-4)))+
  theme_bw()+
  ggtitle("Average estimated intensity")+
  labs(subtitle = "KDE")+
  theme(text = element_text(size = 24),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        title = element_text(colour = "black",face = "bold"),
        legend.position = "right")

ggplot() +
  annotation_spatial(aux_gg) +
  layer_spatial(aux_gg,aes(fill=Dif),col="transparent")+ #,col="transparent"
  scale_fill_gradientn(colors=rev(brewer.pal(8,"RdBu")),name=expression(italic(SR)[hat(Theta)](italic(bold(c)))),limits=c(-8*10^(-5),8*10^(-5)))+
  theme_bw()+
  ggtitle("Residual analysis")+
  theme(text = element_text(size = 24),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        title = element_text(colour = "black",face = "bold"),
        legend.position = "right")
