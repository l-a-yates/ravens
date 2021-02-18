##--------------------------------------------------------------------------------
##
## Code to reproduce analyses from the manuscript:
## “Roadkill islands: carnivore extinction shifts seasonal use of roadside carrion 
## by generalist avian scavenger” by Fielding, Matthew; Buettel, Jessie C; Brook, 
## Barry; Stojanović, Dejan; Yates, Luke (2021)
##
## Prepare raven, roadkill and route data for 1D (linear) point-pattern analysis
##
## Create marked linear-point pattern (spatstat::lpp) objects
## Compute spatial covariate values
## Store all data objects in a spatstat::hyperframe
## Generate quadrature schemes and export the data file
##    ravens_.csv, used in the main analysis
##
## See https://doi.org/10.6084/m9.figshare.14036807.v1 for data metadata
##
## Author: Luke Yates
## Date created: 06/10/2020 (last edit 02/02/2021)
##
##--------------------------------------------------------------------------------

library(spatstat)
library(tidyverse)

rm(list=ls())

#----------
# Functions
#----------

# Creates a 1D linear network (linnet) object.
# Observation window dilated by 1000m to capture all points
createRoute <- function(route.site, dilation = 1000){
  owin.site <- dilation.owin(owin(c(min(route.site$x),max(route.site$x)), c(min(route.site$y),max(route.site$y))),dilation)  
  route.psp <- psp(route.site$x[-length(route.site$x)],route.site$y[-length(route.site$y)], 
                   route.site$x[-1], route.site$y[-1],owin.site)
  as.linnet(route.psp)
}

#------
# Main
#------

# route and habitat data
route_data <- read_csv("data/route_data_2021_02_02.csv")
routes <- route_data %>% split(route_data$route)
linnets <- routes %>% lapply(createRoute)

route_names <- names(routes); names(route_names) <- route_names

# distance to open land type, covariate image (X_F covariate)
dist_fun <- lapply(route_names, function(x){
  routes[[x]] %>% filter(habitat !="FF") %>% select(x,y) %>% lpp(L = linnets[[x]]) %>% distfun.lpp() %>% as.linim
})

# raven and roadkill observations
rav_rk_obs_data <- read_csv("data/rav_rk_obs_2021_02_03.csv")
rav_rk <- rav_rk_obs_data %>% split(rav_rk_obs_data$route) 

# create spatstat linear-point-pattern (lpp) objects
obs_lpp <- lapply(route_names, function(x) {
  select(rav_rk[[x]],x,y,season,type,hab) %>% 
    mutate(across(c(season,type,hab),factor)) %>% 
  lpp(linnets[[x]]) %>% unique
  })

# store objects in a spatstat hyperframe
ravens_hyper <- hyperframe(site = route_names, 
                            linpp = obs_lpp,
                            X_F = dist_fun)

#saveRDS(ravens_hyper,"data/ravens_hyper_2021_02_03.rds")

# create list of 4 hyperframes, by season
seasons <- c(spring = "spring", summer = "summer", autumn = "autumn", winter = "winter")
data <- lapply(seasons, function(seas){
  hf <- ravens_hyper[,1]
  rownames(hf) <- hf$site %>% as.character()
  hf$ravens.lpp <- with(ravens_hyper, subset(linpp , subset = season == seas & type == "FR", select = hab))
  hf$O.dist <- ravens_hyper$dfun
  hf$RK.dist <-  with(ravens_hyper, as.linim(distfun.lpp(subset(linpp , subset = season == seas & type == "RK", select = hab))))
  hf
})

# use spatstat::lppm to fit (full) saturated model
# this generates a quadrat scheme and corresponding model frame for each route-season pairs (32) 
# dummy spacing eps = 15 is equivalent to the default scheme, but has been set explicitly for reproducibility. 
lppm_fits <- lapply(seasons, function(seas) with(data[[seas]], try(lppm(ravens.lpp, ~ O.dist + RK.dist + 1,eps =15))))
model_frames <- lapply(lppm_fits, lapply, function(fit) fit %>% model.frame.lppm %>% rename(y = .mpl.Y, w = `(weights)`))

# combine model frames into single dataframe, define response variables
ravens <- lapply(model_frames, bind_rows, .id = "route") %>% bind_rows(.id = "season") %>% 
  mutate(y1 = y, y2 = as.numeric(y1>0), y = NULL) %>% 
  mutate(across(c(season,route), as.factor))

#write_csv(ravens, "data/ravens_2021_02_03.csv")



