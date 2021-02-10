############################################################
##
## 1-D point pattern analysis of Raven and Roadkill data
##
## glm models of imhomogeneous linear point patterns
## Includes fixed and mixed effects models
## All models estimated via Bayesian MCMC (HMC) methods (rstanarm)
##
## Authors: L Yates, M Fielding 
## Date: 20/10/2020
##
############################################################

#library(spatstat)
library(tidyverse)
library(pbmcapply)
library(ggpubr)
library(rstanarm)
library(bayesplot)

rm(list=ls())


#-----
# Data
#-----
# ravens data are the combined model matrices from (the default) linear quadrat schemes generated 
#   for each route using spatstat's lppm (see RavRK_0_prepare.R)
# Each row is a single quadrat in a given route for a given season 
# y2: 0 or 1, raven observed in quadrat
# w: weights = area of quadrat (proportional to quadrat length)
# y1 = y2/w (intensity/rate = raven/metre in quadrat)
# O.dist (numerical covariate): distance to open farmland from quadrat (X_F)
# RK.dist (numerical covariate): distance to roadkill from quadrat (X_R)
# route (group factor): 8 routes where observations were made
# season (group factor): 4 seasons, observations taken for each route in each season

#ravens <- read_csv("data/ravens_2020_07_01.csv")
ravens <- read_csv("data/ravens_2021_02_03.csv")

# standardise numerical covariates O.dist and RK.dist
ravensScaled <- ravens %>% mutate(O.dist = scale(O.dist), RK.dist = scale(RK.dist))

#--------------------
# Methods
#--------------------

# Point process Poisson models can be fit using standard glms, with Poisson errors and canonical log link
# The quadrat 'sizes' (effort), w,  must incorporated into the fit. 
# Two methods:
# 1) Use w as prior weights with response y1 = 1/w or 0 for presence/absence (non-integer counts may throw warnings but estimates are correct)
# 2) Use log(w) as an offset with response y2 = 1 or 0 presence/absence [we use this]
#
# Viewing the observations for each route in each season as replicated experiments, we use mixed effects models to investigate the role of 
#   O.dist and RK.dist. 
#
# We treat season (4 types) as a fixed effect and route (8 types) as a random effect. 
#
# Here we use Bayesian Monte-Carlo methods to estimate model parameters  

#--------
# Models
#--------

# Note: scaled predictors work much better; i.e. data = ravensScaled 
# using a global intercept (+1) reduces run times

# specify and run models

  models <- list(m.null = y2 ~ 1,
                 m.pool = y2 ~ O.dist + RK.dist + 1,
                 m.fe.seas = y2 ~ O.dist:season + RK.dist:season + season + 1,
                 m.unpool = y2 ~ 1 + O.dist*season*route + RK.dist*season*route, 
                 m.fe.cross = y2 ~ 1 + O.dist*season + RK.dist*season + O.dist*route + RK.dist*route,
                 m.re.cross.cov = y2 ~ 1 + O.dist*season + RK.dist*season + (1 + O.dist + RK.dist|route), 
                 m.re.cross.diag = y2 ~ 1 + O.dist*season + RK.dist*season + (1|route) + (-1+O.dist|route) + (-1+RK.dist|route),
                 m.re.nest.cov = y2 ~ 1 + O.dist*season + RK.dist*season + (1 + O.dist + RK.dist|season:route),
                 m.re.nest.diag = y2 ~ 1 + O.dist*season + RK.dist*season + (1|season:route) + (-1+O.dist|season:route) + (-1+RK.dist|season:route)
                  )
  seed <- 752 # sample(1:1000,1)
  
  # select the model to be run
  set_model <- 9
  
  if(set_model == 1) QR <- F else QR <-  T
  message(paste("Fitting model:", names(models)[[set_model]]))

  fit <- try(stan_glmer(models[[set_model]],   # use stan_glm/stan_glmer for models without/with random effects
                 offset = log(w), 
                 prior = normal(0,0.5,F), 
                 family = poisson, 
                 data = ravensScaled, 
                 init_r = 1,
                 QR = QR,
                 seed = seed,
                 cores = 4))

  #saveRDS(fit,paste0("results/fit_",names(models)[[set_model]],"_2021_02_09.rds"))

