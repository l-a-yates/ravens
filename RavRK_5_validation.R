############################################################
##
## 1-D point pattern analysis of Raven and Roadkill data
##
## glm models of imhomogeneous linear point patterns
## Includes fixed and mixed effects models
## All models estimated via Bayesian MCMC methods (rstanarm)
##
## This script produces posterior predictive checks / model validation
##
## Authors: L Yates, M Fielding 
## Date: 14/07/2020
##
############################################################

library(tidyverse)
library(ggpubr)
library(rstanarm)
library(shinystan)
library(bayesplot)

rm(list=ls())

# load model
file_name <- "m.re.nest.diag_longRun_2021_02_09.rds"
model <- readRDS(paste0("results/",file_name))

# extract data
ravData <- model$data

# weight functions for empirical density estimation
w<- function(w) w/sum(w) # not used
wInv <- function(w) (1/w)/sum(1/w)
getWeights <- wInv

# generate/load posterior predictions
REPS <- 1000 # number of posterior draws and number of resamples for CSR
#p.draws <- posterior_predict(model, draws = REPS, offset = log(ravData$w)) %>% t %>% as_tibble %>% as.list
#saveRDS(p.draws, paste0("results/p.draws_",file_name))
p.draws <- readRDS(paste0("results/p.draws_",file_name))

# Function: validate_models
# Computes model validation summary functions and produces plots
# @p.draws: list of posterior draws
# @ravData: data
# @dcov: the covariate to be summarised, i.e., "RK.dist" or "O.dist"
# @xlabel: label for x-axis
# @route: subset of routes to be included (default is all routes) 
# @season: subset of seasons to be included (default is all seasons) 
validate_models <- function(p.draws, ravData, dcov, xlabel, route = NULL, season = NULL){
  
  # global scaling parameters
  s.centre <- attributes(ravData[[dcov]])$`scaled:center`
  s.scale <- attributes(ravData[[dcov]])$`scaled:scale`
  
  # subset indices, to perform validation by grouping variables
  if(!is.null(route)) ri.route <- ravData$route %in% route else ri.route <- 1 # subset by route
  if(!is.null(season)) ri.season <- ravData$season %in% season else ri.season <- 1 # subset by season
  rindex = ri.route & ri.season
  
  # subset data
  ravData <- ravData[rindex,]
  
  # density parameters
  bw <- 1 # bandwidth
  xStart <- min(ravData[[dcov]])
  xStop <- max(ravData[[dcov]])
  
  # observations
  oData <- ravData %>% filter(y1>0) 
  obs <- density(oData[[dcov]], weights = getWeights(oData$w), bw, n = 512, from = xStart, to = xStop)[['y']] #%>% cumsum()
  x <- density(oData[[dcov]], weights = getWeights(oData$w), bw, n = 512, from = xStart, to = xStop)[['x']]
  data.tib = tibble(x = x, obs = obs)
  
  # null model (complete spatial randomness)
  nrav <- ravData %>% filter(y1>0) %>% nrow
  sims <- replicate(REPS, {sData <-  sample_n(ravData, nrav, replace = T); density(sData[[dcov]], bw, weights = NULL, n = 512, from = xStart, to = xStop)[['y']]})
  null.mean <- apply(sims,1,mean)
  null.tib = as_tibble(sims) %>% bind_cols(x = x, .) %>% gather("sim","y",-1)
  
  # selected model (posterior draws)
  mod.tib = sapply(p.draws, function(x) {
    x <- x[rindex] # subset
    density(ravData[[dcov]][x>0], weights = x[x>0]/sum(x[x>0]), bw, n = 512, from = xStart, to = xStop)[['y']]
  }) %>% 
    as_tibble() %>% bind_cols(x = x, .) %>% gather("mod.sim","y",-1)
  
  # compute intervals
  null.ci <- null.tib %>% group_by(x) %>% summarise(lo = quantile(y, c(0.025)), hi = quantile(y, c(0.975)))
  mod.ci <- mod.tib %>% group_by(x) %>% summarise(lo = quantile(y, c(0.025)), hi = quantile(y, c(0.975)))
  adj <- null.mean # centre plots on the mean values of the null model
  
  null.tib %>% 
    ggplot(aes(x = x*s.scale + s.centre)) +
    geom_line(aes(y = y - adj, group = sim, col = "Null (CSR)"), alpha = 0.02) +
    geom_line(aes(y = y - adj, group = mod.sim, col = "Best Model"), alpha = 0.02, data = mod.tib) +
    geom_line(aes(y = obs - adj, col = "Data"), size =1, data = data.tib) +
    geom_line(aes(y = lo - adj, col = "Null (CSR)"), size = 1, data = null.ci, lty = "dashed") +
    geom_line(aes(y = hi - adj, col = "Null (CSR)"), size = 1, data = null.ci, lty = "dashed") +
    geom_line(aes(y = lo - adj, col = "Best Model"), size = 1, data = mod.ci, lty = "dashed") +
    geom_line(aes(y = hi - adj, col = "Best Model"), size = 1, data = mod.ci, lty = "dashed") +
    theme_classic()+
    scale_colour_manual(name = "Model", values = c("red", "black","blue")) +
    labs(title = paste("Model validation:",paste(route,collapse = ","),paste(season,collapse = ", ")), "  ", subtitle = paste(REPS,"posterior draws", "(centred on null)"), 
         x = xlabel, y = "density")
  
}

# plots
validate_models(p.draws, ravData, dcov = "RK.dist", xlabel = "Distance to roadkill")
validate_models(p.draws, ravData, dcov = "O.dist", xlabel = "Distance to farmland")
