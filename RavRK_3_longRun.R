library(tidyverse)
library(pbmcapply)
library(ggpubr)
library(broom)
library(rstanarm)
library(shinystan)
library(bayesplot)

rm(list=ls())

# load data; standardise numerical covariates; create factors.
ravens <- read_csv("data/ravens_2021_02_03.csv")
ravensScaled <- ravens %>% 
                  mutate(O.dist = scale(O.dist), 
                         RK.dist = scale(RK.dist), 
                         season = factor(season), 
                         route = factor(route))
# fit model
if(F){
    m.re.nest.diag <- stan_glmer(y2 ~ 1 + O.dist*season + RK.dist*season + (1|season:route) + (-1+O.dist|season:route) + (-1+RK.dist|season:route),
                                   offset = log(w), 
                                   family = poisson, # poisson, binomial
                                   prior = normal(0,0.5,F),
                                   data = ravensScaled, 
                                   QR = T, 
                                   iter = 4000,
                                   init_r = 0.1,
                                   seed = 162, # sample(1:1000,1)
                                   cores = 12,
                                   chains = 12)
  #saveRDS(m.re.nest.diag, "results/m.re.nest.diag_longRun_2021_02_09.rds")
}

# Rhat diagnostic
fit <- readRDS("results/m.re.nest.diag_longRun_2021_02_09.rds")
fit$stan_summary[,"Rhat"] %>% unname %>% round(4) # all less than 0.001
