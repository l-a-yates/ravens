############################################################
##
## 1-D point pattern analysis of Raven and Roadkill data
##
## glm models of imhomogenoous linear point patterns
## Includes fixed and mixed effects models
##
## Authors: L Yates, M Fielding 
## Date: 20/10/2020
##
############################################################

library(spatstat)
library(tidyverse)
#library(pbmcapply)
#library(ggpubr)
#library(broom)
library(rstanarm)

rm(list=ls())

# load data hyperframe (see RavRK_0_prepare.R)
data <- readRDS("data/ravens_KF_hyper.rds")
data
KF_linpp <- data$linpp[[1]]
KF_linpp
# choose a site and season
#site <- "K-F"
seas <- "spring"

# subset the linear point pattern and extract/generate spatial covariates
ravens_lpp <- subset.lpp(KF_linpp, season == seas & type == "FR", select = F)
roadkill_lpp <- subset.lpp(KF_linpp, season == seas & type == "RK", select = F)
X_F <- data$X_F[[1]]
X_R <- as.linim(distfun.lpp(roadkill_lpp))

# plots of raven observations and covariate images
ravens_lpp %>% plot
X_F %>% plot
X_R %>% plot

# fit Poisson density models using spatstat
m.spat.0 <- lppm(ravens_lpp ~ 1, eps = 15)
m.spat.1 <- lppm(ravens_lpp ~ 1 + X_F + X_R, eps = 15)

# extract model frame
mf <- m.spat.1 %>% model.frame.lppm %>% rename(y1 = .mpl.Y, w = `(weights)`) %>% mutate(y2 = as.numeric(y1>0))
mf %>% as_tibble

# Point process Poisson models can be fit using glms, with canonical log link
# The quadrat 'sizes', w,  must incorporated into the fit. 
# Two methods:
# 1) Use w as prior weights with response y1 = 1/w or 0 for presence/absence (non-integer counts throw warnings but estimates are correct)
# 2) Use log(w) as an offset with response y2 = 1 or 0 presence/absence
# We use the offset method (2) here.
# The 'dummy' observations, y = 0, must be treated carefully when evaluating goodness of fit, see logLik.adj2 function above

m.glm.a.0 <- glm(y1~1, weights = w, family = poisson, data = mf)  # method a: use weights with non-integer response
m.glm.b.0 <- glm(y2~1, offset = log(w), family = poisson, data = mf) # method b: use offset with integer response
m.glm.a.0 %>%  summary %>% coef; m.glm.b.0 %>%  summary %>% coef; m.spat.0$fit

m.glm.a.1 <- glm(y1~ 1 + X_F + X_R, weights = w, family = poisson, data = mf)
m.glm.b.1 <- glm(y2~ 1 + X_F + X_R, offset = log(w), family = poisson, data = mf)

m.glm.a.1 %>%  summary %>% coef; m.glm.b.1 %>%  summary %>% coef; m.spat.1$fit

# Model scores (AIC) - modified calculations (Baddeley eq(9.56), p347)
# The adjustments depend only on the weights and the values of response variable, so they
#  do not affect model selection using relative score estimates such as delta-AIC.
AIC(m.spat.1) # spatstat value
# adjusted AIC definition for 'weighted' glm
AIC_adj_a <- function(fit) -2*(-1*(deviance(fit)/2 + sum(log(fit$prior.weights[(fit$y>0)])) + sum(fit$y>0)) - length(coef(fit)))
AIC_adj_a(m.glm.a.1)
# adjusted AIC definition for 'offset' glm
AIC_adj_b <- function(fit) -2*(-0.5*(deviance(fit) + 2*sum((fit$offset+1)*fit$y)) - length(coef(fit)))
AIC_adj_b(m.glm.b.1)


# Bayesian estimation using rstanarm
m.stan.0 <- stan_glm(y2~ 1, offset = log(w), family = poisson, data = mf, cores = 4)
m.stan.1 <- stan_glm(y2~ 1 + X_F + X_R, offset = log(w), family = poisson, data = mf, cores = 4)

# estimate approximate leave-one-out cross validation
#loo_c <- loo_compare(loo(m.stan.0, cores = 10),loo(m.stan.1, cores = 10))
matrix(c(0.0,0.0,-6.1,5.4),2,2, byrow = T, dimnames = list(c("m.stan.1","m.stan.0"),c("elpd_diff","se_diff"))) # manually stored results

m.stan.1 %>% plot("mcmc_dens")
m.stan.1$stan_summary


# fixed-effects model by season
seasons <- c(spring = "spring", summer = "summer", autumn = "autumn", winter = "winter")
ravens_seas <- lapply(seasons,  function(s){
                  subset.lpp(data$linpp[[site]], season == s & type == "FR", select = F)
                })
roadkill_seas <- lapply(seasons,  function(s){
  subset.lpp(data$linpp[[site]], season == s & type == "RK", select = F)
})

X_R_seas <- lapply(roadkill_seas, function(x){as.linim(distfun.lpp(x))})
X_F_seas <- data$X_F[[site]]
X_R_seas$spring

lppm_seas <- lapply(seasons, function(s){
  lppm(ravens_seas[[s]] ~ 1 + X_F + X_R, eps = 15, data=list(X_F = X_F_seas, X_R = X_R_seas[[s]]))
})

mf_seas <- lapply(lppm_seas, function(mf) mf %>% 
                    model.frame.lppm %>% 
                    rename(y1 = .mpl.Y, w = `(weights)`) %>% 
                    mutate(y2 = as.numeric(y1>0))) %>% 
  bind_rows(.id = "season") %>% 
  mutate(season = factor(season))

mf_seas %>% as_tibble

m.glm_seas <- glm(y2 ~ season + X_F*season + X_R*season, offset = log(w), family = poisson, data = mf_seas)
m.stan_seas <- stan_glm(y2 ~ season + X_F*season + X_R*season, offset = log(w), family = poisson, data = mf_seas, cores = 4)

m.stan_seas$stan_summary
m.stan_seas %>% plot("mcmc_dens")


