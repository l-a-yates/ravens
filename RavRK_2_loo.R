##---------------------------------------------------------------------------------
##
## Code to reproduce analyses from the manuscript:
## “Roadkill islands: carnivore extinction shifts seasonal use of roadside carrion 
## by generalist avian scavenger” by Fielding, Matthew; Buettel, Jessie C; Brook, 
## Barry; Stojanović, Dejan; Yates, Luke (2021)
##
## 1-D point-pattern analysis of Raven and Roadkill data
##
## glm models of imhomogeneous linear point patterns
##
## LOO script
## Produces: model selection plot (Fig 3); complexity estimates (Table 1)
#
## Authors: L Yates, M Fielding 
## Date: 20/10/2020
##
##---------------------------------------------------------------------------------

library(tidyverse)
library(pbmcapply)
library(ggpubr)
library(rstanarm)
library(bayesplot)

rm(list=ls())

# load model fits
res_dir <- "results/"
run_date <- "2021_02_09"
model_names <- c("m.null", "m.pool", "m.fe.seas", "m.unpool", "m.fe.cross", "m.re.cross.cov", 
                 "m.re.cross.diag", "m.re.nest.cov", "m.re.nest.diag")
names(model_names) <- model_names
model_fits <- lapply(model_names, function(model_name){
  readRDS(paste0(res_dir,"fit_",model_name,"_",run_date,".rds"))
})

# compute and save loo estimates
ncores <- 40
#loo_est <- lapply(model_fits, rstanarm::loo, cores = ncores)
#saveRDS(loo_est,paste0(res_dir,"loo_est_",run_date,".rds"))
loo_est <- readRDS(paste0(res_dir,"loo_est_",run_date,".rds"))

# model comparison
looc <- loo_compare(loo_est)
loo_df <- looc %>% as.data.frame()
loo_tb <- tibble(model = rownames(loo_df), looic = loo_df$elpd_diff*-2, se_diff = loo_df$se_diff*-2)
loo_tb$modelno <- c("8", "9", "6", "7", "4", "5", "3", "2", "1")
loo_tb <- loo_tb %>% arrange(modelno)
loo_tb

# Figure 3 - Estimates of model performance using approximate leave-one-out cross validation
# tiff("RAVRK_ModSel3.tiff", units="in", width=4, height=2, res=300, bg = "transparent")
loo_tb %>% ggplot(aes(x = reorder(modelno, desc(modelno)), y = looic)) +
  geom_linerange(aes(ymin = looic - se_diff, ymax = looic + se_diff)) +
  geom_point() +
  coord_flip() +
  theme_classic() +
  labs(x = "Model", y= "delta-LOO")+
  theme(axis.text.x = element_text(size = 8)) 
# dev.off()    

# complexity estimates / number of estimated parameters
sapply(loo_est[str_detect(names(loo_est),"re")], function(x)x$estimates[2,1]) %>% round(1)
model_fits[1:5] %>% sapply(function(x) x %>% coef %>% length)
