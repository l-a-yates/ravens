##---------------------------------------------------------------------------------
##
## Code to reproduce analyses from the manuscript:
## “Roadkill islands: carnivore extinction shifts seasonal use of roadside carrion 
## by generalist avian scavenger” by Fielding, Matthew; Buettel, Jessie C; Brook, 
## Barry; Stojanović, Dejan; Yates, Luke (2021)
##
##
## This script produce interval plots of MCMC results 
## Figure 5 - Parameter estimates on the log-linear scale. 
##
## Authors: L Yates, M Fielding 
## Date: 01/07/2020
##
##---------------------------------------------------------------------------------

#library(spatstat)
library(tidyverse)
library(pbmcapply)
library(ggpubr)
library(broom)
library(rstanarm)
library(shinystan)
library(bayesplot)

rm(list=ls())

m.re.nest.diag <- readRDS("results/m.re.nest.diag_longRun_2021_02_09.rds")

seasons <- c(spring = "spring", summer = "summer", autumn = "autumn", winter = "winter")

#-------------
# Prepare data
#-------------

m.df <- m.re.nest.diag %>% as.data.frame() %>% as_tibble()

# RK.dist (X_R) data
# account for treatment contrast to calculate seasonal intercepts
RK.df <- m.df %>% select(contains("RK.dist")) %>% 
  mutate(spring = `seasonspring:RK.dist` + RK.dist, `seasonspring:RK.dist` = NULL) %>% 
  mutate(summer = `seasonsummer:RK.dist` + RK.dist, `seasonsummer:RK.dist` = NULL) %>% 
  mutate(winter = `seasonwinter:RK.dist` + RK.dist, `seasonwinter:RK.dist` = NULL) %>% 
  mutate(autumn = RK.dist, RK.dist = NULL) %>% 
  select(-contains("sigma")) %>% 
  rename_at(vars(starts_with("b")), function(x) str_remove(x,fixed("b[RK.dist season:route:")) %>% str_remove(fixed("]")))

RK.sr <- lapply(seasons, function(seas) select(RK.df,contains(paste0(seas,":"))) %>% 
  mutate_all(function(x) x + RK.df[[seas]]) %>% 
  rename_all(~str_remove(., fixed(paste0(seas,":")))) %>% 
  gather("route","estimate")) %>% bind_rows(.id = "season")

# unlogged scale (multiplicative form) 
#RK.sr.exp <- RK.sr %>% mutate(estimate = exp(estimate/attributes(m.re.nest.diag$data$RK.dist)$`scaled:scale`))

# O.dist (X_F) data
# account for treatment contrast to calculate seasonal intercepts
m.df %>% select(contains("O.dist")) %>% names

O.df <- m.df %>% select(contains("O.dist")) %>% 
  mutate(spring = `O.dist:seasonspring` + O.dist, `O.dist:seasonspring` = NULL) %>% 
  mutate(summer = `O.dist:seasonsummer` + O.dist, `O.dist:seasonsummer` = NULL) %>% 
  mutate(winter = `O.dist:seasonwinter` + O.dist, `O.dist:seasonwinter` = NULL) %>% 
  mutate(autumn = O.dist, O.dist = NULL) %>% 
  select(-contains("sigma")) %>% 
  rename_at(vars(starts_with("b")), function(x) str_remove(x,fixed("b[O.dist season:route:")) %>% str_remove(fixed("]")))


O.sr <- lapply(seasons, function(seas) select(O.df,contains(paste0(seas,":"))) %>% 
                  mutate_all(function(x) x + O.df[[seas]]) %>% 
                  rename_all(~str_remove(., fixed(paste0(seas,":")))) %>% 
                  gather("route","estimate")) %>% bind_rows(.id = "season")

# unlogged scale (multiplicative form)
#O.sr.exp <- O.sr %>% mutate(estimate = exp(estimate/attributes(m.re.nest.diag$data$O.dist)$`scaled:scale`)) 

# intercept data
# account for treatment contrast to calculate seasonal intercepts
# adjust for the shift due to the standardisation of the spatial covariates RK.dist and O.dist
int.df <- m.df %>% select(-contains("dist"), -contains("b"), - contains("igma")) %>% 
  mutate(spring = `seasonspring` + `(Intercept)`, `seasonspring` = NULL) %>% 
  mutate(summer = `seasonsummer` + `(Intercept)`, `seasonsummer` = NULL) %>% 
  mutate(winter = `seasonwinter` + `(Intercept)`, `seasonwinter` = NULL) %>% 
  rename(autumn = `(Intercept)`) 

scale.O.dist <- attributes(m.re.nest.diag$data$O.dist)$`scaled:center`/attributes(m.re.nest.diag$data$O.dist)$`scaled:scale`
scale.RK.dist <- attributes(m.re.nest.diag$data$RK.dist)$`scaled:center`/attributes(m.re.nest.diag$data$RK.dist)$`scaled:scale`

int.list <- lapply(seasons, function(seas) lapply(levels(m.re.nest.diag$data$route), function(route){
  tibble(season = seas, route = route, 
         estimate = int.df[[seas]] - (RK.df[[paste0(seas,":",route)]])*scale.RK.dist - (O.df[[paste0(seas,":",route)]])*scale.O.dist)
}))

int.sr <- lapply(int.list, bind_rows) %>% bind_rows()

#------
# PLOT
#------

plot.alt <- function(X.sr) X.sr %>% mutate(route = fct_relevel(route, "F-E", "F-M", "K-F", "K-N", "H-C", "H-N", "T-N", "T-P")) %>% 
  mutate(season = fct_relevel(season, "autumn", "winter", "spring", "summer")) %>% 
  group_by(season,route) %>% summarise(mean = mean(estimate), 
                                            med = median(estimate), 
                                            q05 = quantile(estimate,probs = c(0.05)),
                                            q25 = quantile(estimate,probs = c(0.25)),
                                            q75 = quantile(estimate,probs = c(0.75)) ,
                                            q95 = quantile(estimate,probs = c(0.95))) %>% 
  ggplot() +
  geom_boxplot(aes(x = route, ymin = q05, lower = q25, middle = mean, upper = q75, ymax = q95, fill=season), 
                col = "grey10", width = 0.6, stat = "identity") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=14,face="bold")) +
  scale_fill_manual(values=c("orangered","royalblue4","olivedrab4","gold"))
                

RK.plot <- plot.alt(RK.sr) + 
  labs(subtitle = "Distance to roadkill") +  
  geom_hline(aes(yintercept = 0), lty = "dashed") + 
  ylim(c(-3,1.5)) 

O.plot <- plot.alt(O.sr) + 
  labs(subtitle = "Distance to farmland") +  
  geom_hline(aes(yintercept = 0), lty = "dashed") + 
  ylim(c(-3,1.5))

int.plot <- plot.alt(int.sr) + labs(subtitle = "Intercept")

# tiff("RAVRK_ParaEst.tiff", units="in", width=15, height=10, res=300, bg = "transparent")
ggarrange(int.plot,O.plot,RK.plot, nrow = 3, ncol = 1, common.legend = TRUE)
# dev.off()

