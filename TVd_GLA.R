###############################################################################
## Generalized Laplace Approximation 
## Numerical analysis : Total variation distance

## written by Rosy Oh

## last updated : 29 August 2025

#-----------------------------
# y_distn : "pois", "nb"

# negative binomial 
# Y|R ~ NB(lamb*R, size) , E[Y|R] = lamb*R
# p = size/(lamb*R+size)

# distribution for R :  E[R]=1
# gamma, lognormal, inverse gaussian, weibull

#-----------------------------
## setting 
# Tyear = 1:10
# lamb = c(0.1, 0.5, 0.9)
# var_r = c(0.3, 1.2, 2.0) # var_r = 1 : error  (x0 near zero)

###############################################################################
## package

library(dplyr)
library(tidyverse)
library(actuar)  # for inverse gaussian distn
library(ggplot2)
library(MASS)

source("function_GLA.r")

###############################################################################
## setting 

set.Tyear = 1:10
set.lamb = c(0.1, 0.5, 0.9)
set.var_r = c(0.3, 1.2, 2.0)

II=length(set.Tyear); J=length(set.lamb); K=length(set.var_r); 
mat.param = cbind(rep(set.Tyear, each=J*K),  # Tyear
                  rep(rep(set.lamb, each=K),by=II), #lambda
                  rep(set.var_r, by=II*J)) #Var(R)
colnames(mat.param) = c("Tyear","lamb","var_r")
head(mat.param)
dim(mat.param) # II*J*K, 3  # 90  3

y_distn="nb"; size=0.5;
# y_distn="pois"; size=NULL

################################################################################

B=3000

#-------------------------------------------
# gamma prior 

g0 = TVd.run(B=B, mat.param=mat.param, Rdistn = "gamma", 
             y_distn=y_distn,size=size, method="laplace")
g1 = TVd.run(B=B, mat.param=mat.param, Rdistn = "gamma", 
             y_distn=y_distn,size=size, method="gamma")
g2 = TVd.run(B=B, mat.param=mat.param, Rdistn = "gamma", 
             y_distn=y_distn,size=size, method="invgauss")
g4 = TVd.run(B=B, mat.param=mat.param, Rdistn = "gamma", 
             y_distn=y_distn,size=size, method="lognorm")
g5 = TVd.run(B=B, mat.param=mat.param, Rdistn = "gamma", 
             y_distn=y_distn,size=size, method="gammab")
g6 = TVd.run(B=B, mat.param=mat.param, Rdistn = "gamma", 
             y_distn=y_distn,size=size, method="invgaussb")

res_gamma_TVd = cbind(g0, g1,g2,g4,g5,g6) %>%
  dplyr::select(contains("m1"))%>%
  rename_with(~ sub("^.*/", "", .x))       # "/” 전부 포함 앞부분 삭제

#-------------------------------------------
# lognormal prior 

l0 = TVd.run(B=B, mat.param=mat.param, Rdistn = "lognorm", 
             y_distn=y_distn,size=size, method="laplace")
l1 = TVd.run(B=B, mat.param=mat.param, Rdistn = "lognorm", 
             y_distn=y_distn,size=size, method="gamma")
l2 = TVd.run(B=B, mat.param=mat.param, Rdistn = "lognorm", 
             y_distn=y_distn,size=size, method="invgauss")
l4 = TVd.run(B=B, mat.param=mat.param, Rdistn = "lognorm", 
             y_distn=y_distn, size=size, method="lognorm")
l5 = TVd.run(B=B, mat.param=mat.param, Rdistn = "lognorm", 
             y_distn=y_distn,size=size, method="gammab")
l6 = TVd.run(B=B, mat.param=mat.param, Rdistn = "lognorm", 
             y_distn=y_distn,size=size, method="invgaussb")

res_ln_TVd = cbind(l0,l1,l2,l4,l5,l6) %>%
  dplyr::select(starts_with("m1"))%>% 
  rename_with(~ sub("^.*/", "", .x))

#-------------------------------------------
# inverse gaussian prior 

ig0 = TVd.run(B=B, mat.param=mat.param, Rdistn = "invgauss", 
              y_distn=y_distn,size=size, method="laplace")
ig1 = TVd.run(B=B, mat.param=mat.param, Rdistn = "invgauss", 
              y_distn=y_distn,size=size, method="gamma")
ig2 = TVd.run(B=B, mat.param=mat.param, Rdistn = "invgauss", 
              y_distn=y_distn,size=size, method="invgauss")
ig4 = TVd.run(B=B, mat.param=mat.param, Rdistn = "invgauss", 
              y_distn=y_distn,size=size, method="lognorm")
ig5 = TVd.run(B=B, mat.param=mat.param, Rdistn = "invgauss", 
              y_distn=y_distn,size=size, method="gammab")
ig6 = TVd.run(B=B, mat.param=mat.param, Rdistn = "invgauss", 
              y_distn=y_distn,size=size, method="invgaussb")

res_invgauss_TVd =  cbind(ig0,ig1,ig2, ig4,ig5,ig6) %>%
  dplyr::select(starts_with("m1"))%>% 
  rename_with(~ sub("^.*/", "", .x))

#-------------------------------------------
# weibull prior 

w0 = TVd.run(B=B, mat.param=mat.param, Rdistn = "weibull", 
             y_distn=y_distn,size=size, method="laplace")
w1 = TVd.run(B=B, mat.param=mat.param, Rdistn = "weibull", 
             y_distn=y_distn,size=size, method="gamma")
w2 = TVd.run(B=B, mat.param=mat.param, Rdistn = "weibull", 
             y_distn=y_distn,size=size, method="invgauss")
w4 = TVd.run(B=B, mat.param=mat.param, Rdistn = "weibull", 
             y_distn=y_distn,size=size, method="lognorm")
w5 = TVd.run(B=B, mat.param=mat.param, Rdistn = "weibull", 
             y_distn=y_distn,size=size, method="gammab")
w6 = TVd.run(B=B, mat.param=mat.param, Rdistn = "weibull", 
             y_distn=y_distn,size=size, method="invgaussb")

res_weibull_TVd = cbind(w0,w1,w2,w4,w5,w6)%>%
  dplyr::select(starts_with("m1"))%>% 
  rename_with(~ sub("^.*/", "", .x))

################################################################################
# Table

res.summary(res=res_gamma_TVd, method_name=my_method) 
res.summary(res=res_ln_TVd, method_name=my_method)    
res.summary(res=res_invgauss_TVd, method_name=my_method) 
res.summary(res=res_weibull_TVd, method_name=my_method) 

################################################################################
# Figure

my_method = c("laplace","gamma","invgauss","lognorm", "gammab", "invgaussb")
my_cols = c("#999999","#000000",  "#009E73", "#56B4E9", "#E69F00",  
            "#F0E442",   "#D55E00", "#CC79A7")

get_figure_TVd(res=data.frame(cbind(mat.param, res_gamma_TVd)),
               x_breaks=set.Tyear, cols=my_cols, ylab="TV distance", 
               fig_title="Gamma prior", save=TRUE, 
               save_name=paste0("Figures_GA/sim_TVd_",y_distn,"_gamma.pdf"))
get_figure_TVd(res=data.frame(cbind(mat.param, res_ln_TVd)),
               x_breaks=set.Tyear, cols=my_cols, ylab="TV distance", 
               fig_title="Log-normal prior", save=TRUE, 
               save_name=paste0("Figures_GA/sim_TVd_",y_distn,"_ln.pdf"))

get_figure_TVd(res=data.frame(cbind(mat.param, res_invgauss_TVd)),
               x_breaks=set.Tyear, cols=my_cols, ylab="TV distance", 
               fig_title="Inverse Gaussian prior", save=TRUE, 
               save_name=paste0("Figures_GA/sim_TVd_",y_distn,"_ig.pdf"))
get_figure_TVd(res=data.frame(cbind(mat.param, res_weibull_TVd)),
               x_breaks=set.Tyear, cols=my_cols, ylab="TV distance", 
               fig_title="Weibull prior", save=TRUE, 
               save_name=paste0("Figures_GA/sim_TVd_",y_distn,"_weibull.pdf"))


