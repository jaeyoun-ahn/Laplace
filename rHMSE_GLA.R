###############################################################################
## Generalized Laplace Approximation 
## Numerical analysis : HMSE 

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
d=1

################################################################################

I=10000; B=1

#-------------------------------------------
# gamma prior 

mse.g1 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "gamma",
                 y_distn=y_distn, method="gamma", size=size, d=d)
mse.g2 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "gamma", 
                 y_distn=y_distn, method="invgauss", size=size, d=d)
mse.g4 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "gamma", 
                 y_distn=y_distn, method="lognorm", size=size, d=d)
mse.g5 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "gamma", 
                 y_distn=y_distn, method="gammab", size=size, d=d)
mse.g6 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "gamma", 
                 y_distn=y_distn, method="invgaussb", size=size, d=d)

res_gamma_mse = cbind(mse.g1,mse.g2,mse.g4,mse.g5,mse.g6) #,mse.g.exact
res_gamma_hmse1 = res_gamma_mse%>%dplyr::select(starts_with("hmse1_m1"))%>% 
  rename_with(~ sub("^.*/", "", .x))
res_gamma_hmse2 = res_gamma_mse%>%dplyr::select(contains("hmse2_m1"))%>%
  rename_with(~ sub("^.*/", "", .x))       # "/” 전부 포함 앞부분 삭제

#-------------------------------------------
# lognormal prior 

mse.l1 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "lognorm", 
                 y_distn=y_distn, method="gamma", size=size, d=d)
mse.l2 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "lognorm", 
                 y_distn=y_distn, method="invgauss", size=size, d=d)
mse.l4 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "lognorm", 
                 y_distn=y_distn, method="lognorm", size=size, d=d)
mse.l5 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "lognorm", 
                 y_distn=y_distn, method="gammab", size=size, d=d)
mse.l6 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "lognorm", 
                 y_distn=y_distn, method="invgaussb", size=size, d=d)

res_ln_mse = cbind(mse.l1,mse.l2,mse.l4,mse.l5,mse.l6) #,mse.l.exact
head(res_ln_mse)
res_ln_hmse1 = res_ln_mse%>%dplyr::select(starts_with("hmse1_m1"))%>% 
  rename_with(~ sub("^.*/", "", .x))
res_ln_hmse2 = res_ln_mse%>%dplyr::select(contains("hmse2_m1"))%>%
  rename_with(~ sub("^.*/", "", .x))       # "/” 전부 포함 앞부분 삭제

#-------------------------------------------
# inverse gaussian prior 

mse.ig1 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn="invgauss", 
                  y_distn=y_distn, method="gamma", size=size, d=d)
mse.ig2 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn="invgauss", 
                  y_distn=y_distn, method="invgauss", size=size, d=d)
mse.ig4 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "invgauss", 
                  y_distn=y_distn, method="lognorm", size=size, d=d)
mse.ig5 = mse.run.iter(B=B,I=I, mat.param=mat.param, Rdistn="invgauss", 
                  y_distn=y_distn, method="gammab", size=size, d=d)
mse.ig6 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn="invgauss", 
                  y_distn=y_distn, method="invgaussb", size=size, d=d)

res_invgauss_mse = cbind(mse.ig1,mse.ig2, mse.ig4,mse.ig5,mse.ig6) #,mse.ig.exact
head(res_invgauss_mse)
res_invgauss_hmse1 = res_invgauss_mse%>%dplyr::select(starts_with("hmse1_m1"))%>% 
  rename_with(~ sub("^.*/", "", .x))
res_invgauss_hmse2 = res_invgauss_mse%>%dplyr::select(contains("hmse2_m1"))%>%
  rename_with(~ sub("^.*/", "", .x))       

#-------------------------------------------
# weibull prior 

mse.w1 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "weibull", 
                 y_distn=y_distn, method="gamma", size=size, d=d)
mse.w2 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "weibull", 
                 y_distn=y_distn, method="invgauss", size=size, d=d)
mse.w4 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "weibull", 
                 y_distn=y_distn, method="lognorm", size=size, d=d)
mse.w5 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "weibull", 
                 y_distn=y_distn, method="gammab", size=size, d=d)
mse.w6 = mse.run.iter(B=B, I=I, mat.param=mat.param, Rdistn = "weibull", 
                 y_distn=y_distn, method="invgaussb", size=size, d=d)

res_weibull_mse = cbind(mse.w1,mse.w2,mse.w4,mse.w5,mse.w6) #,mse.w.exact

head(res_weibull_mse)
res_weibull_hmse1 = res_weibull_mse%>%dplyr::select(starts_with("hmse1_m1"))%>% 
  rename_with(~ sub("^.*/", "", .x))
res_weibull_hmse2 = res_weibull_mse%>%dplyr::select(contains("hmse2_m1"))%>%
  rename_with(~ sub("^.*/", "", .x))      

################################################################################
# Table

var_x = mat.param[,2]^2/mat.param[,3]
res.summary(res=res_gamma_hmse1/var_x, method_name=my_method) 
res.summary(res=res_ln_hmse1/var_x, method_name=my_method) 
res.summary(res=res_invgauss_hmse1/var_x, method_name=my_method) 
res.summary(res=res_weibull_hmse1/var_x, method_name=my_method) 

res.summary(res=res_gamma_hmse2/var_x, method_name=my_method) 
res.summary(res=res_ln_hmse2/var_x, method_name=my_method) 
res.summary(res=res_invgauss_hmse2/var_x, method_name=my_method) 
res.summary(res=res_weibull_hmse2/var_x, method_name=my_method) 

################################################################################
# Figure

my_method = c("gamma","invgauss","lognorm", "gammab", "invgaussb")
my_cols = c("#000000",  "#009E73", "#56B4E9", "#E69F00",  
            "#F0E442",   "#D55E00", "#CC79A7")

#-----------------------------------------
var_x = mat.param[,2]^2/mat.param[,3]

all_res = list(gamma = res_gamma_hmse1, lognormal=res_ln_hmse1,
               invgauss = res_invgauss_hmse1, weibull=res_weibull_hmse1 )
all_res=lapply(all_res, function(x){x/var_x})
head(all_res)

get_figure2(all_res=all_res, mat.param=mat.param, cols=my_cols, 
            lamb_vals= set.lamb,
                ylab=quote(rHMSE[1]), fig_title= NULL, save=TRUE,
               save_name=paste0("Figures_GA/sim_",y_distn,"_hmse1_lambda"))

all_res2 = list(gamma = res_gamma_hmse2/EVar2.g[,2], 
                lognormal=res_ln_hmse2/EVar2.ln[,2],
                invgauss = res_invgauss_hmse2/EVar2.ig[,2], 
                weibull=res_weibull_hmse2/EVar2.w[,2] )

get_figure2(all_res=all_res2, mat.param=mat.param, cols=my_cols, 
            lamb_vals= set.lamb,
                ylab=quote(rHMSE[2]), fig_title= NULL, save=TRUE, 
                save_name=paste0("Figures_GA/sim_",y_distn,"_hmse2_lambda"))

