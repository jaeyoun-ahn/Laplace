###############################################################################
## Generalized Laplace Approximation

## written by Rosy Oh

## last updated : 29 Aug 2025
###############################################################################

## package
# install.packages("Deriv")
# install.packages("actuar")
library(Deriv)
library(actuar)

library(dplyr)     # %>%  , arrange
library(tidyr)     # gather()
library(purrr)     # map2_dfr()
##########################################################################

# Y ~ Poisson(lambda)

# Y ~ NB(size=size,prob=size/(R*lambi+size) )
# E(Y|R)=lambi = size/(lambi*R+size)
# yi = rnbinom(Tyeark, size=size, prob=size/(R*lambi+size)) 
# pnbinom(q=d, size=size, prob=size/(r*lambi+size), lower.tail = FALSE)

# R distn : "lognorm",  "gamma", "invgauss"
# method: "gamma", "invgauss", "lognorm" "gammab", "invgaussb" 

# p: f_star, target pdf 
# q: approximated pdf 

##########################################################################
# target pdf f*(x)

logf_star = function(R, lamb, nu, y, distn, y_distn, size= NULL){
  
  inv_c = integrate(f_func, lower=0, upper=Inf, 
                    lamb=lamb, nu=nu, y=y, distn=distn, 
                    y_distn=y_distn, size=size)$value
  
  result = logf_func(R=R, lamb=lamb, nu=nu, y=y, distn=distn, 
                     y_distn=y_distn, size=size) - log(inv_c)
  return(result)
}

#----------------------------------------------------
# non normalized pdf : f(x)

f_func = function(R, lamb, nu, y, distn, y_distn, size= NULL){
  
  exp(logf_func(R, lamb, nu, y, distn, y_distn, size=size))  
}

#----------------------------------------------------
# log non normalized pdf : log f(x)

logf_func = function(R, lamb, nu, y, distn, y_distn, size= NULL){
  tt = length(y)
  
  if(distn=="lognorm"&y_distn=="pois"){
    logfY = -lamb*tt*R + sum(y)*log(R)  
    logfY - log(R)-log(R)/2-(log(R))^2/2/nu  
  }else if(distn=="gamma"&y_distn=="pois"){
    logfY = -lamb*tt*R + sum(y)*log(R)  
    logfY - nu*R+(nu-1)*log(R)
  }else if(distn=="invgauss"&y_distn=="pois"){
    logfY = -lamb*tt*R + sum(y)*log(R)  
    logfY - 1.5*log(R)-nu*R/2 - nu/2/R
  }else if(distn=="weibull"&y_distn=="pois"){
    logfY = -lamb*tt*R + sum(y)*log(R)  
    logfY + (nu-1)*log(R) - (R*gamma(1+1/nu))^nu
  }else if(distn=="lognorm"&y_distn=="nb"){
    logfY = -tt*size*log(lamb*R+size)+sum(y)*log(R)-sum(y)*log(lamb*R+size)
    logfY - log(R)-log(R)/2-(log(R))^2/2/nu  
  }else if(distn=="gamma"&y_distn=="nb"){
    logfY = -tt*size*log(lamb*R+size)+sum(y)*log(R)-sum(y)*log(lamb*R+size)
    logfY - nu*R+(nu-1)*log(R)
  }else if(distn=="invgauss"&y_distn=="nb"){
    logfY = -tt*size*log(lamb*R+size)+sum(y)*log(R)-sum(y)*log(lamb*R+size)
    logfY - 1.5*log(R)-nu*R/2 - nu/2/R
  }else if(distn=="weibull"&y_distn=="nb"){
    logfY = -tt*size*log(lamb*R+size)+sum(y)*log(R)-sum(y)*log(lamb*R+size)
    logfY + (nu-1)*log(R) - (R*gamma(1+1/nu))^nu
  }
}

#-----------------------------------------------------
# g(x) = log(f(x)) - h1(x)

g_func = function(R, lamb, nu, y, distn, y_distn, method, size= NULL){ # log xf(x)
  
  if(method=="invgam"|method=="lognorm"|method=="gamma"|method=="gammab"){ 
    # h1 = -log(R)
    logf_func(R=R, lamb=lamb, nu=nu, y=y, distn=distn, 
              y_distn=y_distn ,size=size) + log(R)
  }else if(method=="invgauss"|method=="invgaussb"){
    # h1 = -log(R)*3/2
    logf_func(R=R, lamb=lamb, nu=nu, y=y, distn=distn,
              y_distn=y_distn,size=size) + log(R)*3/2
  }else{ # laplace
    logf_func(R=R, lamb=lamb, nu=nu, y=y, distn=distn, 
              y_distn=y_distn,size=size)
  }
}

# the second derivative 
g_d2 = Deriv(g_func, "R", nderiv = 2)


################################################################################
# approximation methods

#----------------------------------------------------
## Laplace approximation

l_approx = function(R=NULL, lamb, nu, y, distn, y_distn, size= NULL,
                    out="d", simnum=1e6){
  #--------------
  # find mode
  x0 = optimize(g_func, lower=0, upper=10, maximum=TRUE,
                lamb=lamb, nu=nu, y=y, distn=distn,
                y_distn=y_distn, size=size, method="laplace")$maximum
  
  g2 = g_d2(x0,  lamb, nu, y, distn=distn, y_distn=y_distn, size=size,
            method="laplace")
  
  if(out=="d"){
    result =  dnorm(R, mean=x0, sd = sqrt(-1/g2), log=TRUE)
  }else if(out=="r"){
    result = rnorm(simnum, mean=x0, sd = sqrt(-1/g2))
  }else{ # mean 
    result = x0
  } 
  return(result)
}

#----------------------------------------------------
## Gamma approximation 

g_approx = function(R=NULL, lamb, nu, y, distn, y_distn, size= NULL,
                    out="d", simnum=1e6){
  #--------------
  # find mode 
  x0 = optimize(g_func, lower=0, upper=10, maximum=TRUE,
                lamb=lamb, nu=nu, y=y, distn=distn,
                y_distn=y_distn, size=size, method="gamma")$maximum # 0.7297409
  g2 = g_d2(x0,  lamb, nu, y, distn=distn, y_distn=y_distn, size=size,method="gamma")
  g1=0
  a1 = g1 + g2*x0
  a2 = -g2*x0^2
  
  if(out=="d"){
    result = dgamma(R, shape= a2, rate = -a1, log = TRUE) 
  }else if(out=="r"){
    result = rgamma(simnum, shape= a2, rate = -a1)  
  }else{
    result = a2/(-a1)  
  } 
  return(result)
}

#-------------------------------------------------------------
# Inverse Gaussian approximation

ig_approx = function(R=NULL, lamb, nu, y, distn,y_distn, size= NULL, 
                     out="d", simnum=1e6){ # g(x) = logxf(x)-h1
  #--------------
  # find mode 
  x0 = optimize(g_func, lower=0, upper=10, maximum=TRUE,
                lamb=lamb, nu=nu, y=y, distn=distn,
                y_distn=y_distn, size=size,method="invgauss")$maximum
  
  g2 = g_d2(x0,  lamb, nu, y, distn=distn, y_distn=y_distn, 
            size=size,method="invgauss")
  g1=0
  a1 = g1 + g2*x0/2
  a2 = g2*x0^3/2
  
  if(out=="d"){
    result = dinvgauss(R, mean= x0, shape = -2*a2, log=TRUE)
  }else if(out=="r"){
    result = rinvgauss(simnum, mean= x0, shape = -2*a2)
  }else{
    result = x0
  }
  return(result)
}

#-------------------------------------------------------------
# log-normal approximation

ln_approx = function(R=NULL, lamb, nu, y, distn, y_distn, size= NULL, 
                     out="d", simnum=1e6){ # g(x) = logxf(x)-h1
  
  x0 = optimize(g_func, lower=0, upper=10, maximum=TRUE,
                lamb=lamb, nu=nu, y=y, distn=distn,
                y_distn=y_distn, size=size, method="lognorm")$maximum
  g2 = g_d2(x0,  lamb, nu, y, distn=distn, y_distn=y_distn, size=size, method="lognorm")
  g1=0
  a1 = g1*x0 -x0^2*log(x0)*g2 - x0*log(x0)*g1
  a2 = g2*x0^2/2 + g1*x0/2
  
  if(out=="d"){
    result = dlnorm(R, meanlog= -a1/2/a2, sdlog = sqrt(-1/2/a2), log=TRUE)
  }else if(out=="r"){
    result = rlnorm(simnum, meanlog= -a1/2/a2, sdlog = sqrt(-1/2/a2))
  }else{
    result = exp(-a1/2/a2+1/2/(-2*a2))
  }
  return(result)
}


##########################################################################
## Buhlmann exponential family approximation

cred_func = function(lamb, nu, y, distn, y_distn, size= NULL){ 
  
  tyear = length(y)
  if(distn=="lognorm"){
    var_r = exp(nu)-1
  }else if(distn=="gamma"){
    var_r = 1/nu
  }else if(distn=="invgauss"){
    var_r = 1/nu
  }else if(distn=="weibull"){
    var_r = gamma(1+2/nu)/(gamma(1+1/nu)^2)-1
  }
  
  mu =  lamb   # E[E(Y|R)] = E[lamb*R] = lamb (E[R]=1)
  if(y_distn=="pois"){
    EVar = lamb  # E[Var(Y|R)] = E[lamb*R]
  }else if(y_distn=="nb"){
    EVar = lamb^2/size*(var_r+1)+lamb  # E[Var(Y|R)] = E[lamb^2*R^2/size+lamb*R]
  }
  VarE = lamb^2*var_r  # Var[E(Y|R)] = Var[lamb*R]  
  
  z = tyear*VarE/(tyear*VarE+EVar)
  cred= unname(z*mean(y)+(1-z)*mu)
  
  return(cred/lamb)
}

#----------------------------------------------------
## Gamma BEA

gb_approx = function(R=NULL, lamb, nu, y, distn, y_distn, size= NULL,
                     out="d", simnum=1e6){
  
  x0 = optimize(g_func, lower=0, upper=10, maximum=TRUE,
                lamb=lamb, nu=nu, y=y, distn=distn,
                y_distn=y_distn, size=size,method="gammab")$maximum
  g2 = g_d2(x0,  lamb, nu, y, distn=distn,y_distn=y_distn, 
            size=size, method="gammab")
  
  cred =  cred_func(lamb, nu, y, distn, y_distn, size=size)
  a2 = -g2*x0^2
  a1 = -(a2)/cred 
  
  if(out=="d"){
    result = dgamma(R, shape= a2, rate = -a1, log = TRUE) 
  }else if(out=="r"){
    result = rgamma(simnum, shape= a2, rate = -a1) 
  }else{
    result = cred
  }
  return(result)
}

#----------------------------------------------------
## Inverse Gaussian BEA

igb_approx = function(R=NULL, lamb, nu, y, distn, y_distn, size= NULL, 
                      out="d", simnum=1e6){
  
  x0 = optimize(g_func, lower=0, upper=10, maximum=TRUE,
                lamb=lamb, nu=nu, y=y, distn=distn,
                y_distn=y_distn, size=size,method="invgaussb")$maximum
  
  g2 = g_d2(x0,  lamb, nu, y, distn=distn, y_distn=y_distn, size=size,method="invgaussb")
  cred =  cred_func(lamb, nu, y, distn, y_distn, size=size)
  a2 = g2*x0^3/2
  
  if(out=="d"){
    result = dinvgauss(R, mean= cred, shape = -2*a2, log=TRUE)
  }else if(out=="r"){
    result = rinvgauss(simnum, mean= cred, shape = -2*a2)
  }else{
    result = cred
  }
  return(result)
}


################################################################################
# Total variation distance

TVd.run = function(B, mat.param, Rdistn, y_distn, method, size= NULL){
  # y_distn : "pois", "nb"
  pb <- txtProgressBar(min = 1, max = nrow(mat.param), style = 3)
  m1 = matrix(NA, nrow=nrow(mat.param), ncol=1)
  m2 = matrix(NA, nrow=nrow(mat.param), ncol=1)
  sd = matrix(NA, nrow=nrow(mat.param), ncol=1)
  
  for(ii in 1:nrow(mat.param)){
    
    Tyeark=unname(mat.param[ii,1])
    lambi=unname(mat.param[ii,2])
    var_rj=unname(mat.param[ii,3])
    
    set.seed(639)
    tvd = matrix(NA, nrow=B, ncol=1)
    
    n.error = 0
    for(i in 1:B){
      # distribution for R 
      if(Rdistn=="lognorm"){
        nuj = log(var_rj+1) 
        R = rlnorm(1, meanlog=-nuj/2, sdlog=sqrt(nuj))
      }else if(Rdistn=="gamma"){
        nuj= 1/var_rj
        R= rgamma(1, shape=nuj, rate=nuj)
      }else if(Rdistn=="invgauss"){
        nuj= 1/var_rj
        R= rinvgauss(1, mean=1, shape=nuj)
      }else if(Rdistn=="weibull"){
        nuj= uniroot(function(k, sigsq=var_rj){
          gamma(1+2/k)-(sigsq+1)*(gamma(1+1/k)^2)}, c(0.1,10))$root
        R= rweibull(1, shape=nuj, scale=1/gamma(1+1/nuj))
      }
      #----------------------
      ## generate y
      if(y_distn=="pois"){
        yi =rpois(Tyeark, lambda=R*lambi)  
      }else if(y_distn=="nb"){
        # E(Y|R)=lambi = size/(lambi*R+size)
        yi = rnbinom(Tyeark, size=size, prob=size/(R*lambi+size)) 
      }
      
      l1d_tmp = try(L1d(lambi=lambi, nuj=nuj,  yi=yi, 
                         Rdistn=Rdistn, y_distn=y_distn, 
                         size=size, method=method) )
      if('try-error' %in% class( l1d_tmp)) l1d_tmp = NA
      tvd[i,] = l1d_tmp/2
    }
    
    m1[ii,] = apply(tvd,2,mean)
    m2[ii,] = apply(tvd,2,median) # colMeans
    sd[ii,] = apply(tvd,2,sd)

    setTxtProgressBar(pb, ii) 
  }
  
  res= data.frame(m1,m2,sd) 
  colnames(res) = paste0(colnames(res),"/",method)
  return(res)
}

#---------------------------------------------------------
# L1 distance

L1d = function(lambi, nuj, yi, Rdistn, y_distn, method, size= NULL){
  
  if(method=="laplace"){
    ff = function(r, lambi, nuj, yi, Rdistn, y_distn, size){
      logp = logf_star(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                       y_distn=y_distn, size=size) # log f(x)
      logq = l_approx(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                      y_distn=y_distn, size=size)
      abs(exp(logp)-exp(logq))
    }
  }else if(method=="gamma"){
    ff = function(r, lambi, nuj, yi, Rdistn, y_distn, size){
      logp = logf_star(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                       y_distn=y_distn, size=size) # log f(x)
      logq = g_approx(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                      y_distn=y_distn, size=size)
      abs(exp(logp)-exp(logq))
    }
  }else if(method=="invgam"){
    ff = function(r, lambi, nuj, yi, Rdistn, y_distn, size){
      logp = logf_star(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                       y_distn=y_distn, size=size) 
      logq = invg_approx1(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                          y_distn=y_distn, size=size)
      abs(exp(logp)-exp(logq))
    }
  }else if(method=="invgauss"){
    ff = function(r, lambi, nuj, yi, Rdistn, y_distn, size){
      logp = logf_star(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                       y_distn=y_distn, size=size) 
      logq = ig_approx(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                       y_distn=y_distn, size=size)
      abs(exp(logp)-exp(logq))
    }
  }else if(method=="lognorm"){
    ff = function(r, lambi, nuj, yi, Rdistn, y_distn, size){
      logp = logf_star(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                       y_distn=y_distn, size=size) 
      logq = ln_approx(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                       y_distn=y_distn, size=size)
      abs(exp(logp)-exp(logq))
    }
  }else if(method=="gammab"){
    ff = function(r, lambi, nuj, yi, Rdistn, y_distn, size){
      logp = logf_star(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                       y_distn=y_distn, size=size) 
      logq = gb_approx(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                       y_distn=y_distn, size=size)
      abs(exp(logp)-exp(logq))
    }
  }else if(method=="invgaussb"){
    ff = function(r, lambi, nuj, yi, Rdistn, y_distn, size){
      logp = logf_star(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                       y_distn=y_distn, size=size) 
      logq = igb_approx(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                        y_distn=y_distn, size=size)
      abs(exp(logp)-exp(logq))
    }
  }
  res = integrate(ff, lower=1e-05, upper=Inf, lambi=lambi, nuj=nuj, yi=yi, 
                  Rdistn=Rdistn, y_distn=y_distn, size=size, rel.tol = 1e-8, 
                  stop.on.error=TRUE)$value
  
  names(res) = c(method)
  return(res)
}

################################################################################
# Relative hypothetical mean squared error

# expectation of random effect based on the approximating pdf eta 
# hat_E[X|Y_1:tau]
approx_m = function(lambi, nuj, yi, Rdistn, y_distn, method, size= NULL){
  if(method=="gamma"){
    m = g_approx(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                 y_distn=y_distn, size=size, out="m")
  }else if(method=="invgam"){
    m = invg_approx1(r, lamb=lambi, nu=nuj, y=yi, 
                          distn=Rdistn, size=size, out="m")
  }else if(method=="invgauss"){
    m = optimize(g_func, lower=0, upper=10, maximum=TRUE,
                    lamb=lambi, nu=nuj, y=yi, distn=Rdistn, y_distn=y_distn,
                    method="invgauss", size=size)$maximum 
  }else if(method=="lognorm"){
    m = ln_approx(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                  y_distn=y_distn, size=size, out="m")
  }else if(method=="gammab"){
    m =  cred_func(lamb=lambi, nu=nuj, y=yi, distn=Rdistn, y_distn=y_distn,
                   size=size)
  }else if(method=="invgaussb"){
    m =  cred_func(lamb=lambi, nu=nuj, y=yi, distn=Rdistn, y_distn=y_distn,
                   size=size)
  }
  names(m) = c(method)
  return(m)
}

#----------------------------------------------------------
# deductible adjusted premium
# E[(Y-d)+|X] with target pdf f*(x) 
# E[(Y-d)_+|R]=sum_{y>d} (y-d)*p(y), 

exact_m2_vec = function(R,lambi,y_distn=c("nb", "pois"),size = NULL,d){
  
  y_distn = match.arg(y_distn)      
  yy      = 0:d                     
  
  if(y_distn=="pois"){
    sapply(R, function(Ri) {
      mu  = lambi * Ri
      lev = sum(yy * dpois(yy, lambda = mu)) +
        d * (1 - ppois(d, lambda = mu))
      mu - lev
    })
  }else if(y_distn=="nb"){
    if (is.null(size)) stop("size must be supplied for negative binomial.")
    sapply(R, function(Ri) {
      mu  = lambi * Ri
      lev = sum(yy * dnbinom(yy, size=size, prob=size/(Ri*lambi+size))) +
        d * (1 - pnbinom(d,size=size, prob=size/(Ri*lambi+size)))
      mu - lev
    })
  }
}

#----------------------------------------------------------
approx_m2 = function(lambi, nuj, yi, Rdistn, y_distn, 
                     method, size= NULL, d){
  if(method=="gamma"){
    ff2 = function(r, lambi, nuj, yi, Rdistn, y_distn, size, d){
      logq = g_approx(R=r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                      y_distn=y_distn, size=size)
      sum_ee = exact_m2_vec(R=r, lamb=lambi, y_distn=y_distn, size=size, d=d) 
      sum_ee*exp(logq)
    }
  }else if(method=="invgam"){
    ff2 = function(r, lambi, nuj, yi, Rdistn, y_distn, size, d){
      logq = invg_approx1(R=r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                          y_distn=y_distn, size=size)
      sum_ee = exact_m2_vec(R=r, lamb=lambi, y_distn=y_distn, size=size, d=d) 
      sum_ee*exp(logq)
    }
  }else if(method=="invgauss"){
    ff2 = function(r, lambi, nuj, yi, Rdistn, y_distn, size, d){
      logq = ig_approx(R=r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                       y_distn=y_distn, size=size)
      sum_ee = exact_m2_vec(R=r, lamb=lambi, y_distn=y_distn, size=size, d=d) 
      sum_ee*exp(logq)
    }
  }else if(method=="lognorm"){
    ff2 = function(r, lambi, nuj, yi, Rdistn, y_distn, size, d){
      logq = ln_approx(R=r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                       y_distn=y_distn, size=size)
      sum_ee = exact_m2_vec(R=r, lamb=lambi, y_distn=y_distn, size=size, d=d) 
      sum_ee*exp(logq)
    }
  }else if(method=="gammab"){
    ff2 = function(r, lambi, nuj, yi, Rdistn, y_distn, size, d){
   
      logq = gb_approx(R=r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                       y_distn=y_distn, size=size)
      sum_ee = exact_m2_vec(R=r, lamb=lambi, y_distn=y_distn, size=size, d=d) 
      sum_ee*exp(logq)
    }
  }else if(method=="invgaussb"){
    ff2 = function(r, lambi, nuj, yi, Rdistn, y_distn, size, d){
    
      logq = igb_approx(R=r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                        y_distn=y_distn, size=size)
      sum_ee = exact_m2_vec(R=r, lamb=lambi, y_distn=y_distn, size=size, d=d) 
      sum_ee*exp(logq)
    }
  }

  m2 = integrate(ff2, lower=0, upper=Inf, lambi=lambi, nuj=nuj, yi=yi, 
                 Rdistn=Rdistn, y_distn=y_distn, size=size, d=d,
                 rel.tol = 1e-8, #.Machine$double.eps^0.5, #1e-8, #(23July2024)
                 stop.on.error=TRUE)$value
  
  names(m2) = c(method)
  return(m2)
}


#----------------------------------------------------------
# rHMSE
mse.run.iter = function(B, I, mat.param, Rdistn, y_distn, method, 
                        size= NULL, d=1, seednum=639){
  
  pb <- txtProgressBar(min = 1, max = nrow(mat.param), style = 3)
  
  hmse1_m1 = matrix(NA, nrow=nrow(mat.param), ncol=1)
  hmse1_m2 = matrix(NA, nrow=nrow(mat.param), ncol=1)
  hmse1_sd = matrix(NA, nrow=nrow(mat.param), ncol=1)
 
  hmse2_m1 = matrix(NA, nrow=nrow(mat.param), ncol=1)
  hmse2_m2 = matrix(NA, nrow=nrow(mat.param), ncol=1)
  hmse2_sd = matrix(NA, nrow=nrow(mat.param), ncol=1)
 
  for(ii in 1:nrow(mat.param)){
    
    Tyeark=unname(mat.param[ii,1])
    lambi=unname(mat.param[ii,2])
    var_rj=unname(mat.param[ii,3])
    
    set.seed(seednum)

    hmse1 = matrix(NA, nrow=B, ncol=1)
    hmse2 = matrix(NA, nrow=B, ncol=1)
    n.error = 0
    for(b in 1:B){
      # distribution for R : R : vector length=I
      if(Rdistn=="lognorm"){
        nuj = log(var_rj+1) 
        R_vec = rlnorm(I, meanlog=-nuj/2, sdlog=sqrt(nuj))
      }else if(Rdistn=="gamma"){
        nuj= 1/var_rj
        R_vec= rgamma(I, shape=nuj, rate=nuj)
      }else if(Rdistn=="invgauss"){
        nuj= 1/var_rj
        R_vec= rinvgauss(I, mean=1, shape=nuj)
      }else if(Rdistn=="weibull"){
        nuj= uniroot(function(k, sigsq=var_rj){
          gamma(1+2/k)-(sigsq+1)*(gamma(1+1/k)^2)}, c(0.1,10))$root
        R_vec= rweibull(I, shape=nuj, scale=1/gamma(1+1/nuj))
      }
      
      m = matrix(NA, nrow=I, ncol=1)  # E[X|Y_1:tau] 
      m2 = matrix(NA, nrow=I, ncol=1) # E[(X-d)+|Y_1:tau] 
      for(i in 1:I){
        if(y_distn=="pois"){
          yi =rpois(Tyeark, lambda=R_vec[i]*lambi)  
        }else if(y_distn=="nb"){
          yi = rnbinom(Tyeark, size=size, prob=size/(R_vec[i]*lambi+size)) 
        }
 
        m[i,]=approx_m(lambi=lambi, nuj=nuj, yi=yi,Rdistn=Rdistn, 
                       y_distn=y_distn, size=size, method=method)
        
        m2[i,]=approx_m2(lambi=lambi, nuj=nuj, yi=yi,
                         Rdistn=Rdistn, y_distn=y_distn, method=method,
                         size=size, d=d)
      }
  
      hmse1[b,] = mean((lambi*R_vec-lambi*m)^2) 
      
      exact_E = exact_m2_vec(R=R_vec, lambi=lambi,y_distn=y_distn, size=size, d=d)
      hmse2[b,] = mean((exact_E-m2)^2) 
    } # end iteration B
    
    hmse1_m1[ii,] = apply(hmse1,2,mean) 
    hmse1_m2[ii,] = apply(hmse1,2,median) 
    hmse1_sd[ii,] = apply(hmse1,2,sd)
 
    hmse2_m1[ii,] = apply(hmse2,2,mean) 
    hmse2_m2[ii,] = apply(hmse2,2,median) 
    hmse2_sd[ii,] = apply(hmse2,2,sd)
    setTxtProgressBar(pb, ii) 
  }
  
  res= data.frame(hmse1_m1,hmse1_m2,hmse1_sd,
                  hmse2_m1,hmse2_m2,hmse2_sd) 
  colnames(res) = paste0(colnames(res),"/",method);
  return(res)
}

##############################################################################
# Table

res.summary = function(res, method_name){
  result = cbind(mat.param, res)%>% 
    data.frame()%>%
    gather("method_name","measure",-Tyear, -lamb, -var_r)
  
  lamb_v = unique(result$lamb)
  #(1) lambda = 0.1
  list1=lapply(method_name, function(approx.){
    result%>%filter(lamb==lamb_v[1])%>%
      dplyr::select(-lamb)%>%
      filter(method_name==approx.,
             Tyear%in% c(1,5,10))%>%
      spread(Tyear, measure)
  })
  #(2) lambda = 0.5
  list2=lapply(method_name, function(approx.){
    result%>%filter(lamb==lamb_v[2])%>%
      dplyr::select(-lamb)%>%
      filter(method_name==approx.,
             Tyear%in% c(1,5,10))%>%
      spread(Tyear, measure)
  })
  #(3) lambda = 0.9
  list3=lapply(method_name, function(approx.){
    result%>%filter(lamb==lamb_v[3])%>%
      dplyr::select(-lamb)%>%
      filter(method_name==approx.,
             Tyear%in% c(1,5,10))%>%
      spread(Tyear, measure)
  })
  
  res = lapply(seq(length(list1)),
               function(i) cbind(list1[[i]], list2[[i]],list3[[i]]))
  return(res)
}

##############################################################################
# Figure 

get_figure2 = function(all_res, mat.param, lamb_vals,
                           cols, ylab=NULL, ylims=NULL, fig_title=NULL,
                           save=FALSE, save_name="HMSE_lambda",
                           width=10, height=6.5, dpi=300){
  
  tmp = map2_dfr( .x=all_res,           
                  .y=names(all_res),    
                  .f=\(res, prior_nm){ 
                    data.frame(cbind(mat.param, res)) %>%   
                      gather(key = "method", value = "measure",
                             -c(Tyear, lamb, var_r)) %>%   # wide → long
                      arrange(lamb, var_r, Tyear) %>%
                      mutate(prior = factor(prior_nm, 
                                            levels=c("gamma","lognormal",
                                                     "invgauss","weibull")),
                             method = factor(method, 
                                             levels=c("gamma","invgauss","lognorm",
                                                      "gammab", "invgaussb"),
                                             labels=c("gTEA","igTEA","lnTEA",
                                                      "gBEA","igBEA"))) })
  plots = purrr::map(lamb_vals, \(lv) {
    p=tmp%>%
      dplyr::filter(lamb==lv)%>%
      ggplot(aes(Tyear, measure, colour=method, linetype=method))+
      geom_line(aes(group=method), linewidth=1)+
      facet_grid( rows=vars(var_r), cols=vars(prior),labeller=label_both, 
                  scales='fixed')+ 
      scale_colour_manual(name="Methods",  values = cols)+
      scale_linetype_discrete(name="Methods")+
      scale_x_continuous(breaks=c(1:10))+ 
      labs(x=quote(tau), y=ylab, title=fig_title)+
      theme_bw()+
      theme(text = element_text(size=15))
    
    if (!is.null(ylims)) {
      p = p + lims(y = ylims)
    }
    
    if(save==TRUE){
      file_nm <- sprintf("%s_%.1f.pdf", save_name, lv)
      ggplot2::ggsave(filename = file_nm, plot = p,
                      width = width, height = height, dpi = dpi)
    }
    p
  })
  names(plots) =  paste0("lambda_", lamb_vals)
  invisible(plots) 
}

get_figure_TVd = function(res, cols, ylab=NULL, ylims=NULL , fig_title=NULL,
                          save=FALSE, save_name=NULL, x_breaks=NULL ){
  
  tmp = gather(res, "key"=method, value="measure", 
               -c(Tyear, lamb, var_r))%>%
    arrange(lamb, var_r, Tyear)%>%
    mutate(method = factor(method, 
                           levels=c("laplace","gamma","invgauss","lognorm",
                                    "gammab", "invgaussb"),
                           labels=c("Laplace","gTEA","igTEA","lnTEA",
                                    "gBEA","igBEA")))
  p=tmp%>%
    ggplot(aes(Tyear, measure, colour=method, linetype=method))+
    geom_line(aes(group=method), linewidth=1)+
    facet_grid( cols=vars(lamb), rows=vars(var_r),labeller=label_both, 
                scales='fixed')+ 
    scale_colour_manual(name="Methods",  values = cols)+
    scale_linetype_manual(name="Methods",
                          values=c("twodash","solid","dashed",
                                   "longdash","dotdash","dotted"))+
    scale_x_continuous(breaks=x_breaks)+ 
    labs(x=quote(tau), y=ylab, title=fig_title)+
    theme_bw()+
    theme(text = element_text(size=15))
  
  if (!is.null(ylims)) {
    p = p + lims(y = ylims)
  }
  
  if(save==TRUE){
    ggsave(save_name, width=10,height=6.5)
  }
  return(p)
}

##############################################################################
# Example 8

# exact 
post = function(r, lambi, nuj, yi, Rdistn, y_distn, size, d){
  logp = logf_star(r, lamb=lambi, nu=nuj, y=yi, distn=Rdistn, 
                   y_distn=y_distn, size=size) # log f(x)
  sum_ee = exact_m2_vec(R=r, lambi=lambi, y_distn=y_distn, size=size, d=d) 
  sum_ee*exp(logp)
}

#----------------------------------------------------------

run_ex8 = function(lamb, nu, y, R_distn, y_distn, size=NULL, d){
  
  cred = approx_m(lambi=lamb, nuj=nu, yi=y,Rdistn=R_distn, 
                  y_distn=y_distn, size=size, method="gammab")
  
  exact_m2 = integrate(post, lower=0, upper=Inf, lambi=lamb, nuj=nu, yi=y, 
                       Rdistn=R_distn, y_distn=y_distn, size=size, d=d,
                       rel.tol = 1e-8, stop.on.error=TRUE)$value
  
  GBEA_m2 = approx_m2(lambi=lamb, nuj=nu, yi=y, Rdistn=R_distn, y_distn=y_distn, 
                      method="gammab", size=size, d=d)
  
  res = c(cred*lamb, exact_m2, GBEA_m2, 
          cred*lamb-exact_m2, cred*lamb-GBEA_m2)
  names(res) = c("cred", "exact_m2", "GBEA_m2", "d_exact","d_GBEA")
  return(res)
}
