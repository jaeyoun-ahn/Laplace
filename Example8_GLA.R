###############################################################################
## Generalized Laplace Approximation : Example 8 

## written by Rosy Oh

## last updated : 29 August 2025

###############################################################################

source("function_GLA.R")

lamb_set= c(1,5,10)
varR_set = c(0.1,3,5)

y=0; d=1
y_distn="pois"; R_distn="lognorm"

sc_set = data.frame(lamb =rep(lamb_set, each=length(varR_set)),
                    varR = rep(varR_set, times=length(lamb_set)))
head(sc_set)
dim(sc_set)

res_ex8 = matrix(NA, nrow=nrow(sc_set), ncol=5)
for(i in 1:nrow(sc_set)){
  res_ex8[i,] = run_ex8(lamb=sc_set[i,1], nu=log(sc_set[i,2]+1), 
                        y=y, R_distn=R_distn, y_distn=y_distn, 
                        size=NULL, d=d)
}
colnames(res_ex8) = c("cred", "exact_m2", "GBEA_m2", "d_exact","d_GBEA")

round(cbind(sc_set, res_ex8),3)
