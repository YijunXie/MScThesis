## To westgrid
id = as.integer(Sys.getenv("PBS_ARRAYID"))
cat("This is array index: ", id, "\n")

library("TTR")
library("MASS")
source("new_sv_funs.R")
require('doMC')
require('foreach')
require("fGarch")
require('rugarch')
require("Hmisc")
require('fExtremes')
require('stochvol')
require("ismev")

mh.std = function(o_ht,n_ht,x,b0,b1,sd,df,h1,h2,e_nu){
  u = runif(1)
  pn1 = dstd(x,sd = sqrt(exp(n_ht)),nu = e_nu)
  pn2 = dstd(h2,mean =  b0 + b1*n_ht, sd = sd,nu = df)
  po1 = dstd(x,sd = sqrt(exp(o_ht)),nu = e_nu)
  po2 = dstd(h2,mean =  b0 + b1*o_ht, sd = sd,nu = df)
  a = (pn1*pn2)/(po1*po2)
  a = min(1,a)
  if(a > u){
    res = n_ht
  } else{
    res = o_ht
  }
  return(res)
}


## ARSV
ite = 5000
N = 8000
burn_in = 5000
phi = 0.95
t_nu = 5
epi_nu = 5
omega =  -0.5
delta = 0.35

err1 = rnorm(N)
err2 = rstd(N,mean = 0,sd = delta,nu = t_nu) # eta_t
all_x = rep(0,N)
sig = rep(0,N)
sig[1] = 0.01
for(i in 2:N){
  sig[i] = sqrt(exp(omega+phi * log(sig[i-1]^2) + err2[i]))
  all_x[i] = sig[i]*err1[i]
}

test_x = all_x[3001:5000]

x = all_x[-c(1:burn_in)]

theta_hat = c(-0.502,  0.950,  0.357, 5.37, 4.98)

theta_g = c(2.791076e-06, 2.208389e-01, 7.781607e-01)

alpha = 0.99
w = 1000
v = quantile(test_x,alpha)  # (length(x)-w)

job_per_core <- 15
s = proc.time()
dat = c()
for(i in 1:job_per_core){
  set.seed(1)
  tempdat = x[(i+id*15):(i+id*15+w-1)]
  #v = quantile(tempdat,alpha)
  ## garch
  spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE), 
                    distribution.model="std") # AR(1)-GARCH(1,1) model
  fit = ugarchfit(spec, tempdat)
  s_0 = fit@fit$sigma[w]
  
  gs_1 = sqrt(theta_g[1] + theta_g[3]*s_0^2 + theta_g[2]*tempdat[w]^2)
  cdfs = runif(1e6,alpha,1)
  ep_1 = qnorm(cdfs)
  gs_2 = sqrt(theta_g[1] + theta_g[3]*gs_1^2 + (theta_g[2]*gs_1^2)*ep_1^2)
  x_g = rnorm(1e6) * gs_2
  covar_g = quantile(x_g,alpha)
  
  ## arsv_estimated
  n = 0
  sn = matrix(nrow = 1e3,ncol = 2)
  while(n < 1e3){
    eta_1 = rstd(1,mean = 0,sd = theta_hat[3],nu = theta_hat[4])
    as_1 = sqrt(exp(theta_hat[1] + theta_hat[2] * log(s_0^2) + eta_1))
    ax_1 = as_1 * rstd(1, nu = theta_hat[5])
    if(ax_1 >= v){
      sn[n+1,] = c(as_1,eta_1)
      n = n+1
    }
  }
  eta_2 = rstd(1e3,mean = 0,sd = theta_hat[3],nu = theta_hat[4])
  as_2 = sqrt(exp(theta_hat[1] + theta_hat[2] * log(sn[,1]^2) + eta_2))
  x_a = as_2 * rstd(1e3, nu = theta_hat[5])
  covar_a = quantile(x_a,alpha)
  
  dat = rbind(dat,c(covar_g, covar_a ,i))
}
proc.time()-s

file <- paste("cvar_arsv_99",id,".RDS",sep="")
saveRDS(object = dat_all,file = file)
q("no")