setwd("/nfs/37zfs1-yijun.xie/yijun.xie/MScThesis/codes")
library("TTR")
library("MASS")
source("new_sv_funs.R")
require('doMC')
require('foreach')
require("fGarch")
require('rugarch')
#require('fitdistrplus')
require('fExtremes')
require('stochvol')
s = Sys.time()
set.seed(1)
ite = 8000
# generate simulation data
N = 4500
burn_in = 2500
phi = 0.95
t_nu = 5
epi_nu = 5
omega = -0.5

err1 = rnorm(N)
err2 = rstd(N,mean = 0,sd = 0.35,nu = t_nu)
x = rep(0,N)
sig = rep(0,N)
sig[1] = 0.01
for(i in 2:N){
  sig[i] = sqrt(exp(omega+phi * log(sig[i-1]^2) + err2[i]))
  x[i] = sig[i]*err1[i]
}

x = x[-c(1:burn_in)]
sig = sig[-c(1:burn_in)]
price = rep(0,(N-burn_in+1))
price[1] = 1
for(j in 2:(N-burn_in+1)){
  price[j] = price[j-1]*exp(x[j-1])
}

# initial guess
rv = TTR::volatility(price,calc="close", n = 5) / sqrt(260)

dx = x[251:2000]
dsig = sig[251:2000]
drv = rv[252:2001]

#initial guess
hs = log(drv^2)

afspec = arfimaspec(mean.model = list(armaOrder = c(1, 0)),
                    distribution.model = 'std')
afm = arfimafit(afspec,hs)
tb1 = afm@model$pars[2,1] # beta_1
tb0 = (1-afm@model$pars[2,1])*afm@model$pars[1,1] # beta_0
tetasd = afm@model$pars[7,1]
tdf = afm@model$pars[17,1]

ths = hs
ths_g = hs

intmr = data.frame(b0 = numeric(),b1 = numeric(),
                   sd = numeric(),df = numeric())

intmr_g = data.frame(b0 = numeric(),b1 = numeric(),
                   sd = numeric(),df = numeric())
tb1_g = tb1
tb0_g = tb0
tetasd_g = tetasd
tdf_g = tdf
# with M-H step
for(l in 1:ite){
  n = length(ths)
  for(j in 2:(n-1)){
    # with M-H
    o_ht = ths[j]
    n_ht = tb0 + tb1 * ths[j-1] + rstd(1,mean = 0,sd = tetasd,nu = tdf)
    # Normal-Std
    t_ht = mh.ht(o_ht,n_ht,dx[j],
                 tb0,tb1,tetasd,tdf,ths[j-1],ths[j+1])
    ths[j] = t_ht
  }
  
  ths1 = ths[201:1500] # truncated log(sigma^2)
  # fitting AR(1) with t innovations
  afspec = arfimaspec(mean.model = list(armaOrder = c(1, 0)),
                      distribution.model = 'std')
  # with M-H steps
  afm = arfimafit(afspec,ths1)
  tb1 = afm@model$pars[2,1] # beta_1
  tb0 = (1-afm@model$pars[2,1])*afm@model$pars[1,1] # beta_0
  tetasd = afm@model$pars[7,1] # delta
  tdf = afm@model$pars[17,1] # eta_df
  intmr = rbind(intmr,c(tb0,tb1,tetasd,tdf))
  
  n_g = length(ths_g)
  for(j in 2:(n_g-1)){
    # without M-H
    g_ht = tb0_g + tb1_g * ths_g[j-1] + rstd(1,mean = 0,sd = tetasd_g,nu = tdf_g)
    ths_g[j] = g_ht
  }
  
  # without M-H steps
  afm_g = arfimafit(afspec,ths_g)
  tb1_g = afm_g@model$pars[2,1] # beta_1
  tb0_g = (1-afm_g@model$pars[2,1])*afm_g@model$pars[1,1] # beta_0
  tetasd_g = afm_g@model$pars[7,1] # delta
  tdf_g = afm_g@model$pars[17,1] # eta_df
  intmr_g = rbind(intmr_g,c(tb0_g,tb1_g,tetasd_g,tdf_g))
}

Sys.time()-s

plot(intmr[,1],ylim = c(-20,0),type = 'l',xlab = "",ylab = "")
lines(intmr_g[,1],type = 'l',col = 'red')
abline(h = -0.5)
title(main=" ",ylab = "beta_0",xlab = "iterations")
legend("bottomleft",c("With M-H step","Without M-H step"), lty = c(1,1),
       lwd = c(2,2),col = c("black","red"))

plot(intmr[,2],ylim = c(-1,1),type = 'l',xlab = "",ylab = "")
lines(intmr_g[,2],type = 'l',col = 'red')
abline(h = 0.95)
title(main=" ",ylab = "beta_1",xlab = "iterations")
legend("bottomleft",c("With M-H step","Without M-H step"), lty = c(1,1),
       lwd = c(2,2),col = c("black","red"))
