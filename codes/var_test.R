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

set.seed(1)
ite = 5000

# simulate data
N = 8000
burn_in = 1000
phi = 0.9
t_nu = 5
omega = -1

err1 = rnorm(N)
err2 = rstd(N,mean = 0,sd = 0.35,nu = t_nu) # eta_t
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

# use first 2500 days to estimate variables
x_sample = x[1:2500]
p_sample = price[1:2501]
rv = TTR::volatility(p_sample,calc="close", n = 5) / sqrt(260)

dx = x[251:2250]
dsig = sig[251:2250]
drv = rv[252:2251]

#initial guess
hs = log(drv^2) # sigma^2
m1 = lm(hs[-1]~hs[-length(hs)])
tb0 = m1$coefficients[1] # beta_0
tb1 = m1$coefficients[2] # beta_1
tresid = resid(m1)
tetasd = sd(tresid)# delta
tdq = fitdist(distribution = 'std',tresid)
tdf = tdq$pars[3] # eta_df
ths = hs
intmr = data.frame(b0 = numeric(),b1 = numeric(),
                   sd = numeric(),df = numeric())

for(l in 1:ite){
  n = length(ths)
  for(j in 2:(n-1)){
    o_ht = ths[j]
    n_ht = tb0 + tb1 * ths[j-1] + rstd(1,mean = 0,sd = tetasd,nu = tdf)
    # Normal-Std
    t_ht = mh.ht(o_ht,n_ht,dx[j],
                tb0,tb1,tetasd,tdf,ths[j-1],ths[j+1])
    ths[j] = t_ht
  }
  ths[1] = mean(ths)
  ths[length(ths)] = tb0 + tb1 * ths[n-1] + rstd(1,mean = 0,sd = tetasd,nu = tdf)
  ths1 = ths[251:1750] # truncated log(sigma^2)
  # fitting AR(1) with t innovations
  afspec = arfimaspec(mean.model = list(armaOrder = c(1, 0)),
                      distribution.model = 'std')
  afm = arfimafit(afspec,ths1)
  tb1 = afm@model$pars[2,1] # beta_1
  tb0 = (1-afm@model$pars[2,1])*afm@model$pars[1,1] # beta_0
  tetasd = afm@model$pars[7,1] # delta
  tdf = afm@model$pars[17,1] # eta_df
  intmr = rbind(intmr,c(tb0,tb1,tetasd,tdf))
}
intmr = intmr[-c(1:2000),]

theta_hat = c(median(intmr[,1]),median(intmr[,2]),
              median(intmr[,3]), median(intmr[,4]))
