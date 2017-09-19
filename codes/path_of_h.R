setwd("/home/xyj/ownCloud/Thesis/MScThesis")
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

ite = 3e3
burn_in = 2000
# ite = 10
# burn_in = 2

N = 4500

phi = 0.95
t_nu = 5
epi_nu = 5
omega = -0.5

set.seed(2)

# generate simulated datasets
err1 = rnorm(N,mean = 0, sd = 1)# epsilon_t
#err2 = rnorm(N, mean = 0, sd = 1)
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

#estimate realized volatility
rv = TTR::volatility(price,calc="close", n = 5) / sqrt(260)

dx = x[251:2250]
dsig = sig[251:2250]
drv = rv[252:2251]

#initial guess
hs = log(drv^2)
afspec = arfimaspec(mean.model = list(armaOrder = c(1, 0)),
                    distribution.model = 'std')
afm = arfimafit(afspec,hs)
tb1 = afm@model$pars[2,1] # beta_1
tb0 = (1-afm@model$pars[2,1])*afm@model$pars[1,1] # beta_0
tetasd = afm@model$pars[7,1] # delta
tdf = afm@model$pars[17,1] # eta_df
ths = hs
intmr = data.frame(b0 = numeric(),b1 = numeric(),
                   sd = numeric(),df = numeric())
hsample = c()
gbh = c()

for(l in 1:ite){
  n = length(ths)
  for(j in 2:(n-1)){
    o_ht = ths[j]
    n_ht = tb0 + tb1 * ths[j-1] + rstd(1,mean = 0,sd = tetasd,nu = tdf)
    if(j == 1750){gbh = c(gbh, n_ht)}
    # Normal-Std
    t_ht = mh.ht(o_ht,n_ht,dx[j],
                 tb0,tb1,tetasd,tdf,ths[j-1],ths[j+1])
    ths[j] = t_ht
  }
  ths[1] = mean(ths)
  ths[length(ths)] = mean(ths)
  ths1 = ths[251:1750]
  hsample = c(hsample,ths1[1500])
  afspec = arfimaspec(mean.model = list(armaOrder = c(1, 0)),
                      distribution.model = 'std')
  afm = arfimafit(afspec,ths1)
  tb1 = afm@model$pars[2,1]
  tb0 = (1-afm@model$pars[2,1])*afm@model$pars[1,1]
  tetasd = afm@model$pars[7,1]
  tdf = afm@model$pars[17,1]
  #tx = dx[251:1750]
  # Std-Std
  #z = tx / sqrt(exp(ths1)) # back-out first innovation
  #epn = fitdist(distribution = "std",z)$pars[3]
  #intmr = rbind(intmr,c(tb0,tb1,tetasd,tdf,epn))
  #if(l%%50 == 0){print(paste(l/50,"% finished"))}
  intmr = rbind(intmr,c(tb0,tb1,tetasd,tdf))
}

par(mfrow = c(2,2))
plot(intmr[,1],type = 'p',pch = ".",
     xlab = "",ylab = expression(beta[0]))
abline(h = -0.5,col = 'red')

plot(intmr[,2],type = 'p',pch = ".",xlab = "",ylab = expression(beta[1]))
abline(h = .95,col = 'red')

plot(intmr[,3],type = 'p',pch = ".",xlab = "",ylab = expression(delta))
abline(h = 0.35,col = 'red')

plot(intmr[,4],type = 'p',pch = ".",xlab = "",ylab = expression(nu))
abline(h = 5,col = 'red')