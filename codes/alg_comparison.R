setwd("/nfs/37zfs1-yijun.xie/yijun.xie/MScThesis/codes")
library("TTR")
library("MASS")
source("sv_funs.R")
require('doMC')
require('foreach')
require("fGarch")
require('rugarch')
#require('fitdistrplus')
require('fExtremes')
require('stochvol')
require("skewt")

s = Sys.time()

# Metropolis-Hastings for Normal-Std
mh.ht.old = function(o_ht,n_ht,x,b0,b1,sd,df,h2){
  u = runif(1)
  pn1 = dnorm(x,sd = sqrt(exp(n_ht)))
  pn2 = dstd(h2,mean =  b0 + b1*n_ht, sd = sd,nu = df)
  po1 = dnorm(x,sd = sqrt(exp(o_ht)))
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

mh.ht = function(o_ht,n_ht,x,b0,b1,sd,df,h2,oh2,old_theta){
  o_b0 = old_theta[1]
  o_b1 = old_theta[2]
  o_sd = old_theta[3]
  o_df = old_theta[4]
  u = runif(1)
  pn1 = dnorm(x,sd = sqrt(exp(n_ht)))
  pn2 = dstd(h2,mean =  b0 + b1*n_ht, sd = sd,nu = df)
  po1 = dnorm(x,sd = sqrt(exp(o_ht)))
  po2 = dstd(oh2,mean =  o_b0 + o_b1*o_ht, sd = o_sd,nu = o_df)
  a = (pn1*pn2)/(po1*po2)
  a = min(1,a)
  if(a > u){
    res = n_ht
  } else{
    res = o_ht
  }
  return(res)
}

ite = 3e3
burn_in = 2000
# ite = 10
# burn_in = 2

N = 4500

phi = 0.95
t_nu = 5
epi_nu = 5
omega = -0.5

registerDoMC(8)
s = Sys.time()
results = foreach(q = 1:150,.errorhandling = 'remove',.combine = rbind)%dopar%{
  set.seed(q)
  # generate simulated datasets
  #err1 = rnorm(N)# epsilon_t
  err1 = rnorm(N,mean = 0, sd = 1)
  # err2 = rskt(N,df = t_nu, gamma = 1.5) 
  # err2 = 0.35*(err2 - mean(err2))/sd(err2)
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
  old_hs = ths
  old_hs2 = ths
  old_theta = c(tb0,tb1,tetasd,tdf)
  intmr = data.frame(b0 = numeric(),b1 = numeric(),
                     sd = numeric(),df = numeric())
  for(l in 1:ite){
    n = length(ths)
    for(j in 2:(n-1)){
      o_ht = ths[j]
      n_ht = tb0 + tb1 * ths[j-1] + rstd(1,mean = 0,sd = tetasd,nu = tdf)
      # Normal-Std
      t_ht = mh.ht(o_ht,n_ht,dx[j], tb0,tb1,tetasd,tdf,
                   ths[j+1],ths[j+1],old_theta)
      ths[j] = t_ht
    }
    ths[1] = mean(ths)
    ths[length(ths)] = mean(ths)
    ths1 = ths[251:1750]
    afspec = arfimaspec(mean.model = list(armaOrder = c(1, 0)),
                        distribution.model = 'std')
    afm = arfimafit(afspec,ths1)
    old_theta = c(tb0,tb1,tetasd,tdf)
    tb1 = afm@model$pars[2,1]
    tb0 = (1-afm@model$pars[2,1])*afm@model$pars[1,1]
    tetasd = afm@model$pars[7,1]
    tdf = afm@model$pars[17,1]
    #if(l%%50 == 0){print(paste(l/50,"% finished"))}
    old_theta = c(tb0,tb1,tetasd,tdf)
    old_hs2 = old_hs
    old_hs = ths
    intmr = rbind(intmr,c(tb0,tb1,tetasd,tdf))
  }
  intmr = intmr[-c(1:burn_in),]
  c(median(intmr[,1]),median(intmr[,2]),median(intmr[,3]),
    median(intmr[,4]))
}

write.csv(results,'alg_com_1.csv')

results2 = foreach(q = 1:150,.errorhandling = 'remove',.combine = rbind)%dopar%{
  set.seed(q)
  # generate simulated datasets
  #err1 = rnorm(N)# epsilon_t
  err1 = rnorm(N,mean = 0, sd = 1)
  # err2 = rstd(N,df = t_nu, gamma = 1.5) 
  # err2 = 0.35*(err2 - mean(err2))/sd(err2)
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
  for(l in 1:ite){
    n = length(ths)
    for(j in 2:(n-1)){
      o_ht = ths[j]
      n_ht = tb0 + tb1 * ths[j-1] + rstd(1,mean = 0,sd = tetasd,nu = tdf)
      # Normal-Std
      # t_ht = mh.ht.old(o_ht,n_ht,dx[j], tb0,tb1,tetasd,tdf,
      #              ths[j+1])
      ths[j] = n_ht
    }
    ths[1] = mean(ths)
    ths[length(ths)] = mean(ths)
    ths1 = ths[251:1750]
    afspec = arfimaspec(mean.model = list(armaOrder = c(1, 0)),
                        distribution.model = 'std')
    afm = arfimafit(afspec,ths1)
    tb1 = afm@model$pars[2,1]
    tb0 = (1-afm@model$pars[2,1])*afm@model$pars[1,1]
    tetasd = afm@model$pars[7,1]
    tdf = afm@model$pars[17,1]
    if(l%%50 == 0){print(paste(l/50,"% finished"))}
    intmr = rbind(intmr,c(tb0,tb1,tetasd,tdf))
  }
  intmr = intmr[-c(1:burn_in),]
  c(median(intmr[,1]),median(intmr[,2]),median(intmr[,3]),
    median(intmr[,4]))
}

write.csv(results2,'alg_com_2.csv')
Sys.time()-s
