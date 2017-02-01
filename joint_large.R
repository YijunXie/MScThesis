setwd("/nfs/37zfs1-yijun.xie/yijun.xie/MScThesis/codes")
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

set.seed(2)
qs = seq(0.9,0.99,1e-4)
#res = matrix(0,ncol = 4, nrow = length(qs))
registerDoMC(8)
res = foreach(j = 1:length(qs), .combine = rbind)%dopar%{
  q = qs[j]
  # simulate data
  N = 80000
  burn_in = 5000
  phi = 0.3#0.95
  t_nu = 5
  epi_nu = 5
  omega = -5 # -0.5
  delta = 2 # 0.35
  
  err1 = rnorm(N)
  err2 = rstd(N,mean = 0,sd = delta,nu = t_nu) # eta_t
  x = rep(0,N)
  sig = rep(0,N)
  sig[1] = 0.01
  for(i in 2:N){
    sig[i] = sqrt(exp(omega+phi * log(sig[i-1]^2) + err2[i]))
    x[i] = sig[i]*err1[i]
  }
  
  x = x[-c(1:burn_in)]
  
  thred = quantile(x,q)
  w = which(x > thred)
  dw = diff(w)
  tsv = which(dw == 1)
  
  # simulate data
  
  err1 = rnorm(N)
  err2 = rnorm(N,mean = 0,sd = delta)
  #err2 = rstd(N,mean = 0,sd = 0.35,nu = t_nu) # eta_t
  x2 = rep(0,N)
  sig2 = rep(0,N)
  sig2[1] = 0.01
  for(i in 2:N){
    sig2[i] = sqrt(exp(omega+phi * log(sig2[i-1]^2) + err2[i]))
    x2[i] = sig2[i]*err1[i]
  }
  
  x2 = x2[-c(1:burn_in)]
  
  thred2 = quantile(x2,q)
  w2 = which(x2 > thred2)
  dw2 = diff(w2)
  nsv = which(dw2 == 1)
  # for GARCH
  
  
  phi = 0.3
  beta = 0.1
  omega = 1e-5
  
  err1 = rnorm(N)
  
  gx = rep(0,N)
  gsig = rep(0,N)
  gsig[1] = 0.01
  gx[1] = 0
  for(i in 2:N){
    gsig[i] = sqrt(omega + beta * gsig[i-1]^2 + phi * gx[i-1]^2)
    gx[i] = gsig[i]*err1[i]
  }
  
  gx = gx[-c(1:burn_in)]
  
  thred_g = quantile(gx,q)
  gw = which(gx > thred_g)
  dgw = diff(gw)
  gn = which(dgw == 1)
  
  n = length(x)
  c((1-q)*(1-q),length(tsv)/(n),
              length(nsv)/(n),length(gn)/(n))
}

plot(qs,res[,1],ylim = c(0,0.05),type = 'l')
lines(qs,res[,2],type = 'l',col = 'red')
lines(qs,res[,3],type = 'l',col = 'blue')
lines(qs,res[,4],type = 'l',col = 'green')

plot(qs,res[,3]/res[,3],ylim = c(0.5,5),type = 'l')
lines(qs,res[,2]/res[,3],type = 'l',col = 'red')
lines(qs,res[,4]/res[,3],type = 'l',col = 'green')

