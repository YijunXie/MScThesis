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

s = Sys.time()

# Metropolis-Hastings for Normal-Std
mh.ht = function(o_ht,n_ht,x,b0,b1,sd,df,h1,h2){
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

# Metropolis-Hastings for Std-Std
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

ite = 5e3
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
results = foreach(q = 1:200,.errorhandling = 'remove',.combine = rbind)%dopar%{
  set.seed(q)
  
  # generate simulated datasets
  #err1 = rnorm(N)# epsilon_t
  err1 = rstd(N,mean = 0,sd = 1,nu = epi_nu)
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
  m1 = lm(hs[-1]~hs[-length(hs)])
  tb0 = m1$coefficients[1]
  tb1 = m1$coefficients[2]
  tresid = resid(m1)
  tetasd = sd(tresid)#fitdist(distribution = 'std',tresid)$pars[2]
  tdq = fitdist(distribution = 'std',tresid)
  tdf = tdq$pars[3]
  z = dx / sqrt(exp(hs))
  epn = fitdist(distribution = "std",z)$pars[3]
  ths = hs
  intmr = data.frame(b0 = numeric(),b1 = numeric(),
                       sd = numeric(),df = numeric(),
                     epn = numeric())
  for(l in 1:ite){
    n = length(ths)
    for(j in 2:(n-1)){
      o_ht = ths[j]
      n_ht = tb0 + tb1 * ths[j-1] + rstd(1,mean = 0,sd = tetasd,nu = tdf)
      # Normal-Std
      #t_ht = mh.ht(o_ht,n_ht,dx[j],
      #             tb0,tb1,tetasd,tdf,ths[j-1],ths[j+1])
      # Std-Std
      t_ht = mh.std(o_ht,n_ht,dx[j],tb0,tb1,tetasd,
                   tdf,ths[j-1],ths[j+1],epn)
      ths[j] = t_ht
    }
    ths[1] = mean(ths)
    ths[length(ths)] = tb0 + tb1 * ths[n-1] + rstd(1,mean = 0,sd = tetasd,nu = tdf)
    ths1 = ths[251:1750]
#     s1 = ths1[-1]
#     s2 = ths1[-length(ths1)]
#     m1 = lm(s1~s2)
#     tb0 = m1$coefficients[1]
#     tb1 = m1$coefficients[2]
#     tresid = resid(m1)
#     tdq = fitdist(distribution = 'std',tresid)
#     tetasd  = tdq$pars[2] 
#     tdf = tdq$pars[3]
    afspec = arfimaspec(mean.model = list(armaOrder = c(1, 0)),
                        distribution.model = 'std')
    afm = arfimafit(afspec,ths1)
    tb1 = afm@model$pars[2,1]
    tb0 = (1-afm@model$pars[2,1])*afm@model$pars[1,1]
    tetasd = afm@model$pars[7,1]
    tdf = afm@model$pars[17,1]
    tx = dx[251:1750]
    #################################################
    # Std-Std
    z = tx / sqrt(exp(ths1)) # back-out first innovation
    epn = fitdist(distribution = "std",z)$pars[3]
    intmr = rbind(intmr,c(tb0,tb1,tetasd,tdf,epn))
    ###################################################
    #if(l%%50 == 0){print(paste(l/50,"% finished"))}
    #intmr = rbind(intmr,c(tb0,tb1,tetasd,tdf))
  }
  intmr = intmr[-c(1:burn_in),]
  
  # draws = svsample(dx, draws = 3000,priornu = c(2,100),burnin = burn_in)
  # 
  # b1 = draws$summary$para[2,1]
  # b0 = draws$summary$para[1,1] * (1-b1)
  # svg = draws$summary$para[3,1]
  # snu = draws$summary$para[4,1]
  c(median(intmr[,1]),median(intmr[,2]),median(intmr[,3]),
    median(intmr[,4]),median(intmr[,5]))#,b0,b1,svg,snu,q)
}

write.csv(results,'nn_tt1.csv')

Sys.time()-s

r1 = read.csv("multi_mcmc_nt2.csv")
r1 = cbind(r1,"i"=seq(1:193))
#res = rbind(r1,r2)
res = rbind(nt1,nt2,nt3[,1:7])

par(mfrow = c(2,2))
hist(res[,1],col=rgb(1,0,0,0.5), main= 'b0', 
     xlim = c(-1.3,-0.05),
     breaks = 15,ylim = c(0,80),xlab = '')
hist(res[,5],col=rgb(0,0,1,0.5), add = T, breaks = 15)
abline(v = omega)
legend('topleft',c('MCMC','stochvol'),
       col = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),lwd = c(5,5))
hist(res[,2],col=rgb(1,0,0,0.5), main = 'b1', 
     xlim = c(0.87,0.99),
     breaks = 15,ylim = c(0,80),xlab = '')
hist(res[,6],col=rgb(0,0,1,0.5), add = T, breaks = 15)
abline(v = phi)
legend('topleft',c('MCMC','stochvol'),
       col = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),lwd = c(5,5))
hist(res[,3],col=rgb(1,0,0,0.5), main = 'sig', 
     xlim = c(0.15,0.66),
     breaks = 15,ylim = c(0,80),xlab = '')
hist(res[,7],col=rgb(0,0,1,0.5), add = T, breaks = 15)
abline(v = 0.35)
legend('topright',c('MCMC','stochvol'),
       col = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),lwd = c(5,5))
hist(res[,4],col=rgb(1,0,0,0.5), main = 'nu',
     xlim = c(0,25),
     breaks = 15,ylim = c(0,100),xlab = '')
#hist(res[,9],col=rgb(0,0,1,0.5), add = T, breaks = 200)
#abline(v = 0.35)
#hist(res[,5],col=rgb(0,1,1,0.5), add = T, breaks = 5)
abline(v = epi_nu)
legend('topright',c('MCMC'),
       col = c(rgb(1,0,0,0.5)),lwd = c(5,5))
legend('topright',c('MCMC','stochvol'),
       col = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),lwd = c(5,5))

# i = 1
# median(res[,i]);sd(res[,i]);i = i+1

Sys.time()-s

for(i in c(1,6,2,7,3,8,5,9)){
  print(median(res[,i]))
  print(sd(res[,i]))
}