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
setwd("/nfs/37zfs1-yijun.xie/yijun.xie/MScThesis/codes")
set.seed(1)
ite = 5000
s = Sys.time()
sp <- read.csv("sp.csv", row.names=1)
dat = apply(sp,2,rev)
p = dat[,6]
all_x = diff(log(p))

test_x = all_x[1:2021]
price = p[1:2022]

rv = TTR::volatility(price,calc="close", n = 5) / sqrt(260)
drv = rv[-c(1:101)]
dx = test_x[-c(1:100)]
#initial guess
hs = log(drv^2) # sigma^2
# m1 = lm(hs[-1]~hs[-length(hs)])
# tb0 = m1$coefficients[1] # beta_0
# tb1 = m1$coefficients[2] # beta_1
# tresid = resid(m1)
# tetasd = sd(tresid)# delta
# tdq = fitdist(distribution = 'std',tresid)
# tdf = tdq$pars[3] # eta_df
registerDoMC(8)
est_paras = foreach(i = 129:200,.combine = rbind, .errorhandling = 'remove')%dopar%{
  set.seed(i)
  afspec = arfimaspec(mean.model = list(armaOrder = c(1, 0)),
                      distribution.model = 'norm')
  afm = arfimafit(afspec,hs)
  tb1 = afm@model$pars[2,1] # beta_1
  tb0 = (1-afm@model$pars[2,1])*afm@model$pars[1,1] # beta_0
  tetasd = afm@model$pars[7,1] # delta

  intmr = data.frame(b0 = numeric(),b1 = numeric(),
                     sd = numeric())
  ths = hs
  for(l in 1:ite){
    n = length(ths)
    for(j in 2:(n-1)){
      o_ht = ths[j]
      n_ht = tb0 + tb1 * ths[j-1] + rnorm(1,mean = 0,sd = tetasd)
      # Normal-Std
      t_ht = mh.dn(o_ht,n_ht,dx[j],
                   tb0,tb1,tetasd,ths[j-1],ths[j+1])
      ths[j] = t_ht
    }
    #ths[1] = mean(ths)
    #ths[length(ths)] = tb0 + tb1 * ths[n-1] + rnorm(1,mean = 0,sd = tetasd)
    ths1 = ths[201:1700] # truncated log(sigma^2)
    # fitting AR(1) with t innovations
    afspec = arfimaspec(mean.model = list(armaOrder = c(1, 0)),
                        distribution.model = 'norm')
    afm = arfimafit(afspec,ths1)
    tb1 = afm@model$pars[2,1] # beta_1
    tb0 = (1-afm@model$pars[2,1])*afm@model$pars[1,1] # beta_0
    tetasd = afm@model$pars[7,1] # delta
    intmr = rbind(intmr,c(tb0,tb1,tetasd))
  }
  intmr = intmr[-c(1:15000),]

  # estimated parameters for sv(nn) model
  theta_hat1 = c(median(intmr[,1]),median(intmr[,2]),
                 median(intmr[,3]))
  
  set.seed(i)
  ## ARSV-NT
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
      t_ht = mh.ht(o_ht,n_ht,dx[j],
                   tb0,tb1,tetasd,tdf,ths[j-1],ths[j+1])
      ths[j] = t_ht
    }
    #ths[1] = mean(ths)
    #ths[length(ths)] = tb0 + tb1 * ths[n-1] + rstd(1,mean = 0,sd = tetasd,nu = tdf)
    ths1 = ths[201:1700] # truncated log(sigma^2)
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
  intmr = intmr[-c(1:15000),]
  
  # estimated parameters for sv(nt) model
  theta_hat2 = c(median(intmr[,1]),median(intmr[,2]),
                 median(intmr[,3]), median(intmr[,4]))
  
  
  # stochvol package
  
  d = svsample(test_x)
  b1 = d$summary$para[2,1]
  b0 = d$summary$para[1,1] * (1-b1)
  svg = d$summary$para[3,1]
  
  theta_hat3 = c(b0, b1, svg)
  
  c(theta_hat1, theta_hat2, theta_hat3 ,i)
}

Sys.time()-s
write.csv(est_paras,'sp_paras_4c.csv')

paras = paras[,-11]
w1 = paras[order(paras[,1]),1:3]
np1 = w1[c(-c(1:15),-c(112:126)),]
w2 = paras[order(paras[,4]),4:7]
np2 = w2[c(-c(1:15),-c(112:126)),]
