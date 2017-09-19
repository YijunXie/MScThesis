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
burn_in = 2000
s = Sys.time()
sp <- read.csv("sp.csv", row.names=1)
dat = apply(sp,2,rev)
p = dat[,6]
all_x = -diff(log(p))

test_x = all_x[1:2021]
price = p[1:2022]

rv = TTR::volatility(price,calc="close", n = 5) / sqrt(260)
drv = rv[-c(1:101)]
dx = test_x[-c(1:100)]

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

mh.tn = function(o_ht,n_ht,x,b0,b1,sd,h1,h2,e_nu){
  u = runif(1)
  pn1 = dstd(x,sd = sqrt(exp(n_ht)),nu = e_nu)
  pn2 = dnorm(h2,mean =  b0 + b1*n_ht, sd = sd)
  po1 = dstd(x,sd = sqrt(exp(o_ht)),nu = e_nu)
  po2 = dnorm(h2,mean =  b0 + b1*o_ht, sd = sd)
  a = (pn1*pn2)/(po1*po2)
  a = min(1,a)
  if(a > u){
    res = n_ht
  } else{
    res = o_ht
  }
  return(res)
}
#initial guess
hs = log(drv^2) # sigma^2
# m1 = lm(hs[-1]~hs[-length(hs)])
# tb0 = m1$coefficients[1] # beta_0
# tb1 = m1$coefficients[2] # beta_1
# tresid = resid(m1)
# tetasd = sd(tresid)# delta
# tdq = fitdist(distribution = 'std',tresid)
# tdf = tdq$pars[3] # eta_df
ths = hs
registerDoMC(8)
est_paras = foreach(i = 1:200,.combine = rbind, .errorhandling = 'remove')%dopar%{
  set.seed(i)
  #ths1 = hs[251:1750]
  afspec = arfimaspec(mean.model = list(armaOrder = c(1, 0)),
                      distribution.model = 'norm')
  afm = arfimafit(afspec,hs)
  tb1 = afm@model$pars[2,1] # beta_1
  tb0 = (1-afm@model$pars[2,1])*afm@model$pars[1,1] # beta_0
  tetasd = afm@model$pars[7,1] # delta
  z = dx / sqrt(exp(hs))
  epn = fitdist(distribution = "std",z)$pars[3]

  intmr = data.frame(b0 = numeric(),b1 = numeric(),
                     sd = numeric(),epn = numeric())
  ths = hs
  for(l in 1:ite){
    n = length(ths)
    for(j in 2:(n-1)){
      o_ht = ths[j]
      n_ht = tb0 + tb1 * ths[j-1] + rnorm(1,mean = 0,sd = tetasd)
      # Normal-Std
      #t_ht = mh.ht(o_ht,n_ht,dx[j],
      #             tb0,tb1,tetasd,tdf,ths[j-1],ths[j+1])
      # Std-Std
      t_ht = mh.tn(o_ht,n_ht,dx[j],tb0,tb1,tetasd,
                    ths[j-1],ths[j+1],epn)
      ths[j] = t_ht
    }
    ths[1] = mean(ths)
    ths[length(ths)] = mean(ths)#tb0 + tb1 * ths[n-1] + rstd(1,mean = 0,sd = tetasd,nu = tdf)
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
    #tdf = afm@model$pars[17,1]
    tx = dx[251:1750]
    #################################################
    # Std-Std
    z = tx / sqrt(exp(ths1)) # back-out first innovation
    epn = fitdist(distribution = "std",z)$pars[3]
    intmr = rbind(intmr,c(tb0,tb1,tetasd,epn))
    ###################################################
    #if(l%%50 == 0){print(paste(l/50,"% finished"))}
    #intmr = rbind(intmr,c(tb0,tb1,tetasd,tdf))
  }
  intmr = intmr[-c(1:burn_in),]
  
  
  ## stochvol
  draws = svsample(dx, draws = 3000,priornu = c(2,100),burnin = burn_in)

  b1 = draws$summary$para[2,1]
  b0 = draws$summary$para[1,1] * (1-b1)
  svg = draws$summary$para[3,1]
  snu = draws$summary$para[4,1]
  # estimated parameters for sv(nn) model
  c(median(intmr[,1]),median(intmr[,2]),median(intmr[,3]),median(intmr[,4]),
    b0,b1,svg,snu)
}

Sys.time()-s
write.csv(est_paras,'sp_paras_tn.csv')

paras = paras[,-11]
w1 = paras[order(paras[,1]),1:3]
np1 = w1[c(-c(1:15),-c(112:126)),]
w2 = paras[order(paras[,4]),4:7]
np2 = w2[c(-c(1:15),-c(112:126)),]
