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

sp <- read.csv("sp.csv", row.names=1)
dat = apply(sp,2,rev)
p = dat[,6]
all_x = diff(log(p))

theta_hat1 = c(-0.09979134,0.9888713,0.07071118)

theta_hat2 = c(-0.1146881,0.987974,0.0845055,2.840286)

x_test = all_x[2022:6000]
w = 1000 # window si
p1 = 0.9
p2 = 0.95
alpha1 = 0.95
alpha2 = 0.99
# setting up comparison 
registerDoMC(8)
res = foreach(i = 1:(length(x_test)-w),.errorhandling = 'remove',
              .combine = rbind)%dopar%{
  temp_dat = x_test[i:(i+w-1)]
  spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE), 
                    distribution.model="std") # AR(1)-GARCH(1,1) model
  fit = ugarchfit(spec, temp_dat)
  resid=residuals(fit) 
  sigma_garch=sigma(fit)
  z0=resid/sigma_garch
  frcst = ugarchforecast(fit,n.ahead=1)
  sigt1 = sigma(frcst) # expected sigma for next day
  
  # sv(nn) var 0.95
  sigma_t=as.numeric(sigma_garch[w])
  tmu1 = theta_hat1[1] + theta_hat1[2] * log(sigma_t^2)#log(testsig[i+w-1]^2)
  msig1 = theta_hat1[3]
  var_sv_a1 = varx_nn(tmu = tmu1, sig = msig1, 
                 lev = 1-alpha1, accuracy = 10e-15)
  
  # sv(nn) var 0.99
  var_sv_a2 = varx_nn(tmu = tmu1, sig = msig1, 
                      lev = 1-alpha2, accuracy = 10e-15)
  
  # sv(nt) var 0.95
  tmu2 = theta_hat2[1] + theta_hat2[2] * log(sigma_t^2)#log(testsig[i+w-1]^2)
  msig2 = theta_hat2[3]
  tnu = theta_hat2[4]
  var_sv_b1 = varx2(tmu = tmu2, sig = msig2, tnu = tnu,
                 lev = 1-alpha1, accuracy = 10e-15)
  
  # sv(nt) var 0.99
  var_sv_b2 = varx2(tmu = tmu2, sig = msig2, tnu = tnu,
                    lev = 1-alpha2, accuracy = 10e-15)
  
  #garch-evt var 0.95
  z0 = as.numeric(z0)
  u0 = sort(z0,decreasing=T)[p1*length(z0)+1]
  bo_gpd = gpd.fit(z0,threshold=u0,show=F)
  bo_beta=bo_gpd$mle[1]; bo_xi=bo_gpd$mle[2] # MLE's of scale and shape parameters
  bo_fq = u0 +(bo_beta/bo_xi) * (((1-alpha1)/p1)^(-bo_xi) -1) # upper quantile of GPD
  var_bo1 =sigt1*bo_fq
  
  #garch-evt var 0.99
  u0 = sort(z0,decreasing=T)[p2*length(z0)+1]
  bo_gpd = gpd.fit(z0,threshold=u0,show=F)
  bo_beta=bo_gpd$mle[1]; bo_xi=bo_gpd$mle[2] # MLE's of scale and shape parameters
  bo_fq = u0 +(bo_beta/bo_xi) * (((1-alpha2)/p2)^(-bo_xi) -1) # upper quantile of GPD
  var_bo2 =sigt1*bo_fq
  
  
  
  # true value
  true_x = x_test[i+w]
  
  c(true_x,var_sv_a1, var_sv_b1,var_bo1,var_sv_a2,var_sv_b2,var_bo2,i)
}

write.csv(res,"sp_var.csv")
