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
set.seed(1)
s = Sys.time()
## conditional VaR
## Xt > VaR(Xt): not realistic settings
ite = 5000
N = 9000
burn_in = 5000
t_nu = 5
phi = 0.95
omega = -0.5


err1 = rnorm(N)
err2 = rstd(N,mean = 0,sd = 1,nu = t_nu) # eta_t

x = rep(0,N)
sig = rep(0,N)
sig[1] = 0.01
x[1] = 0.01
for(i in 2:N){
  sig[i] = sqrt(exp(omega+phi * log(sig[i-1]^2) + err2[i]))
  x[i] = sig[i]*err1[i]
}

x = x[-c(1:3000)]
sig = sig[-c(1:3000)]

theta_hat = c( -0.506,   0.950,  0.357,  6.54)
#theta_hat1 = c(-0.5918670,  0.9393397,  0.3703698)

x_test = x[2501:6000]
sig_test = sig[2501:6000]
w = 1000 # window si
p1 = 0.85
p2 = 0.95
alpha1 = 0.9
alpha2 = 0.99
# setting up comparison
registerDoMC(8)
res = foreach(i = 1:(length(x_test)-w),.errorhandling = 'remove',   
              .combine = rbind)%dopar%{
  temp_dat = x_test[i:(i+w-1)]
  #temp_dat = x_test[i:(i+w-1)]
  spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE), 
                    distribution.model="std") # AR(1)-GARCH(1,1) model
  fit = ugarchfit(spec, temp_dat)
  resid=residuals(fit) 
  sigma_garch=sigma(fit)
  z0=resid/sigma_garch
  frcst = ugarchforecast(fit,n.ahead=1)
  sigt1 = sigma(frcst) # expected sigma for next day
  
  # arsv(1) var 0.95
  sigma_t=as.numeric(sigma_garch[w])
  tmu1 = theta_hat[1] + theta_hat[2] * log(sigma_t^2)#log(testsig[i+w-1]^2)
  msig1 = theta_hat[3]
  tnu = theta_hat[4]
  var_sv_a = varx2(tmu = tmu1, sig = msig1, tnu = tnu,
                    lev = 1-alpha1, accuracy = 10e-15)
  
  # arsv(1) var 0.99
  var_sv_b = varx2(tmu = tmu1, sig = msig1, tnu = tnu,
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
  
  arsv_s_1 = pwls(true_x,var_sv_a,1-alpha1)
  arsv_s_2 = pwls(true_x,var_sv_b,1-alpha2)
  
  garch_s_1 = pwls(true_x,var_bo1,1-alpha1)
  garch_s_2 = pwls(true_x,var_bo2,1-alpha2)
  
  c(true_x,var_sv_a, arsv_s_1, var_sv_b,arsv_s_2, var_bo1,garch_s_1,
    var_bo2,garch_s_2,i)
}

Sys.time()-s
write.csv(res,"sim_var_90.csv")
