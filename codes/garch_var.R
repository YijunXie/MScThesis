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
s = Sys.time()
## conditional VaR
## Xt > VaR(Xt): not realistic settings
ite = 5000
N = 9000
burn_in = 5000
alpha  = 0.85
beta = 0.1
omega =  5e-6

#err = rstd(N,mean = 0,sd = 1,nu = t_nu) # eta_t
err = rnorm(N)
x = rep(0,N)
sig = rep(0,N)
sig[1] = 0.01
x[1] = 0.01
for(i in 2:N){
  sig[i] = sqrt(omega + alpha*sig[i-1]^2 + beta*x[i-1]^2)
  x[i] = sig[i]*err[i]
}

x = x[-c(1:3000)]
sig = sig[-c(1:3000)]

theta_hat = c(-0.5980024,  0.9395961,  0.3793616, 10.3933544)

x_test = x[2501:6000]
sig_test = sig[2501:6000]
w = 1000 # window si
p1 = 0.9
p2 = 0.95
alpha1 = 0.95
alpha2 = 0.99
# setting up comparison
registerDoMC(8)
res = foreach(i = 1:(length(x_test)-w),.errorhandling = 'remove',                 .combine = rbind)%dopar%{
  temp_dat = x_test[i:(i+w-1)]
  #temp_dat = x_test[i:(i+w-1)]
  spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE), 
                    distribution.model="norm") # GARCH(1,1) model
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
  nob = length(temp_dat)*(1-alpha1)
  z0 = as.numeric(z0)
  #u0 = sort(z0,decreasing=F)[p1*length(z0)+1]
  u0 = quantile(z0,p1)
  bo_gpd = gpd.fit(z0,threshold=u0,show=F,method = "BFGS")
  bo_beta=bo_gpd$mle[1]; bo_xi=bo_gpd$mle[2] # MLE's of scale and shape parameters
  bo_fq = u0 +(bo_beta/bo_xi) * ((nob*(1-p1))^(-bo_xi) -1) # upper quantile of GPD
  #var_bo1 =sigt1*bo_fq
  var_bo1 = quantile(sample(z0,1e5,replace = T),alpha1)*sigt1
  
  #garch-evt var 0.99
  #u0 = sort(z0,decreasing=F)[p2*length(z0)+1]
  u0 = quantile(z0,p2)
  nob = length(temp_dat)*(1-alpha2)
  bo_gpd = gpd.fit(z0,threshold=u0,show=F,method = "BFGS")
  bo_beta=bo_gpd$mle[1]; bo_xi=bo_gpd$mle[2] # MLE's of scale and shape parameters
  bo_fq = u0 +(bo_beta/bo_xi) * ((nob*(1-p2))^(-bo_xi) -1) # upper quantile of GPD
  #var_bo2 =sigt1*bo_fq
  var_bo2 = quantile(sample(z0,1e5,replace = T),alpha2)*sigt1
  
  # true value
  true_x = x_test[i+w]
  
  arsv_s_1 = pwls(var_sv_a,true_x,alpha1)
  arsv_s_2 = pwls(var_sv_b,true_x,alpha2)
  
  garch_s_1 = pwls(var_bo1,true_x,alpha1)
  garch_s_2 = pwls(var_bo2,true_x,alpha2)
  
  c(true_x,var_sv_a, arsv_s_1, var_sv_b,arsv_s_2, var_bo1,garch_s_1,
    var_bo2,garch_s_2,i)
}

res<- read.csv("~/MScThesis/codes/sim_var_garch.csv")
res = (res[,-1])
ind = res[,10]+1000
plot(x_test[ind],type = 'l', 
     ylim = c(-0.05,0.07),
     ylab = expression('X'[t]))
lines(res[,6],type = 'l',col = 'red')
lines(res[,2],type = 'l',col = 'blue')
legend("topright", c("VaR under GARCH","VaR under ARSV"), lty = c(1,1),
       col = c("red","blue"))

plot(x_test[ind],type = 'l', 
     ylim = c(-0.05,0.08),
     ylab = expression('X'[t]))
lines(res[,8],type = 'l',col = 'red')
lines(res[,4],type = 'l',col = 'blue')
legend("topright", c("VaR under GARCH","VaR under ARSV"), lty = c(1,1),
       col = c("red","blue"))

length(which(res[,2] < res[,1])) / nrow(res)
dat1 = (res[,2] < res[,1])
lrtest(dat1,0.05)
length(which(res[,4] < res[,1])) / nrow(res)
dat2 = (res[,4] < res[,1])
lrtest(dat2,0.01)

length(which(res[,6] < res[,1])) / nrow(res)
dat3 = (res[,6] < res[,1])
lrtest(dat3,0.05)
length(which(res[,8] < res[,1])) / nrow(res)
dat4 = (res[,8] < res[,1])
lrtest(dat4,0.01)

ecpa.test(res[,3],res[,7])
ecpa.test(res[,5],res[,9])

Sys.time()-s
write.csv(res,"sim_var_garch.csv")


