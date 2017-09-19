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
all_x = -diff(log(p))
test_x = all_x[1:2021]


#theta_hat = c(-0.1146881,0.987974,0.0845055,2.840286)
theta_hat = c(-0.1717536,0.98212,0.1152672,2.49915,10.46629)

spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE), 
                  distribution.model="norm") # AR(1)-GARCH(1,1) model
fit = ugarchfit(spec, test_x)

omega_g = fit@model$pars[7,1]
alpha_g = fit@model$pars[8,1]
beta_g = fit@model$pars[9,1]
#nu_g = fit@model$pars[17,1]

theta_g = c(omega_g, alpha_g, beta_g)

x = all_x[2023:7000]

alpha = 0.99
w = 1000
v = quantile(test_x,alpha)
s = Sys.time()
# 
registerDoMC(8)
res = foreach(i = 1:(length(x)-w),.errorhandling = 'remove', .combine = rbind)%dopar%{
  set.seed(i)
  tempdat = x[i:(i+w-1)]
  ## garch
  spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE), 
                    distribution.model="norm") # AR(1)-GARCH(1,1) model
  fit = ugarchfit(spec, tempdat)
  s_0 = fit@fit$sigma[w]
  gs_1 = sqrt(theta_g[1] + theta_g[3]*s_0^2 + theta_g[2]*tempdat[w]^2)
  cdfs = runif(1e6,alpha,1)
  ep_1 = qnorm(cdfs)
  gs_2 = sqrt(theta_g[1] + theta_g[3]*gs_1^2 + (theta_g[2]*gs_1^2)*ep_1^2)
  x_g = rnorm(1e6) * gs_2
  covar_g = quantile(x_g,alpha)
  
  ## arsv_estimated
  n = 0
  counter = 0
  sn = matrix(nrow = 1e3,ncol = 3)
  while(n < 1e3){
    counter = counter+1
    set.seed(counter)
    eta_a = rstd(1,mean = 0,sd = theta_hat[3],nu = theta_hat[4])
    eta_b = rstd(1,mean = 0,sd = theta_hat[3],nu = theta_hat[4])
    eta_1 = min(eta_a, eta_b)
    as_1 = sqrt(exp(theta_hat[1] + theta_hat[2] * log(s_0^2) + eta_1))
    ax_1 = as_1 * rstd(1,nu = theta_hat[5])
    if(ax_1 >= v){
      sn[n+1,] = c(as_1,eta_1,counter)
      n = n+1
    }
  }
  eta_2 = rstd(1e3,mean = 0,sd = theta_hat[3],nu = theta_hat[4])
  as_2 = sqrt(exp(theta_hat[1] + theta_hat[2] * log(sn[,1]^2) + eta_2))
  x_a = as_2 * rstd(1e3,nu = theta_hat[5])
  covar_a = quantile(x_a,alpha)

  c(covar_g, covar_a ,i)

}
Sys.time()-s
write.csv(res,"new_covar_sp_99.csv")
ind = res[,3]+1000
plot(x[ind],type = 'l',ylim = c(-0.06,0.08),
     ylab = expression('X'[t]))#,ylim = c(-0.05,0.2))
lines(res[,1],type = 'l',col = 'red')
lines(res[,2],type = 'l',col = 'blue')
legend("topleft", c("CoVaR under GARCH","CoVaR under ARSV"), lty = c(1,1),
       col = c("red","blue"))
diff = res[,1] - res[,2]
qtd = quantile(diff,0.95)
ins = which(diff>qtd)
nx = x[ind]^2
plot(nx,type = 'l', 
     ylab = expression('X'[t]^2))
for(i in 1:length(ins)){
  rect((ins[i]-2),0,(ins[i]),.1,border = rgb(0,0,1,0.2), col = rgb(0,0,1,0.2))
}
abline(h = quantile(nx,0.95))

