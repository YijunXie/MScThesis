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

set.seed(1)
s = Sys.time()
## conditional VaR
## Xt > VaR(Xt): not realistic settings
ite = 5000
N = 4500
burn_in = 2000
phi = 0.95
t_nu = 5
epi_nu = 5
omega =  -0.5
delta = 0.35
set.seed(1)
err1 = rnorm(N)
err2 = rstd(N,mean = 0,sd = delta,nu = t_nu) # eta_t
all_x = rep(0,N)
sig = rep(0,N)
sig[1] = 0.01
for(i in 2:N){
  sig[i] = sqrt(exp(omega+phi * log(sig[i-1]^2) + err2[i]))
  all_x[i] = sig[i]*err1[i]
}

test_x = all_x[2001:4500]

x = all_x[-c(1:burn_in)]

theta_hat = c(-0.502,  0.950,  0.357, 5.37, 4.98)

theta_g = c(2.791076e-06, 2.208389e-01, 7.781607e-01)

alpha = 0.99
w = 1000
v = quantile(test_x,alpha)  # (length(x)-w)
registerDoMC(8)
res = foreach(i = 1:(length(x)-w),.errorhandling = 'remove', .combine = rbind)%dopar%{
  #s1 = Sys.time()
  tempdat = x[i:(i+w-1)]
  #v = quantile(tempdat,alpha)
  ## garch
  spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE), 
                    distribution.model="std") # AR(1)-GARCH(1,1) model
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
  sn = matrix(nrow = 1e3,ncol = 2)
  while(n < 1e3){
    eta_1 = rstd(1,mean = 0,sd = theta_hat[3],nu = theta_hat[4])
    as_1 = sqrt(exp(theta_hat[1] + theta_hat[2] * log(s_0^2) + eta_1))
    ax_1 = as_1 * rstd(1, nu = theta_hat[5])
    if(ax_1 >= v){
      sn[n+1,] = c(as_1,eta_1)
      n = n+1
    }
  }
  eta_2 = rstd(1e3,mean = 0,sd = theta_hat[3],nu = theta_hat[4])
  as_2 = sqrt(exp(theta_hat[1] + theta_hat[2] * log(sn[,1]^2) + eta_2))
  x_a = as_2 * rstd(1e3, nu = theta_hat[5])
  covar_a = quantile(x_a,alpha)

  c(covar_g, covar_a ,i)

}
Sys.time()-s
ind = res[,3]+1000
plot(x[ind],type = 'l')#,ylim = c(-0.05,0.2))
lines(res[,1],type = 'l',col = 'red')
lines(res[,2],type = 'l',col = 'blue')
qqplot(res[,2],res[,1])
abline(0,1)
save(res,file = 'cvar_sv_99.Rda')



### GARCH
set.seed(1)
s = Sys.time()
## conditional VaR
## Xt > VaR(Xt): not realistic settings
ite = 5000
N = 4500
burn_in = 2000
alpha  = 0.85
beta = 0.1
t_nu = 5
omega =  5e-6
set.seed(1)
err = rstd(N,mean = 0,sd = 1,nu = t_nu) # eta_t
all_x = rep(0,N)
sig = rep(0,N)
sig[1] = 0.01
all_x[1] = 0.01
for(i in 2:N){
  sig[i] = sqrt(omega + alpha*sig[i-1]^2 + beta*all_x[i-1]^2)
  all_x[i] = sig[i]*err[i]
}

test_x = all_x[2001:4500]

# estimate parameters
price = rep(0,2501)
price[1] = 1
for(j in 2:2501){
  price[j] = price[j-1]*exp(test_x[j-1])
}

rv = TTR::volatility(price,calc="close", n = 5) / sqrt(260)
dx = test_x[251:2250]
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
intmr = intmr[-c(1:2000),]

# estimated parameters for sv(nt) model
theta_hat = c(median(intmr[,1]),median(intmr[,2]),
              median(intmr[,3]), median(intmr[,4]), median(intmr[,5]))

# parameter for GARCH
spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE), 
                  distribution.model="norm") # AR(1)-GARCH(1,1) model
fit = ugarchfit(spec, dx)

omega_g = fit@model$pars[7,1]
alpha_g = fit@model$pars[8,1]
beta_g = fit@model$pars[9,1]

theta_g = c(omega_g, alpha_g, beta_g)


x = all_x[-c(1:burn_in)]

alpha = 0.99
w = 1000
v = quantile(test_x,alpha)
registerDoMC(8)
res = foreach(i = 1:(length(x)-w),.errorhandling = 'remove', .combine = rbind)%dopar%{
  #s1 = Sys.time()
  tempdat = x[i:(i+w-1)]
  #v = quantile(tempdat,alpha)
  ## garch
  spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE), 
                    distribution.model="std") # AR(1)-GARCH(1,1) model
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
  sn = matrix(nrow = 1e3,ncol = 2)
  while(n < 1e3){
    eta_1 = rstd(1,mean = 0,sd = theta_hat[3],nu = theta_hat[4])
    as_1 = sqrt(exp(theta_hat[1] + theta_hat[2] * log(s_0^2) + eta_1))
    ax_1 = as_1 * rstd(1, theta_hat[5])
    if(ax_1 >= v){
      sn[n+1,] = c(as_1,eta_1)
      n = n+1
    }
  }
  eta_2 = rstd(1e3,mean = 0,sd = theta_hat[3],nu = theta_hat[4])
  as_2 = sqrt(exp(theta_hat[1] + theta_hat[2] * log(sn[,1]^2) + eta_2))
  x_a = as_2 * rstd(1e3, nu = theta_hat[5])
  covar_a = quantile(x_a,alpha)
  
  c(covar_g, covar_a ,i)
  
}


Sys.time()-s
ind = res[,3]+1000
plot(x[ind],type = 'l')#,ylim = c(-0.05,0.2))
lines(res[,1],type = 'l',col = 'red')
lines(res[,2],type = 'l',col = 'blue')
qqplot(res[,2],res[,1])
abline(0,1)
save(res,file = 'cvar_garch_99.Rda')

