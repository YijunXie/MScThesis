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
ite = 6000
N = 8000
burn_in = 5000
phi = 0.3#0.95
t_nu = 3
epi_nu = 5
omega = -5 # -0.5
delta = 1 # 0.35

err1 = rnorm(N)
err2 = rstd(N,mean = 0,sd = delta,nu = t_nu) # eta_t
x = rep(0,N)
sig = rep(0,N)
sig[1] = 0.01
for(i in 2:N){
  sig[i] = sqrt(exp(omega+phi * log(sig[i-1]^2) + err2[i]))
  x[i] = sig[i]*err1[i]
}

test_x = x[3001:5000]
# estimate parameters
price = rep(0,2001)
price[1] = 1
for(j in 2:2001){
  price[j] = price[j-1]*exp(test_x[j-1])
}

rv = TTR::volatility(price,calc="close", n = 5) / sqrt(260)
drv = rv[-c(1:100)]
dx = test_x[-c(1:100)]
#initial guess
hs = log(drv^2) # sigma^2
m1 = lm(hs[-1]~hs[-length(hs)])
tb0 = m1$coefficients[1] # beta_0
tb1 = m1$coefficients[2] # beta_1
tresid = resid(m1)
tetasd = sd(tresid)# delta
tdq = fitdist(distribution = 'std',tresid)
tdf = tdq$pars[3] # eta_df
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
  ths[1] = mean(ths)
  ths[length(ths)] = tb0 + tb1 * ths[n-1] + rstd(1,mean = 0,sd = tetasd,nu = tdf)
  ths1 = ths[251:1750] # truncated log(sigma^2)
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
intmr = intmr[-c(1:2000),]

# estimated parameters for sv(nt) model
theta_hat = c(median(intmr[,1]),median(intmr[,2]),
              median(intmr[,3]), median(intmr[,4]))

# estimated parameters for sv(nn) model
draws = svsample(dx, draws = 3000,burnin = 2000)

b1 = draws$summary$para[2,1]
b0 = draws$summary$para[1,1] * (1-b1)
svg = draws$summary$para[3,1]

theta_hat_n = c(b0,b1,svg)

Sys.time()-s
x = x[-c(1:burn_in)]



w = 1000
registerDoMC(8)
res = foreach(i = 1:(length(x)-w),.errorhandling = 'remove', .combine = rbind)%dopar%{
  tempdat = x[i:(i+w-1)]
  v = quantile(tempdat,0.95)
  ## garch
  spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE), 
                    distribution.model="norm") # AR(1)-GARCH(1,1) model
  fit = ugarchfit(spec, tempdat)
  s_0 = fit@fit$sigma[w]
  omega_g = fit@model$pars[7,1]
  alpha_g = fit@model$pars[8,1]
  beta_g = fit@model$pars[9,1]
  gs_1 = sqrt(omega_g + beta_g*s_0^2 + alpha_g*tempdat[w]^2)
  cdfs = runif(1e6,0.95,1)
  ep_1 = qnorm(cdfs)
  gs_2 = sqrt(omega_g + beta_g*gs_1^2 + (alpha_g*gs_1^2)*ep_1^2)
  x_g = rnorm(1e6) * gs_2
  covar_g = quantile(x_g,0.95)
  
  ## arsv_estimated
  n = 0
  sn = matrix(nrow = 1e3,ncol = 2)
  while(n < 1e3){
    eta_1 = rstd(1,mean = 0,sd = theta_hat[3],nu = theta_hat[4])
    as_1 = sqrt(exp(theta_hat[1] + theta_hat[2] * log(s_0^2) + eta_1))
    ax_1 = as_1 * rnorm(1)
    if(ax_1 >= v){
      sn[n+1,] = c(as_1,eta_1)
      n = n+1
    }
  }
  eta_2 = rstd(1e3,mean = 0,sd = theta_hat[3],nu = theta_hat[4])
  as_2 = sqrt(exp(theta_hat[1] + theta_hat[2] * log(sn[,1]^2) + eta_2))
  x_a = as_2 * rnorm(1e3)
  covar_a = quantile(x_a,0.95)
  
  ## arsv_oracle
  n = 0
  os = matrix(nrow = 1e3,ncol = 2)
  while(n < 1e3){
    eta_1 = rstd(1,mean = 0,sd = delta,nu = t_nu)
    os_1 = sqrt(exp(-5 + 0.3 * log(s_0^2) + eta_1))
    ox_1 = os_1 * rnorm(1)
    if(ox_1 >= v){
      os[n+1,] = c(os_1,eta_1)
      n = n+1
    }
  }
  oeta_2 = rstd(1e3,mean = 0,sd = delta,nu = t_nu)
  os_2 = sqrt(exp(-5 + 0.3 * log(os[,1]^2) + oeta_2))
  x_o = os_2 * rnorm(1e3)
  covar_o = quantile(x_o,0.95)
  
  ## arsv_normal
  n = 0
  ns = matrix(nrow = 1e3,ncol = 2)
  while(n < 1e3){
    eta_1 = rnorm(1,sd = theta_hat_n[3])
    ns_1 = sqrt(exp(theta_hat_n[1] + theta_hat_n[2] * log(s_0^2) + eta_1))
    nx_1 = ns_1 * rnorm(1)
    if(nx_1 >= v){
      ns[n+1,] = c(ns_1,eta_1)
      n = n+1
    }
  }
  neta_2 = rnorm(1e3,sd = theta_hat_n[3])
  ns_2 = sqrt(exp(theta_hat_n[1] + theta_hat_n[2] * log(ns[,1]^2) + neta_2))
  x_n = ns_2 * rnorm(1e3)
  covar_n = quantile(x_n,0.95)
  
  c(covar_g, covar_a, covar_o, covar_n ,i)

}
Sys.time()-s
ind = res[,5]+1000
plot(x[ind],type = 'l')
lines(res[,1],type = 'l',col = 'red')
lines(res[,2],type = 'l',col = 'blue')
lines(res[,3],type = 'l',col = 'green')

s_g = s_a = s_o = rep(0,nrow(res))
for(j in 1:nrow(res)){
  s_g[j] = pwls(res[j,1],x[res[j,4]+1000],0.95)
  s_a[j] = pwls(res[j,2],x[res[j,4]+1000],0.95)
  s_o[j] = pwls(res[j,3],x[res[j,4]+1000],0.95)
}
efp.test(s_g,s_a)
efp.test(s_g,s_o)
efp.test(s_o,s_a)
mean(s_g)
mean(s_a)
mean(s_o)
