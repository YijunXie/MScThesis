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
N = 8000
burn_in = 5000
phi = 0.95
t_nu = 5
epi_nu = 5
omega =  -0.5
delta = 0.35

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

draws = svsample(dx, draws = 3000,burnin = 2000)

b1 = draws$summary$para[2,1]
b0 = draws$summary$para[1,1] * (1-b1)
svg = draws$summary$para[3,1]

theta_hat_n = c(b0,b1,svg)

# parameters for GARCH 
spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE), 
                  distribution.model="std") # AR(1)-GARCH(1,1) model
fit = ugarchfit(spec, dx)

omega_g = fit@model$pars[7,1]
alpha_g = fit@model$pars[8,1]
beta_g = fit@model$pars[9,1]

theta_g = c(omega_g, alpha_g, beta_g)

Sys.time()-s
x = x[-c(1:burn_in)]

#theta_hat = c(-0.5980024,  0.9395961,  0.3793616, 10.3933544)
#theta_hat_n = c(-0.5918670,  0.9393397,  0.3703698)
#theta_g = c(2.791076e-06, 2.208389e-01, 7.781607e-01)

alpha = 0.95
w = 1000
v = quantile(test_x,alpha)
registerDoMC(8)
res = foreach(i = 1:(length(x)-w),.errorhandling = 'remove', .combine = rbind)%dopar%{
  #s1 = Sys.time()
  tempdat = x[i:(i+w-1)]
  #v = quantile(tempdat,alpha)
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
  covar_a = quantile(x_a,alpha)
  
  # ## arsv_oracle
  # n = 0
  # os = matrix(nrow = 1e3,ncol = 2)
  # while(n < 1e3){
  #   eta_1 = rstd(1,mean = 0,sd = delta,nu = t_nu)
  #   os_1 = sqrt(exp(omega + phi * log(s_0^2) + eta_1))
  #   ox_1 = os_1 * rnorm(1)
  #   if(ox_1 >= v){
  #     os[n+1,] = c(os_1,eta_1)
  #     n = n+1
  #   }
  # }
  # oeta_2 = rstd(1e3,mean = 0,sd = delta,nu = t_nu)
  # os_2 = sqrt(exp(omega + phi * log(os[,1]^2) + oeta_2))
  # x_o = os_2 * rnorm(1e3)
  # covar_o = quantile(x_o,alpha)
  # 
  # ## arsv_normal
  # n = 0
  # ns = matrix(nrow = 1e3,ncol = 2)
  # while(n < 1e3){
  #   eta_1 = rnorm(1,sd = theta_hat_n[3])
  #   ns_1 = sqrt(exp(theta_hat_n[1] + theta_hat_n[2] * log(s_0^2) + eta_1))
  #   nx_1 = ns_1 * rnorm(1)
  #   if(nx_1 >= v){
  #     ns[n+1,] = c(ns_1,eta_1)
  #     n = n+1
  #   }
  # }
  # neta_2 = rnorm(1e3,sd = theta_hat_n[3])
  # ns_2 = sqrt(exp(theta_hat_n[1] + theta_hat_n[2] * log(ns[,1]^2) + neta_2))
  # x_n = ns_2 * rnorm(1e3)
  # covar_n = quantile(x_n,alpha)
  c(covar_g, covar_a ,i)
  #c(covar_g, covar_a, covar_o, covar_n ,i)
  #Sys.time()-s1
}
Sys.time()-s
ind = res[,3]+1000
plot(x[ind],type = 'l')#,ylim = c(-0.05,0.2))
lines(res[,1],type = 'l',col = 'red')
lines(res[,2],type = 'l',col = 'blue')
qqplot(res[,2],res[,1])
abline(0,1)
