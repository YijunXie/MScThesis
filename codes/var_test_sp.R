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

dat <- read.csv("~/MScThesis/codes/sp.csv", row.names=1)
dat = apply(dat,2,rev)
OHLC = as.numeric(dat[,c(1,2,3,4)])
OHLC = matrix(OHLC,ncol = 4)
price = dat[,6]

ret = diff(log(price))
x = ret

x_sample = x[1:2500]
p_sample = price[1:2501]
rv = TTR::volatility(p_sample,calc="close", n = 5) / sqrt(260)

dx = x[251:2250]
drv = rv[252:2251]

ite = 5000
burn_in = 2000
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

x_test = x[2501:7000]
w = 1000 # window si
p = 0.9
alpha = 0.95
# setting up comparison
registerDoMC(8)
res = foreach(i = 1:(length(x_test)-w),.errorhandling = 'remove',              .combine = rbind)%dopar%{
  temp_dat = x_test[i:(i+w-1)]
  spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE), 
                    distribution.model="norm") # AR(1)-GARCH(1,1) model
  fit = ugarchfit(spec, temp_dat)
  resid=residuals(fit) 
  sigma_garch=sigma(fit)
  z0=resid/sigma_garch
  frcst = ugarchforecast(fit,n.ahead=1)
  sigt1 = sigma(frcst) # expected sigma for next day
  
  # sv var
  sigma_t=as.numeric(sigma_garch[w])
  tmu = theta_hat[1] + theta_hat[2] * log(sigma_t^2)#log(testsig[i+w-1]^2)
  msig = theta_hat[3]
  tnu = theta_hat[4]
  var_sv = varx2(tmu = tmu, sig = msig, tnu = tnu,
                 lev = 0.05, accuracy = 10e-15)
  
  # garch-evt var
  z0 = as.numeric(z0)
  u0 = sort(z0,decreasing=T)[p*length(z0)+1]
  bo_gpd = gpd.fit(z0,threshold=u0,show=F)
  bo_beta=bo_gpd$mle[1]; bo_xi=bo_gpd$mle[2] # MLE's of scale and shape parameters
  bo_fq = u0 +(bo_beta/bo_xi) * (((1-alpha)/p)^(-bo_xi) -1) # upper quantile of GPD 
  var_bo =mut1 + sigt1*bo_fq
  
  # garch-mb var
  var_mb = mut1 + sigt1 * qnorm(alpha)

  
  # true value
  true_x = x_test[i+w]
  
  
  bo_backtest = 0
  sv_backtest = 0
  mb_backtest = 0
  
  if(var_bo <= true_x){
    bo_backtest = 1 # if violate, indicate 1
  } 
  
  if(var_mb <= true_x){
    mb_backtest = 1 # if violate, indicate 1
  }
  
  if(var_sv <= true_x){
    sv_backtest = 1 # if violate, indicate 1
  }
  
  # scoring 
  s_bo = pwls(var_bo,true_x,alpha)
  s_sv = pwls(var_sv,true_x,alpha)
  s_mb = pwls(var_mb,true_x,alpha)
  
  return(c(true_x,var_sv, sv_backtest, s_sv, var_bo, bo_backtest, s_bo,
           var_mb, mb_backtest, s_mb))
}