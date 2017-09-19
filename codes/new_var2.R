setwd("/nfs/37zfs1-yijun.xie/yijun.xie/MScThesis/codes")
require(ismev)
require(foreach)
require(doMC)
require(fGarch)
require(rugarch)
require(stochvol)
source('sv_funs.R')

s = Sys.time()
# simulate ARSV
set.seed(1)
s = Sys.time()
## conditional VaR
## Xt > VaR(Xt): not realistic settings
ite = 5000
N = 6500
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

test_x = all_x[2001:4000]

x = all_x[4001:6500]

# price = rep(0,2001)
# price[1] = 1
# for(j in 2:(N-burnin+1)){
#   price[j] = price[j-1]*exp(x[j-1])
# }

# testx = x
# testsig = sig
# trv = TTR::volatility(price,calc="close", n = 9) / sqrt(260)
# testrv = trv[2001:4500]


tempf = function(x,mu,sig,tnu,enu){
  f = function(y){
    d1 = 1 - pstd(log((x/y)^2)-mu,mean = 0,sd = sig, nu = tnu)
    d2 = dstd(y,nu = enu)
    return(d1 * d2)
  }
  return(f)
}

varx = function(tmu,sig,tnu,enu,lev){
  x = seq(0,0.07,10e-6)
  ps = rep(0,length(x))
  for(k in 1:length(x)){
    ps[k] = integrate(tempf(x = x[k],mu = tmu, sig = sig, 
                            tnu = tnu,enu = enu),0,Inf)$value
  }
  ind = which.min(abs(ps - lev))
  if(ps[ind] > lev){ind = ind+1}
  return(x[ind])
}

varx2 = function(tmu,sig,tnu,enu,lev,accuracy){
  left = 0; right = 1
  diff = right-left
  while(diff > accuracy){
    mid = (left + right)/2
    mp = integrate(tempf(x = mid,mu = tmu, sig = sig, 
                         tnu = tnu, enu = enu),0,Inf)$value
    if(mp > lev){
      left = mid
    } else{
      right = mid
    }
    diff = right - left
  }
  return(mid)
}

spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE), 
                  distribution.model="norm") # AR(1)-GARCH(1,1) model
fit = ugarchfit(spec, test_x)

omega_g = fit@model$pars[7,1]
alpha_g = fit@model$pars[8,1]
beta_g = fit@model$pars[9,1]

theta_g = c(omega_g, alpha_g, beta_g)

lva = 0.01
w = 1000
registerDoMC(8)
rest = foreach(i = 1:(length(x)-w),.errorhandling = 'remove', 
               .combine = rbind)%dopar%{
  tx = x[i:(i+w-1)]
  #out_sig = testrv[i+w-1]
  t_spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE),
                      variance.model = list(garchOrder = c(1,1)),
                      distribution.model='norm')
  t_fit = ugarchfit(t_spec, tx)
  sigma.t=as.numeric(sigma(t_fit)[w])
  ret.t=as.numeric(residuals(t_fit)[w])
  sig.t2 = sqrt(theta_g[1] + theta_g[2] * ret.t^2 + 
                  theta_g[3] *sigma.t^2)
  var_garch = sig.t2 * qnorm(1-lva)

  
  tmu = -0.561 + 0.946 * log(sigma.t^2)#log(testsig[i+w-1]^2)
  msig = 0.363
  tnu = 5.37
  enu = 5.44
  var_arsv = varx2(tmu = tmu, sig = msig, tnu = tnu, enu = enu,
              lev = lva, accuracy = 10e-15)
  
  tmu_o = -0.5 + 0.95 * log(sigma.t^2) #log(testsig[i+w-1]^2)
  msig_o = 0.35
  tnu_o = 5
  enu_o = 5
  var_o = varx2(tmu = tmu_o, sig = msig_o, tnu = tnu_o,enu = enu_o,
             lev = lva, accuracy = 10e-15)
  c(var_garch,var_arsv,var_o,i)
}
Sys.time()-s
write.csv(rest,'var_arsv_99.csv')
ind = rest[,4]+1000
w1 = (x[ind] > rest[,1])
w2 = (x[ind] > rest[,2])
LRuc(0.01,w1)
sum(w1)/length(w1)
LRuc(0.01,w2)
sum(w2)/length(w2)
ind = 1001:2500
plot(x[ind],type = 'l', 
     ylim = c(-0.05,0.1),
     ylab = expression('X'[t]))
lines(rest[,1],type = 'l',col = 'red')
lines(rest[,2],type = 'l',col = 'blue')
legend("topright", c("VaR under GARCH","VaR under ARSV"), lty = c(1,1),
       col = c("red","blue"))
### GARCH
set.seed(1)
s = Sys.time()
ite = 5000
N = 7000
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
x = all_x[4501:7000]

# # estimate parameters
# price = rep(0,2501)
# price[1] = 1
# for(j in 2:2501){
#   price[j] = price[j-1]*exp(test_x[j-1])
# }
# 
# rv = TTR::volatility(price,calc="close", n = 5) / sqrt(260)
# dx = test_x[251:2250]
# dsig = sig[251:2250]
# drv = rv[252:2251]
# #initial guess
# hs = log(drv^2)
# m1 = lm(hs[-1]~hs[-length(hs)])
# tb0 = m1$coefficients[1]
# tb1 = m1$coefficients[2]
# tresid = resid(m1)
# tetasd = sd(tresid)#fitdist(distribution = 'std',tresid)$pars[2]
# tdq = fitdist(distribution = 'std',tresid)
# tdf = tdq$pars[3]
# z = dx / sqrt(exp(hs))
# epn = fitdist(distribution = "std",z)$pars[3]
# ths = hs
# intmr = data.frame(b0 = numeric(),b1 = numeric(),
#                    sd = numeric(),df = numeric(),
#                    epn = numeric())
# for(l in 1:ite){
#   n = length(ths)
#   for(j in 2:(n-1)){
#     o_ht = ths[j]
#     n_ht = tb0 + tb1 * ths[j-1] + rstd(1,mean = 0,sd = tetasd,nu = tdf)
#     # Normal-Std
#     #t_ht = mh.ht(o_ht,n_ht,dx[j],
#     #             tb0,tb1,tetasd,tdf,ths[j-1],ths[j+1])
#     # Std-Std
#     t_ht = mh.std(o_ht,n_ht,dx[j],tb0,tb1,tetasd,
#                   tdf,ths[j-1],ths[j+1],epn)
#     ths[j] = t_ht
#   }
#   ths[1] = mean(ths)
#   ths[length(ths)] = tb0 + tb1 * ths[n-1] + rstd(1,mean = 0,sd = tetasd,nu = tdf)
#   ths1 = ths[251:1750]
#   afspec = arfimaspec(mean.model = list(armaOrder = c(1, 0)),
#                       distribution.model = 'std')
#   afm = arfimafit(afspec,ths1)
#   tb1 = afm@model$pars[2,1]
#   tb0 = (1-afm@model$pars[2,1])*afm@model$pars[1,1]
#   tetasd = afm@model$pars[7,1]
#   tdf = afm@model$pars[17,1]
#   tx = dx[251:1750]
#   #################################################
#   # Std-Std
#   z = tx / sqrt(exp(ths1)) # back-out first innovation
#   epn = fitdist(distribution = "std",z)$pars[3]
#   intmr = rbind(intmr,c(tb0,tb1,tetasd,tdf,epn))
#   ###################################################
# }
# intmr = intmr[-c(1:2000),]
# 
# # estimated parameters for sv(nt) model
# theta_hat = c(median(intmr[,1]),median(intmr[,2]),
#               median(intmr[,3]), median(intmr[,4]), median(intmr[,5]))

theta_hat = c(-0.4304139, 0.9551008,  0.1778457,
              14.2585930,10.1077356)
# parameter for GARCH
spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE), 
                  distribution.model="norm") # AR(1)-GARCH(1,1) model
fit = ugarchfit(spec, dx)

omega_g = fit@model$pars[7,1]
alpha_g = fit@model$pars[8,1]
beta_g = fit@model$pars[9,1]

theta_g = c(omega_g, alpha_g, beta_g)


#x = all_x[-c(1:burn_in)]

rest = foreach(i = 1:(length(x)-w),.errorhandling = 'remove',
               .combine = rbind)%dopar%{
  tx = x[i:(i+w-1)]
  #out_sig = testrv[i+w-1]
  t_spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE),
                      variance.model = list(garchOrder = c(1,1)),
                      distribution.model='norm')
  t_fit = ugarchfit(t_spec, tx)
  sigma.t=as.numeric(sigma(t_fit)[w])
  ret.t=as.numeric(residuals(t_fit)[w])
  sig.t2 = sqrt(theta_g[1] + theta_g[2] * ret.t^2 + 
                  theta_g[3] *sigma.t^2)
  var_garch = sig.t2 * qnorm(1-lva)
  
  
  tmu = theta_hat[1] + theta_hat[2] * log(sigma.t^2)#log(testsig[i+w-1]^2)
  msig = theta_hat[3]
  tnu = theta_hat[4]
  enu = theta_hat[5]
  var_arsv = varx2(tmu = tmu, sig = msig, tnu = tnu, enu = enu,
                   lev = lva, accuracy = 10e-15)
  
  # tmu_o = -0.5 + 0.95 * log(sigma.t^2) #log(testsig[i+w-1]^2)
  # msig_o = 0.35
  # tnu_o = 5
  # enu_o = 5
  # var_o = varx2(tmu = tmu_o, sig = msig_o, tnu = tnu_o,enu = enu_o,
  #               lev = lva, accuracy = 10e-15)
  c(var_garch,var_arsv,i)
}
Sys.time()-s
write.csv(rest,'var_garch_99.csv')

# ind = as.vector(var_garch_95[,4])+1000
# ind = as.numeric(as.matrix(ind))
ind = rest[,3]+1000
w1 = (x[ind] > rest[,1])
sum(w1)/length(w1)
w2 = (x[ind] > rest[,2])
sum(w2)/length(w2)
LRuc(0.01,w1)
LRuc(0.01,w2)
plot(x[ind],type = 'l', 
     ylim = c(-0.05,0.09),
     ylab = expression('X'[t]))
lines(rest[,1],type = 'l',col = 'red')
lines(rest[,2],type = 'l',col = 'blue')
legend("topright", c("VaR under GARCH","VaR under ARSV"), lty = c(1,1),
       col = c("red","blue"))
