require(ismev)
require(foreach)
require(doMC)
require(fGarch)
require(rugarch)
require(stochvol)
source('sv_funs.R')

s = Sys.time()
# simulate data
set.seed(2)
N = 8000
burnin = 1500
phi = 0.95
t_nu = 5
epi_nu = 5
omega = -0.5
err1 = rstd(N,mean = 0,sd = 1,nu = epi_nu)
err2 = rstd(N,mean = 0,sd = 0.35,nu = t_nu)

x = rep(0,N)
sig = rep(0,N)
sig[1] = 0.01
for(i in 2:N){
  sig[i] = sqrt(exp(omega+phi * log(sig[i-1]^2) + err2[i]))
  x[i] = sig[i]*err1[i]
}

x = x[-c(1:burnin)]
sig = sig[-c(1:burnin)]

price = rep(0,(N-burnin+1))
price[1] = 1
for(j in 2:(N-burnin+1)){
  price[j] = price[j-1]*exp(x[j-1])
}

testx = x[2501:6500]
testsig = sig[2501:6500]
trv = TTR::volatility(price,calc="close", n = 9) / sqrt(260)
testrv = trv[2501:6500]

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
                         tnu = tnu,enu = enu),0,Inf)$value
    if(mp > lev){
      left = mid
    } else{
      right = mid
    }
    diff = right - left
  }
  return(mid)
}

w = 1000
registerDoMC(10)
rest = foreach(i = 1:(length(testx)-w),.errorhandling = 'remove',
              .combine = rbind)%dopar%{
  tx = testx[i:(i+w-1)]
  t_spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE),
                      variance.model = list(garchOrder = c(1,1)),
                      distribution.model='std')
  t_fit = ugarchfit(t_spec, tx)
  sigma.t=as.numeric(sigma(t_fit)[w])
#   draws = svsample(tx, draws = 5000,burnin = 2000,
#                    priornu = c(2,10))
#   ht = draws$summary$latent[1000,1]
  #lev = 0.001 * which(cs<0.08)[2]
#   lev = 0.011
#   tdc = hill.est(testsig,lev)
  
  tmu = -0.9106 + 0.9084 * log(sigma.t^2)#log(testsig[i+w-1]^2)
  msig = 0.5547
  tnu = 4.58
  enu = 6.02
  va = varx2(tmu = tmu, sig = msig, tnu = tnu,
               enu = enu, lev = 0.01, accuracy = 10e-15)
  
  tmu_o = -0.5 + 0.95 * log(sigma.t^2) #log(testsig[i+w-1]^2)
  msig_o = 0.35
  tnu_o = 5
  enu_o = 5
  va_o = varx2(tmu = tmu_o, sig = msig_o, tnu = tnu_o,
            enu = enu_o, lev = 0.01, accuracy = 10e-15)
#   c = integrate(tempf(x = 0.00745,mu = tmu, sig = msig, 
#                       tnu = tnu,enu = enu),0,Inf)$value
#   nal = 1 - (1 - 0.95) / c
#   mu = -0.5 + 0.95 * log(testsig[i+w]^2)
#   ls = qstd(nal,mean = mu, sd = 0.35,nu = 5)
#   va = sqrt(exp(ls))
#   k = floor(0.12*length(tsig)+1)
#   u = unname(sort(tsig,decreasing = T)[k])
#   GPDfit = gpd.fit(tsig,threshold=u,show=F)
#   beta=GPDfit$mle[1]; xi=GPDfit$mle[2]
#   qPOT = u + (beta/xi) * (((1-nal)/0.12)^(-xi) -1)
  c(va,va_o,i)
}
Sys.time()-s
plot(testx[1001:4000],type = 'l')
lines(rest[,1],type = 'l',col = 'red')
lines(rest[,2],type = 'l',col = 'blue')
length(which(testx[1001:4000] > rest[,1])) / nrow(rest)
length(which(testx[1001:4000] > rest[,2])) / nrow(rest)
legend('topright',c('estimated','optimal'),
       col = c('red','blue'),lwd = c(2,2))
write.csv(rest,'var_bisec_99_2.csv')

dat = data.frame(estimated = rest[,1],
                 optimal = rest[,2],
                 `return` = testx[1001:4000],
                 time = rest[,3])

dat %>% ggplot(mapping = aes(time))+ 
  geom_line(aes(y = return, color = "return"),data = dat)+
  geom_line(aes(y = estimated, color = "estimated"),data = dat,  alpha=0.4) + 
  geom_line(aes(y = optimal, color = "optimal"),data = dat,  alpha=0.2)

