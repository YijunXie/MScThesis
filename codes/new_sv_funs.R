chi.bar = function(dat, q = 0.95, n = 1, conf = 0.95){
  c1 = dat[-c((length(dat)-n+1):length(dat))]
  c2 = dat[-c(1:n)]
  temp_dat = cbind(c1,c2)
  z = rowMins(temp_dat)
  u = quantile(dat,q)
  
  nu = length(which(z > u))/length(dat)
  chibar = 2*log(1-q) / log(nu) -1
  chibarvar = (((4 * log(1 - q)^2)/(log(nu)^4 * nu^2)) * 
                  nu * (1 - nu))/length(dat)
  chibarlb = chibar - qnorm((1+ conf)/2)*sqrt(chibarvar)
  chibarub = chibar + qnorm((1+ conf)/2)*sqrt(chibarvar)
  return(c(chibarlb, chibar, chibarub))
}

mh.ht = function(o_ht,n_ht,x,b0,b1,sd,df,h1,h2){
  u = runif(1)
  pn1 = dnorm(x,sd = sqrt(exp(n_ht)))
  pn2 = dstd(h2,mean =  b0 + b1*n_ht, sd = sd,nu = df)
  po1 = dnorm(x,sd = sqrt(exp(o_ht)))
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

tempf = function(x,mu,sig,tnu){
  f = function(y){
    d1 = 1 - pstd(log((x/y)^2)-mu,mean = 0,sd = sig, nu = tnu)
    d2 = dnorm(y)
    return(d1 * d2)
  }
  return(f)
}

varx = function(tmu,sig,tnu,lev){
  x = seq(0,0.07,10e-6)
  ps = rep(0,length(x))
  for(k in 1:length(x)){
    ps[k] = integrate(tempf(x = x[k],mu = tmu, sig = sig, 
                            tnu = tnu),0,Inf)$value
  }
  ind = which.min(abs(ps - lev))
  if(ps[ind] > lev){ind = ind+1}
  return(x[ind])
}

varx2 = function(tmu,sig,tnu,lev,accuracy){
  left = 0; right = 1
  diff = right-left
  while(diff > accuracy){
    mid = (left + right)/2
    mp = integrate(tempf(x = mid,mu = tmu, sig = sig, 
                         tnu = tnu),0,Inf)$value
    if(mp > lev){
      left = mid
    } else{
      right = mid
    }
    diff = right - left
  }
  return(mid)
}

pwls = function(x,t,alpha){
  q = alpha
  if(x <= t){
    s = q*abs(x-t)
  } else {
    s = (1-q)*abs(x-t)
  }
  return(s)
}

lrtest = function(dat,p){
  y = sum(dat)
  n = length(dat)
  phat = y/n
  test_stat = 2*(y*log(phat/p) + (n-y)*log((1-phat)/(1-p)))
  z = test_stat
  p_value = 1 - pchisq(z,df = 1)
  return(list(z,p_value))
}

to.unitfrechet = function(dat){
  ind = rep(0,length(dat))
  cdfs = Ecdf(dat,pl = F)
  tempx = cdfs$x[-1]
  for(i in 1:length(dat)){
    ind[i] = which(dat == tempx[i])
  }
  fx = cdfs$y[ind]
  res = -1 / log(fx)
  return(res)
}

from.unitfrechet = function(ufdat,odat){
  cd = exp(-1 / ufdat)
  x = quantile(odat,cd)
  return(x)
}