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

mh.dn = function(o_ht,n_ht,x,b0,b1,sd,h1,h2){
  u = runif(1)
  pn1 = dnorm(x,sd = sqrt(exp(n_ht)))
  pn2 = dnorm(h2,mean =  b0 + b1*n_ht, sd = sd)
  po1 = dnorm(x,sd = sqrt(exp(o_ht)))
  po2 = dnorm(h2,mean =  b0 + b1*o_ht, sd = sd)
  a = (pn1*pn2)/(po1*po2)
  if(is.na(a>0)){
    a = 0
  }
  a = min(1,a)
  if(a > u){
    res = n_ht
  } else{
    res = o_ht
  }
  return(res)
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

tempf_nn = function(x,mu,sig,tnu){
  f = function(y){
    d1 = 1 - pnorm(log((x/y)^2)-mu,mean = 0,sd = sig)
    d2 = dnorm(y)
    return(d1 * d2)
  }
  return(f)
}

tempf2 = function(x,tmu,tsig,tnu,p,alpha,eta_xi,eta_beta,u0){
  f = function(y){
    d1 = p*(1+eta_xi*(log((x/y)^2)-u0) / eta_beta)^(-1/eta_xi)
    #d1 = u0 +(eta_beta/eta_xi) * (((1-alpha)/p)^(-eta_xi) -1) # upper quantile of GPD 
    d2 = dnorm(y)
    return(d1 * d2)
  }
  return(f)
}

varx = function(tmu,sig,tnu,lev,accuracy,p){
  temp = rstd(2000,mean = tmu,sd = sig, nu = tnu)
  u0 = sort(temp,decreasing=F)[p*length(temp)+1]
  eta_gpd = gpd.fit(temp,threshold=u0,show=F)
  eta_beta=eta_gpd$mle[1]; eta_xi = eta_gpd$mle[2] # MLE's of scale and shape parameters
  left = 0; right = 1
  diff = right-left
  while(diff > accuracy){
    mid = (left + right)/2
    if(eta_xi>0){
      up1 = mid / exp(u0/2)
      up2 = mid / exp((u0 - eta_beta / eta_xi)/2)
      up = max(up1,up2)
      mp = integrate(tempf2(x = mid,tmu = tmu, tsig = tsig, 
                            tnu = tnu,p,1-lev,eta_xi,eta_beta,u0),0,Inf)$value
    } else{
      up = mid / exp(u0/2)
      low = mid / exp((u0 - eta_beta / eta_xi)/2)
      mp = integrate(tempf2(x = mid,tmu = tmu, tsig = tsig, 
                            tnu = tnu,p,1-lev,eta_xi,eta_beta,u0),low,Inf)$value
    }
    if(mp > lev){
      left = mid
    } else{
      right = mid
    }
    diff = right - left
  }
  return(mid)
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

varx_nn = function(tmu,sig,lev,accuracy){
  left = 0; right = 1
  diff = right-left
  while(diff > accuracy){
    mid = (left + right)/2
    mp = integrate(tempf_nn(x = mid,mu = tmu, sig = sig),0,Inf)$value
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


# ===========================================================
# Test of equal forecasting performance between two forecasts
# Reference: Gneiting & Ranjan (2011), 
# J. Business and Economic Statistics 29(3), pp.411-422

# k=1 # lead time for the forecast; only 1-step-ahead forecasts are considered

# Inputs: 
# s1 and s2: scoring functions for method 1 and 2
# ===========================================================


efp.test <- function(s1,s2,type="two-sided")
{
  d=s1-s2 # differences in scores
  n=length(d)
  
  # HAC estimator for the asymptotic variance of the average relative scores
  # HAC: heteroskedasticity and autocorrelation-consistent variance estimator
  # using Parzen window as the lag window (set of weights)
  # M: trancation point, set to approx. 2sqrt(n)	
  m=ceiling(2*sqrt(n))
  gam = acf(d,lag.max=m,type="covariance",plot=F)$acf
  k1 = 1:ceiling(m/2)
  k2 = (ceiling(m/2)+1):m
  
  lam = c(1, 2*(1-6*(k1/m)^2+6*(k1/m)^3),2*2*(1-k2/m)^3)
  sn = sqrt(gam %*% lam)
  
  #	sn = sqrt(sum(d^2)/n) # assumes zero correlation
  
  # test statistic
  tn=sqrt(n)*(mean(d))/sn
  if (type == "two-sided")
  {
    pv=2*ifelse(tn<0,pnorm(tn),(1-pnorm(tn))) #p-value for two-sided hypothesis
  }
  if (type == "one-sided-ge")
  {
    pv= pnorm(tn)
  }
  if (type == "one-sided-le")
  {
    pv= 1-pnorm(tn)
  }
  #	cat("mean=", mean(d),"; var=", var(d), ";sn=",sn, "\n")
  return(list(tn=tn,pvalue=pv,sn=sn))	
}

# Scoring function for VaR (1-homog. and 0-homog.).
# Inputs: 
# r: forecast
# x: verifying observation
# a: risk measure confidence level 
# h: parameter to identify which homogeneous score function to use. There are two choices: h=1 (b-homog. for a certain b>0), h=0 (zero-homog.)

# VaR-alpha

# h=1: "standard" choice, 1-homogeneous
sfVaR <- function(r,x,a,h=1)
{
  s <- vector(mode="numeric",length=length(r))
  ind1 = (x>r) & (r>0)
  ind2 = (x<r) & (r>0)
  if(h==1){
    s[ind1] <- -a*r[ind1] + x[ind1]
    s[ind2] <- (1-a)*r[ind2]
  }
  if(h==0){
    s[ind1] <- -a*log(r[ind1]) + log(x[ind1])
    s[ind2] <- (1-a)*log(r[ind2])
  }
  return(s)
}

