require("Hmisc")

to.stdnorm = function(dat){
  ind = rep(0,length(dat))
  cdfs = Ecdf(dat, pl = F)
  tempx = cdfs$x[-1]
  for(i in 1:length(dat)){
    ind[i] = which(tempx == dat[i])
  }
  fx = cdfs$y[ind]
  fx[which.min(fx)]=fx[which.min(fx[-which.min(fx)])]
  res = qnorm(fx)
  return(res)
}

find.beta =function(su,dat,b0){
  diff = 1
  b = b0
  while (diff > 1e-3) {
    b_old = b
    d = dgamma(b,1)+(b-1)*sum(log(dat))-sum(dat)-su
    n = trigamma(b)+sum(log(dat))
    b = b - d/n
    diff = abs(b_old - b)
  }
  return(b)
}

alpha.mcmc = function(dat){
  N = 10000
  a = rep(0,N)
  a[1] = runif(1)
  b = rep(0,N)
  for(i in 2:N){
    #     ta = 2
    #     while(ta > 1 | ta < 0){
    #       ta = rnorm(1, mean = a[i-1],sd = 1)
    #     }
    ta = rbeta(1,1,a[i-1]*10)
    ta = runif(1)
    #b_old = (1- a[i-1]) * mean(dat)
    b_old = beta.ite(a[i-1],dat)
    l_old = dgamma(mean(dat),shape = b_old, rate = (1- a[i-1]))*dbeta(a[i-1],2,2)
    #b_new = (1 - ta) * mean(dat)
    b_new = beta.ite(ta,dat)
    l_new = dgamma(mean(dat),shape = b_new, rate = (1- ta))*dbeta(ta,2,2)
    u = runif(1)
    #rho = (l_new / l_old) * (dnorm(a[i-1],mean =ta) / dnorm(ta, mean =a[i-1]))
    rho = (l_new / l_old) * (dbeta(a[i-1],1,ta*10) / dbeta(ta,1,a[i-1]*10))
    r = min(rho,1)
    if(u < r){
      a[i] = ta
      b[i] = b_new
    } else{
      a[i] = a[i-1]
      b[i] = b_old
    }
  }
  return(c(a,b))
}

beta.ite = function(alpha,dat){
  diff = 1
  ite = 0
  while(diff > 10e-5 | ite > 10e3){
    alpha_old = alpha
    beta = (1- alpha^2) * var(dat)
    alpha = 1 - beta / mean(dat)
    if(alpha < 0){alpha = 0.1}
    diff = abs(alpha - alpha_old)
  }
  return(beta)
}

block.boot = function(dat,l){
  blocklen = floor(length(dat) / l)
  nblock = floor(length(dat)/blocklen)
  ind = sample(1:nblock,nblock,replace = T)
  res = c()
  for(i in 1:nblock){
    res = c(res,dat[((ind[i]-1)*400+1):(ind[i]*400)])
  }
  return(res)
}

est.vol = function(dat,k = 7){
  est_vol = rep(0,length(dat))
  for(j in (k+1):length(est_vol)){
    est_vol[j] = sd(x[(j-k):j])
  }
  est_vol = est_vol[-(1:k)]
  return(est_vol)
}

rsymgamma = function(n,sh,sc = 1){
  u1 = 2*(runif(n)>0.5)-1
  u2 = runif(n)
  g = qgamma(u2,shape = sh, scale = sc)
  rands = u1*g
  #if(n != 1){rands = (rands - mean(rands))/sd(rands)}
  return(rands)
}

dsymgamma = function(x,beta,sc = 1){
  if(x >= 0 ){
    r = 0.5*dgamma(x,shape = beta, scale = sc)
  } else{
    r = 0.5*dgamma(-x,shape = beta, scale = sc)
  }
  return(r)
}

ec2 = function(df){
  es1 = dstd(qstd(0.95,nu = df),nu = df)/(1-0.95)
  es2 = (df - 2 + qstd(0.95,nu = df)^2)/(df-1)
  es = es1*es2
  return(es)
}

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
  
  # sn = sqrt(sum(d^2)/n) # assumes zero correlation
  
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

ecpa.test <- function(s1, s2)
{
  d=s1-s2 # score differences
  n=length(d) # out-of-sample size
  
  h=rbind(rep(1,n-1),d[-n])
  q=dim(h)[1]
  
  d1=d[-1] # score differences used in test statistic computation
  
  beta = coef(lm(d1 ~ d[-n])) # regression coefficients to be used in decision rule if the test rejects H0
  
  zbar=(h%*%d1)/(n-1) # q x 1
  
  Omega = matrix(0, nrow=q,ncol=q)
  for (t in 1:(n-1))
  {		
    Omega = Omega + d1[t]^2*(h[,t]%*%t(h[,t]))	
  }
  Omega = Omega/(n-1)
  Omega.inv = solve(Omega)
  
  #test statistic
  tn <- (n-1)*t(zbar)%*% Omega.inv %*% zbar
  pv <- pchisq(tn, df=q, lower.tail=FALSE)
  
  ind=0 # indicator which method (1=row, 2=col) performs better based on the decision rule beta.hat %*% h[T]>/<0
  # prop: proportion when method 1 outperforms method 2
  
  ev = t(beta) %*% h 
  ind = (ev[1,n-1]<0) + 2*(ev[1,n-1]>0)
  prop = sum(ev<0)/(n-1)
  
  return(list(tn=tn,pvalue=pv, prop=prop, ind=ind))
}

prd = function(x){
  return(besselK(abs(x),0) / pi)
}

sump = function(ub){
  v = integrate(prd,lower = -Inf, upper = ub)$value
  return(v)
}

sparsePOT = function(x,n = 1000,p = 0.95,alpha = 0.1){
  k = alpha*length(x)
  u=sort(x,decreasing=T)[k+1] 
  sx = sample(x,size = n,replace = F)
  GPDfit = gpd.fit(sx,threshold=u,show=F)
  beta=GPDfit$mle[1]; xi=GPDfit$mle[2]
  qPOT = u + (beta/xi) * (((1-p)/alpha)^(-xi) -1)
  return(qPOT)
}

LRuc <- function(lam, h)
{
  n1 = sum(h); lamh=n1/length(h)
  mllh = sum(log(lamh*h+(1-lamh)*(1-h))) # maximum log-likelihood
  llh0 = sum(log(lam*h+(1-lam)*(1-h))) # log-likelihood under null
  LR = 2*(mllh - llh0)
  pval = pchisq(q=LR, df=1, lower.tail=F)
  return(list(LR=LR,pvalue=pval))
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

prdt = function(z){
  n = length(z)
  res = rep(0,n)
  for(i in 1:n){
    tf = function(x){
      dstd(x,nu = 5)*dstd(z/x,nu = 5)/(abs(x))
    }
    res[i] = 2*integrate(tf,lower = 1e-5, upper = Inf,subdivisions = 1e5)$value
  }
  return(res)
}

sump.t = function(ub){
  v = integrate(prdt,lower = -Inf, upper = ub)$value
  return(v)
}


###### realized volatility for SV model ##########
sv.var = function(dat,l,p=0.95,al=0.1){
  #rv1 = TTR::volatility(dat,calc = "yang.zhang", n = l) / sqrt(260)
  rv1 = TTR::volatility(dat,calc="close", n = l) / sqrt(260)
  trv1 = rv1[-c(1:l)]
  
  t1 = log(trv1[-length(trv1)])
  t2 = log(trv1[-1])
#   cp = dat[,4]
#   lor = diff(log(cp))
#   adlor = lor[l:length(lor)]
#   rop = adlor / trv1
#   df = fitdist(distribution = 'std',rop)$pars[3]
  
  m1 = lm(t2~t1)
  b1 = unname(coefficients(m1)[2])
  b0 = unname(coefficients(m1)[1])
  resid = m1$fitted.values - t2
  beta = mean(abs(resid))
  sig2 = trv1[length(trv1)]
  re = ana.var(p_level = (1-p),sig_t = sig2,xmin = 0,xmax = 0.2,xacc = 1e-4,
               sh = beta,alpha = b1,omega = b0, ncore = 0)
  return(re)
}

sv.vol = function(dat,l){
  #rv1 = TTR::volatility(dat,calc = "yang.zhang", n = l) / sqrt(260)
  rv1 = TTR::volatility(dat,calc="close", n = l) / sqrt(260)
  trv1 = rv1[-c(1:l)]
  
  t1 = log(trv1[-length(trv1)])
  t2 = log(trv1[-1])
  m1 = lm(t2~t1)
  b1 = unname(coefficients(m1)[2])
  b0 = unname(coefficients(m1)[1])
  plsig = b0 + b1*t2[length(t2)]
  psig = sqrt(exp(plsig))
  return(psig)
}

hs.var = function(dat,p){
  lt = diff(log(dat[,4]))
  hist_sim = sample(lt,1e6,replace = T)
  #hist_sim = lt
  hs_var = quantile(hist_sim,p)
  return(hs_var)
}

inner.std = function(xi,df){
  f = function(y){
    return(y^xi * dstd(y,nu = df))
  }
  return(f)
}

inner.ep = function(xi,ep){
  df = approxfun(density(ep))
  f = function(y){
    return(y^xi * df(y))
  }
  return(f)
}

inner.norm = function(xi){
  f = function(y){
    return(y^xi * dnorm(y))
  }
  return(f)
}

pot.var = function(dat,l,p,al = 0.1){
  rv1 = TTR::volatility(dat,calc="close", n = l) / sqrt(260)
  trv1 = rv1[-c(1:l)]
  z0 = trv1#log(trv1)
  ti = hill.est(z0,al)
  mu = mean(z0); sig = sd(z0)
  nz = (z0-mu) / sig
  z = nz#sample(nz,1e6,replace = T)
  k = al*length(z)
  u=sort(z,decreasing=T)[ceiling(k+1)]
  GPDfit = gpd.fit(z,threshold=u,show=F)
  beta=GPDfit$mle[1]; xi=GPDfit$mle[2] # MLE's of scale and shape parameters
  c = integrate(inner.std(ti),0,Inf)$value
  p2 = 1- (1-p)/c
  qPOT = u + (beta/xi) * (((1-p2)/al)^(-xi) -1)
  svar = qPOT*sig +mu
  return(c(svar,xi))
}

ep.var = function(dat,l,p,al = 0.1){
  rv1 = TTR::volatility(dat,calc="close", n = l) / sqrt(260)
  trv1 = rv1[-c(1:l)]
  z0 = trv1#log(trv1)
  ti = hill.est(z0,al)
  mu = mean(z0); sig = sd(z0)
  nz = (z0-mu) / sig
  z = sample(nz,1e6,replace = T)
  #xi = hillestimator(trv1,al)
  c = integrate(inner.std(ti),0,Inf)$value
  p2 = 1- (1-p)/c
  qr = quantile(z,p2)
  sqr = qr*sig + mu
  #r = exp(sqr)
  return(sqr)
}

mb.var = function(dat,l,p,al = 0.1){
  rv1 = TTR::volatility(dat,calc = "close", n = l) / sqrt(252)
  trv1 = rv1[-c(1:l)]
  
  t1 = log(trv1[-length(trv1)])
  t2 = log(trv1[-1])

  m1 = lm(t2~t1)
  b1 = unname(coefficients(m1)[2])
  b0 = unname(coefficients(m1)[1])
  resid = resid(m1)
  mr = mean(resid); sr = sd(resid)
  nr = (resid - mr) / sr
  df = fitdist(distribution = 'std',nr)$pars[3]
  #beta = mean(abs(resid))
  sig2 = trv1[length(trv1)]
  
  z0 = trv1#log(trv1)
  ti = hill.est(z0,al)
  c = integrate(inner.norm(ti),0,Inf)$value
  
  p2 = 1- (1-p)/c
  eq = qstd(p2,nu = df)*sr+mr
  ls2 = b0+b1*log(sig2) + eq
  qs = exp(ls2)
  return(qs)
}

mbexp.var = function(dat,l,p,al = 0.1){
  rv1 = TTR::volatility(dat,calc = "close", n = l) / sqrt(260)
  trv1 = rv1[-c(1:l)]
  
  t1 = log(trv1[-length(trv1)])
  t2 = log(trv1[-1])
  
  m1 = lm(t2~t1)
  b1 = unname(coefficients(m1)[2])
  b0 = unname(coefficients(m1)[1])
  resid = resid(m1)
  pr = resid[which(resid > 0)]
  df = unname(MASS::fitdistr(pr,"exponential")$estimate)
  #beta = mean(abs(resid))
  sig2 = trv1[length(trv1)]
  
  z0 = trv1#log(trv1)
  ti = hill.est(z0,al)
  c = integrate(inner.norm(ti),0,Inf)$value
  
  p2 = 1- (1-p)/c
  p3 = 1 - 2*(1-p2)
  eq = qexp(p3,rate = df)
  ls2 = b0+b1*log(sig2) + eq
  qs = exp(ls2)
  return(qs)
}

mbpot.var = function(dat,l,p,al = 0.1){
  rv1 = TTR::volatility(dat,calc = "close", n = l) / sqrt(260)
  trv1 = rv1[-c(1:l)]
  
  t1 = log(trv1[-length(trv1)])
  t2 = log(trv1[-1])
  
  m1 = lm(t2~t1)
  b1 = unname(coefficients(m1)[2])
  b0 = unname(coefficients(m1)[1])
  resid = resid(m1)
  sig2 = trv1[length(trv1)]
  
  z0 = trv1#log(trv1)
  ti = hill.est(z0,al)
  c = integrate(inner.norm(ti),0,Inf)$value
  
  #rs = sd(resid)
  nr = resid
  k = al*length(nr)
  u=sort(nr,decreasing=T)[ceiling(k+1)]
  GPDfit = gpd.fit(nr,threshold=u,show=F)
  beta=GPDfit$mle[1]; xi=GPDfit$mle[2] # MLE's of scale and shape parameters
  p2 = 1- (1-p)/c
  eq = u + (beta/xi) * (((1-p2)/al)^(-xi) -1)
  ls2 = b0+b1*log(sig2) + eq
  qs = exp(ls2)
  return(qs)
}

g.var.2 = function(dat,p,l){
  tret = diff(log(dat))
  vols = TTR::volatility(dat,calc="close", n = l) / sqrt(260)
  vols = vols[-c(1:(l-1))]
  h1 = vols[-1]^2
  h2 = vols[-length(vols)]^2
  r2 = (tret[(l-1):(length(tret)-1)])^2
  m = lm(h1~h2+r2)
  nd = data.frame(h2 = vols[length(vols)]^2,r2 = tret[length(tret)]^2)
  vh1 = predict(m,nd)
  var = sqrt(vh1)*qstd(p)
  return(var)
}

garch.var = function(dat,p,al = 0.1){
  if(is.matrix(dat)){
    lt = diff(log(dat[,4]))
  } else{
    lt = diff(log(dat))
  }
  t_spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE),
                      variance.model = list(garchOrder = c(1,1)),
                      distribution.model='std')
  t_fit = ugarchfit(t_spec, lt)
  resid=residuals(t_fit) 
  sigma.t=sigma(t_fit)
  z=as.numeric(resid/sigma.t) # residuals
  frcst = ugarchforecast(t_fit,n.ahead=1)
  sigt1 = sigma(frcst) # expected sigma for next day
  mut1 = fitted(frcst)
  k = al*length(z)
  u=sort(z,decreasing=T)[ceiling(k+1)]
  GPDfit = gpd.fit(z,threshold=u,show=F)
  beta1=GPDfit$mle[1]; xi1=GPDfit$mle[2] # MLE's of scale and shape parameters
  qPOT = u + (beta1/xi1) * (((1-p)/al)^(-xi1) -1)
  return(mut1+sigt1*qPOT)
}

garch.vol = function(dat){
  if(ncol(dat) == 4){
    lt = diff(log(dat[,4]))
  } else{
    lt = diff(log(dat))
  }
  
  t_spec = ugarchspec(mean.model=list(armaOrder=c(0,0),include.mean=FALSE),
                      variance.model = list(garchOrder = c(1,1)),
                      distribution.model='norm')
  t_fit = ugarchfit(t_spec, lt)
  frcst = ugarchforecast(t_fit,n.ahead=1)
  sigt1 = as.numeric(sigma(frcst))
  return(sigt1)
}

condemp = function(x,sig_t,sh,phi,alpha){
  f = function(y){
    d2 = dnorm(y)
    ds = rsymgamma(1000,sh = sh,sc = 1)
    po_s = exp(alpha*log(sig_t)+phi+ds)
    d1 = (length(which(po_s > x/y)) / length(po_s))
    return(d1*d2)
  }
  return(f)
}

int.y2 = function(x_val,sig_t,sh,phi,alpha){
  r = integrate(condemp(x_val,sig_t,sh,phi,alpha),0,Inf)$value
  return(r)
} 

int.y = function(x_val,sig_t,sh,phi,alpha,seed = 1,acc = 0.001){
  set.seed(seed)
  d1 = rsymgamma(1000,sh = sh,sc = 1)
  po_s = exp(alpha*log(sig_t)+phi+d1)
  y = seq(0.01,10,by = acc)
  r = 0
  for(i in 1:length(y)){
    r = r + (length(which(po_s > x_val/y[i])) / length(po_s)) * dnorm(y[i]) * acc
  }
  return(r)
}


int.var = function(p_level,xmin = 0, xmax = 5, xacc = 0.05,
                   sig_t,sh,phi,alpha,seed = 1,acc = 0.001,ncore = 10){
  x = seq(xmin,xmax,by = xacc)
  if(ncore == 0){
    rs = rep(0,length(x))
    for(j in 1:length(x)){
      rs[j]=int.y(x[j],sig_t,sh,phi,alpha,seed,acc)
    }
  } else{
    registerDoMC(ncore)
    rs = foreach(j = 1:length(x))%dopar%{
      int.y(x[j],sig_t,sh,phi,alpha,seed,acc)
    }
    rs = as.numeric(rs)
  }
  ind = which.min(abs(rs-p_level))
  if(rs[ind]>p_level){ind = ind+1}
  return(x[ind])
}

condepi.normal = function(x,sig_t,sh,alpha,omega){
  f = function(y){
    d1 = dnorm(y)
    t = log(x/y)-alpha*log(sig_t)-omega
    d2 = 0.5-0.5*pgamma(t,shape = sh,rate = 1)
    return(d1*d2)
  }
  return(f)
}

condepi.std = function(x,sig_t,sh,alpha,omega){
  f = function(y){
    d1 = dstd(y,nu = 3)
    t = log(x/y)-alpha*log(sig_t)-omega
    d2 = 0.5-0.5*pgamma(t,shape = sh,rate = 1)
    return(d1*d2)
  }
  return(f)
}

ana.y = function(x_val,sig_t,sh,phi,alpha){
  r = integrate(condepi.normal(x_val,sig_t,sh,phi,alpha,omega),0,Inf)$value
  return(r)
}


ana.var = function(p_level,xmin = 0, xmax = 5, xacc = 0.05,
                   sig_t,sh,alpha,omega,seed = 1,acc = 0.001,ncore = 10){
  x = seq(xmin,xmax,by = xacc)
  if(ncore == 0){
    rs = rep(0,length(x))
    for(j in 1:length(x)){
      rs[j] = integrate(condepi.normal(x[j],sig_t,sh,alpha,omega),0,Inf)$value
    }
  }else{
    registerDoMC(ncore)
    rs = foreach(j = 1:length(x),.combine = rbind)%dopar%{
      integrate(condepi.std(x[j],sig_t,sh,alpha),0,Inf)$value
    }
  }
  
  #rs = as.numeric(rs)
  ind = which.min(abs(rs-p_level))
  if(rs[ind]>p_level){ind = ind+1}
  return(x[ind])
}

qlike.score = function(x,y){
  return(log(x)+y/x)
}

mse.score = function(x,y){
  return((x-y)^2)
}

condepi.t = function(x,sig_t,df,alpha,omega){
  f = function(y){
    d1 = dstd(y,nu = df)
    t = log(x/y)-alpha*log(sig_t)-omega
    d2 = 1 - pstd(t,nu = df)
    return(d1*d2)
  }
  return(f)
}

ana.t = function(p_level,xmin = 0, xmax = 5, xacc = 0.05,
                   sig_t,df,alpha,omega,seed = 1,acc = 0.001,ncore = 10){
  x = seq(xmin,xmax,by = xacc)
  if(ncore == 0){
    rs = rep(0,length(x))
    for(j in 1:length(x)){
      rs[j] = integrate(condepi.t(x[j],sig_t,df,alpha,omega),0,Inf)$value
    }
  }else{
    registerDoMC(ncore)
    rs = foreach(j = 1:length(x),.combine = rbind)%dopar%{
      integrate(condepi.t(x[j],sig_t,df,alpha,omega),0,Inf)$value
    }
  }
  #rs = as.numeric(rs)
  ind = which.min(abs(rs-p_level))
  if(rs[ind]>p_level){ind = ind+1}
  return(x[ind])
}

hill.est = function(x,al){
  u = quantile(x,1-al)
  exc = x[which(x>u)]
  hill = 1/(sum(log(exc)-log(u))/length(exc))
  return(hill)
}

tf = function(epi,re,x,sig2,al,b0,b1){
  etat = exp(re)
  k = floor(al*length(etat)+1)
  u = unname(sort(etat,decreasing = T)[k])
  GPDfit = gpd.fit(etat,threshold=u,show=F)
  beta=GPDfit$mle[1]; xi=GPDfit$mle[2]
  epif = approxfun(density(epi))
  eptat = approxfun(density(etat))
  f = function(y){
    dk = exp(b0 + b1*log(sig2))
    df = x / (dk*y)
#     p1 = rep(0,length(y))
#     for(z in 1:length(y)){
#       if(df[z] <u){
#         p1[z] = 0
#       } else{
#         q = 1+(xi/beta)*(df[z]-u)
#         p1[z] = (k/length(etat))*(q)^(-1/xi)
#       }
#     }
    integrate(eptat,df)
    p2 = epif(y)
    return(p1*p2)
  }
  return(f)
}

ana.v = function(xmin,xmax,xacc,
                 epi,re,x,sig2,al,b0,b1,p){
  x = seq(xmin,xmax,xacc)
  ps = rep(0,length(x))
  for(l in 1:length(x)){
    tx = x[l]
    y= tf(epi = epi,re = re,x[l],
               sig2 = sig2,al = al,b0 = b0,b1 = b1)
    ps[l] = integrate(y,0.01,3,stop.on.error = F)$value
    #print(l)
  }
}

new.inner = function(resid,K,df,dfr){
  rs = rstd(1e5,mean= mean(resid),sd = sd(resid),
            nu = dfr)
  drs = approxfun(density(exp(rs)))
  lower = round(min(exp(rs)),5)+1e-5
  upper = round(max(exp(rs)),5)-1e-5
  f = function(y){
    pr = drs(y)
    pe = (K*y)^(df)
    return(pr*pe)
  }
  return(list(f,lower,upper))
}

vol.pred = function(dat,l){
  rvs = TTR::volatility(dat,calc="close", n = l) / sqrt(260)
  rvs = rvs[-c(1:l)]
  h1 = log(rvs^2)[-1]
  h2 = log(rvs^2)[-length(rvs)]
  m1 = lm(h1~h2)
  pred = predict(m1,newdata = data.frame(h2 = h1[length(h1)]))
  return(sqrt(exp(pred)))
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