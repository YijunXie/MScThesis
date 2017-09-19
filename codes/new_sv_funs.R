"bisearch" <-
  function (x1, x2, f, tol = 1e-07, niter = 25, upcross.level = 0) 
  {
    f1 <- f(x1) - upcross.level
    f2 <- f(x2) - upcross.level
    if (f1 > f2) 
      stop(" f1 must be < f2 ")
    iter <- niter
    for (k in 1:niter) {
      xm <- (x1 + x2)/2
      fm <- f(xm) - upcross.level
      if (fm < 0) {
        x1 <- xm
        f1 <- fm
      }
      else {
        x2 <- xm
        f2 <- fm
      }
      if (abs(fm) < tol) {
        iter <- k
        break
      }
    }
    xm <- (x1 + x2)/2
    fm <- f(xm) - upcross.level
    list(x = xm, fm = fm, iter = iter)
  }


"gpd.parameterCI" <-
  function (z, m, conf = 0.95, nint = 100, rl.xup=NULL, rl.xlow=NULL, xi.xup=NULL, xi.xlow=NULL, rl.only=FALSE, xi.only=FALSE,
            make.plot=FALSE)
  {
    ###
    ## Function to try to find confidence limits (via profile likelihoods) for the GPD
    ## parameters (shape and return level).
    #
    # 'z' is an object of class "gpd.fit" created by the 'gpd.fit' function of the 'ismev' package.
    # 'm' is the 'm'-year return level to find CIs for.
    # 'conf' is the alpha confidence level.
    # 'nint' are the number of points with which to find the likelihood value (and CIs for).
    # 'rl.only' if TRUE only the return level parameter will be estimated.
    # 'xi.only' if TRUE only the shape parameter will be estimated.
    # 'make.plot' if TRUE will plot the profile likelihoods.
    #
    # Returns a list with components:
    # 'rl' and 'xi', each of which are lists containing info on the respective confidence intervals.
    # 	Including 'up' and 'dn' (the upper and lower limits), and 'sfun', the spline fit function
    #	used to create profile plot.
    #
    lmts <- list()
    xdat <- z$data
    ma <- -z$nllh
    u <- z$threshold
    la <- z$rate
    npy <- z$npy
    v <- numeric(nint)
    eps <- 1e-6
    if( make.plot) {
      if( !rl.only & !xi.only) par( mfrow=c(2,1))
      else par( mfrow=c(1,1))
    }
    
    if( !xi.only) {
      est.rl.xup <- is.null( rl.xup)
      est.rl.xlow <- is.null( rl.xlow)
      mstar <- m
      m <- m*npy
      rl.mle <- gpd.ret( z, mstar)
      q <- gpdq2( z$mle, u, la, m)
      d1 <- rep(0, length(q))
      d2 <- (gpdq2(c( z$mle[1]+eps, z$mle[2]), u, la, m) - q)/eps
      d3 <- (gpdq2(c( z$mle[1], z$mle[2]), u, la, m) - q)/eps
      d <- cbind(d1, d2, d3)
      mat <- z$cov
      mat <- matrix(c((la * (1 - la))/z$n, 0, 0, 0, mat[1, 1], mat[1, 2], 0, mat[2, 1], mat[2, 2]), ncol = 3)
      vv <- apply(d, 1, q.form, m = mat)
      if( est.rl.xlow) rl.xlow <- rl.mle - 1.5*qnorm((1-conf)/2, lower.tail=FALSE)*sqrt(vv)
      if( est.rl.xup) rl.xup <- rl.mle + 1.5*qnorm((1-conf)/2, lower.tail=FALSE)*sqrt(vv)
      x <- seq(rl.xlow, rl.xup, length = nint)
      sol <- z$mle[2]
      gpd.plik <- function(a) {
        if (m != Inf) 
          sc <- (a * (xp - u))/((m * la)^a - 1)
        else sc <- (u - xp)/a
        if (abs(a) < 10^(-4)) 
          l <- length(xdat) * log(sc) + sum(xdat - u)/sc
        else {
          y <- (xdat - u)/sc
          y <- 1 + a * y
          if (any(y <= 0) || sc <= 0) 
            l <- 10^6
          else l <- length(xdat) * log(sc) + sum(log(y)) * 
            (1/a + 1)
        }
        l
      }
      for (i in 1:nint) {
        xp <- x[i]
        opt <- optim(sol, gpd.plik, method = "BFGS")
        sol <- opt$par
        v[i] <- opt$value
      }
      lmts$upcross.level <- ma-0.5*qchisq(conf,1)
      lmts$rl$mle <- rl.mle
      sfun <- splinefun(x, -v)
      lmts$rl$sfun <- sfun
      # if( !est.rl.xup) {
      x1 <- order( sfun( c(rl.mle, rl.xup)))
      if( x1[1] == 1) lmts$rl$up <- bisearch( rl.mle, rl.xup, f=sfun, upcross.level=ma-0.5*qchisq(conf,1))$x
      else lmts$rl$up <- bisearch( rl.xup, rl.mle, f=sfun, upcross.level=ma-0.5*qchisq(conf,1))$x
      # } else lmts$rl$up <- findparlims( sfun, upcross.level=ma-0.5*qchisq(conf,1), guess=rl.xup)
      # if( !est.rl.xlow) {
      x2 <- order( sfun( c(rl.mle, rl.xlow)))
      if( x2[1] == 1) lmts$rl$dn <- bisearch( rl.mle, rl.xlow, f=sfun, upcross.level=ma-0.5*qchisq(conf,1))$x
      else lmts$rl$dn <- bisearch( rl.xlow, rl.mle, f=sfun, upcross.level=ma-0.5*qchisq(conf,1))$x
      # } else lmts$rl$dn <- findparlims( sfun, upcross.level=ma-0.5*qchisq(conf,1), guess=rl.xlow)
      
      if( make.plot) {
        plot(x, -v, type = "l", xlab = "Return Level", ylab = "Profile Log-likelihood")
        abline(h = ma, col="green")
        abline(h = ma - 0.5 * qchisq(conf, 1), col="green")
        abline(v=c(lmts$rl$dn, lmts$rl$up),lty=2)
      } # end of if make.plot stmt.
    } # end of if !xi.only stmt.
    
    if( !rl.only) {
      # Now do the shape parameter...
      # cat("If routine fails (shape parameter), try changing plotting interval", fill = TRUE)
      est.xi.xup <- is.null( xi.xup)
      est.xi.xlow <- is.null( xi.xlow)
      xdat <- z$data
      xi.mle <- z$mle[2]
      if( est.xi.xlow) xi.xlow <- xi.mle - 1.5*qnorm((1-conf)/2, lower.tail=FALSE)*z$se[2]
      if( est.xi.xup) xi.xup <- xi.mle + 1.5*qnorm((1-conf)/2, lower.tail=FALSE)*z$se[2]
      u <- z$threshold
      v <- numeric(nint)
      x <- seq(xi.xup+qnorm((1-conf)/2, lower.tail=FALSE)*z$se[2],
               xi.xlow-qnorm((1-conf)/2, lower.tail=FALSE)*z$se[2], length = nint)
      sol <- z$mle[1]
      gpd.plikxi <- function(a) {
        if (abs(xi) < 10^(-4)) 
          l <- length(xdat) * log(a) + sum(xdat - u)/a
        else {
          y <- (xdat - u)/a
          y <- 1 + xi * y
          if (any(y <= 0) || a <= 0) 
            l <- 10^6
          else l <- length(xdat) * log(a) + sum(log(y)) * (1/xi + 1)
        }
        l
      }
      for (i in 1:nint) {
        xi <- x[i]
        opt <- optim(sol, gpd.plikxi, method = "BFGS")
        sol <- opt$par
        v[i] <- opt$value
      }
      
      sfun <- splinefun(x, -v)
      lmts$xi$sfun <- sfun
      x1 <- order( sfun( c(xi.mle, xi.xup)))
      x2 <- order( sfun( c(xi.mle, xi.xlow)))
      if( x1[1] == 1) lmts$xi$up <- bisearch( xi.mle, xi.xup, f=sfun, upcross.level=ma-0.5*qchisq(conf,1))$x
      else lmts$xi$up <- bisearch( xi.xup, xi.mle, f=sfun, upcross.level=ma-0.5*qchisq(conf,1))$x
      if( x2[1] == 1) lmts$xi$dn <- bisearch( xi.mle, xi.xlow, f=sfun, upcross.level=ma-0.5*qchisq(conf,1))$x
      else lmts$xi$dn <- bisearch( xi.xlow, xi.mle, f=sfun, upcross.level=ma-0.5*qchisq(conf,1))$x
      
      if( make.plot) {
        plot(x, -v, type = "l", xlab = "Shape Parameter", ylab = "Profile Log-likelihood")
        abline(h = ma, lty = 1, col="green")
        abline(h = ma - 0.5 * qchisq(conf, 1), col="green")
        abline(v=c(lmts$xi$dn, lmts$xi$up), lty=2)
      }
    } # end of if !rl.only stmt
    class( lmts) <- "gpd.parameterCI.obj"
    invisible(lmts)
  }

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

mh.dn = function(o_ht,n_ht,x,b0,b1,sd,h2){
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


# mh.ht = function(o_ht,n_ht,x,b0,b1,sd,df,h2,oh2,old_theta){
#   o_b0 = old_theta[1]
#   o_b1 = old_theta[2]
#   o_sd = old_theta[3]
#   o_df = old_theta[4]
#   u = runif(1)
#   pn1 = dnorm(x,sd = sqrt(exp(n_ht)))
#   pn2 = dstd(h2,mean =  b0 + b1*n_ht, sd = sd,nu = df)
#   po1 = dnorm(x,sd = sqrt(exp(o_ht)))
#   po2 = dstd(oh2,mean =  o_b0 + o_b1*o_ht, sd = o_sd,nu = o_df)
#   a = (pn1*pn2)/(po1*po2)
#   a = min(1,a)
#   if(a > u){
#     res = n_ht
#   } else{
#     res = o_ht
#   }
#   return(res)
# }

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
    d1 = 1 - pstd((log((x/y)^2)-mu)/sig,mean = 0,sd = 1, nu = tnu)
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


# ===================================================
# Test of equal conditional predictive ability (ECPA test)
# Reference: Giacomini and White (2006, Econometrica)

# tau=1 # lead time for the forecast; only 1-step-ahead forecasts are considered
# ECPA test statistic is given in Eqn (4)

# Inputs: 
# s1 and s2: scoring functions for method 1 and 2 (n x 1 vectors)
# h: test function (in q x n matrix)
# currently, h[t]=(1, d[t])' where d is score difference

# Outputs:
# tn: test statistic
# pvalue
# prop: proportion of the times in the out-of-sample when method 1 outperforms method 2
# ind: indicator which method (1=row, 2=col) performs better based on the decision rule beta.hat %*% h[T]>/<0
# ===================================================
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

