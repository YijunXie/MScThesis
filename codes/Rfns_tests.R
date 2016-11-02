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

