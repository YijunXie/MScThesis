library(evd)
library(fGarch)
library(ismev)
library(extRemes)
source("new_sv_funs.R")

foo <- function(u){
  # fit<-gpd.fit(tvec, threshold=quantile(tvec,u), show=FALSE)
  fit<-gpd.fit(tvec, threshold=quantile(tvec,u), show=F)
  xi.pro <- gpd.parameterCI(fit, conf=.95, xi.xlow=0.15, xi.xup=1.5,xi.only=TRUE)$xi
  c(xi.pro$dn, fit$mle[2], xi.pro$up)
}

eta <- function(x){
  bar <- data.frame(t(sapply(x,foo)))
  names(bar) <- c("dn","mle","up")
  bar	
}

set.seed(8) ## t innovation
#set.seed(9) ## N innovation
N = 10000
burnin = 2500
phi = 0.95
t_nu = 3
epi_nu = 5
omega = -0.5

err1 = rnorm(N)
#err1 = rstd(N,mean = 0,sd = 1,nu = epi_nu) # epsilon_t
#err2 = rnorm(N,sd = 0.35)
err2 = rstd(N,mean = 0,sd = 0.5,nu = t_nu) # eta_t
x = rep(0,N)
sig = rep(0,N)
sig[1] = 0.01
for(i in 2:N){
  sig[i] = sqrt(exp(omega+phi * log(sig[i-1]^2) + err2[i]))
  x[i] = sig[i]*err1[i]
}

x = x[-c(1:burnin)]
sig = sig[-c(1:burnin)]

smdat = cbind(x[-length(x)],x[-1])
#par(mfcol = c(2,3))
chiplot(smdat,nq = 1000,which = 1,xlim = c(0.8,1),
        ylim1 = c(-0.1,1),ylim2 = c(-0.1,1),conf = 0.95, 
        main1 = c(expression(chi)," Plot for ARSV Process"))
abline(h = 0,lty = 'dotdash')
abline(h = 1,lty = 'dotdash')
#par(mfrow = c(1,1))
chiplot(smdat,nq = 1000,which = 2,xlim = c(0.8,1),
        ylim1 = c(-0.1,1),ylim2 = c(-0.1,1),conf = 0.95,
        main2 = expression(bar(chi)))
abline(h = 0,lty = 'dotdash')
abline(h = 1,lty = 'dotdash')
#par(mfrow = c(1,1))
# calculate chi_bar
# chi.bar(x,q = 0.999,n = 1)


# n=nrow(smdat)
# 
# u1=rank(smdat[,1])/(n+1); u2=rank(smdat[,2])/(n+1)
# z1=-1/log(u1); z2=-1/log(u2)
# tvec=apply(cbind(z1,z2), 1, min)
# 
# urange <- seq(0.5,0.97, length=200)
# plotdat <- cbind(u=urange, eta(urange))
# 
# plot(plotdat$u, plotdat$mle, type="l",ylim=c(0,1.4),xlab="threshold quantile", ylab=expression(hat(eta)))
# lines(plotdat$u, plotdat$dn, lty=2)
# lines(plotdat$u, plotdat$up, lty=2)
# abline(h=0.5, col=2)
# abline(h=1,col = 2)

# for GARCH
set.seed(1)
N = 10000
burnin = 2500
phi = 0.8
beta = 0.1
t_nu = 5
epi_nu = 5
omega = 1e-6

err1 = rnorm(N)

gx = rep(0,N)
gsig = rep(0,N)
gsig[1] = 0.01
gx[1] = 0
for(i in 2:N){
  gsig[i] = sqrt(omega + beta * gsig[i-1]^2 + phi * gx[i-1]^2)
  gx[i] = gsig[i]*err1[i]
}

gx = gx[-c(1:burnin)]
gsig = gsig[-c(1:burnin)]

gsim = cbind(gx[-length(gx)],gx[-1])

#par(mfrow = c(1,1))
#par(mfrow = c(1,2))
chiplot(gsim,nq = 1000,which = 1,xlim = c(0.8,1),
        ylim1 = c(-0.1,1),ylim2 = c(-0.1,1),conf = 0.95, 
        main1 = expression(chi))
abline(h = 0,lty = 'dotdash')
abline(h = 1,lty = 'dotdash')
#par(mfrow = c(1,1))
chiplot(gsim,nq = 1000,which = 2,xlim = c(0.8,1),
        ylim1 = c(-0.1,1),ylim2 = c(-0.1,1),conf = 0.95,
        main2 = expression(bar(chi)))
abline(h = 0,lty = 'dotdash')
abline(h = 1,lty = 'dotdash')
#par(mfrow = c(1,1))
# # calculate chi_bar
# chi.bar(gx,q = 0.999,n = 1)
# 
# # Estimation of eta
# n=nrow(gsim)
# w = rep(0,n)
# for(i in 1:n){
#   w[i] = min(gsim[,i])
# }
# 


u1=rank(gsim[,1])/(n+1); u2=rank(gsim[,2])/(n+1)
z1=-1/log(u1); z2=-1/log(u2)
tvec=apply(cbind(z1,z2), 1, min)

urange <- seq(0.5,0.97, length=1000)
plotdat <- cbind(u=urange, eta(urange))

plot(plotdat$u, plotdat$mle, type="l",ylim=c(.2,1.4),xlab="threshold quantile", ylab=expression(hat(eta)))
lines(plotdat$u, plotdat$dn, lty=2)
lines(plotdat$u, plotdat$up, lty=2)
abline(h=1, col=2)
abline(h = 0.5,col = 2)
