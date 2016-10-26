library(evd)
library(fGarch)
source("new_sv_funs.R")

set.seed(2)
N = 50000
burnin = 5000
phi = 0.95
t_nu = 5
epi_nu = 5
omega = -0.5

err1 = rnorm(N)
#err1 = rstd(N,mean = 0,sd = 1,nu = epi_nu) # epsilon_t
#err2 = rnorm(N)
err2 = rstd(N,mean = 0,sd = 0.35,nu = t_nu) # eta_t
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
#par(mfrow = c(1,2))
chiplot(smdat,nq = 200,which = 2,main2 = "Chi Bar Plot for ARSV Process")
abline(h = 1,lty = 'dotdash')
#par(mfrow = c(1,1))

# calculate chi_bar
chi.bar(x,q = 0.999,n = 1)

# for GARCH
set.seed(2)
N = 50000
burnin = 5000
phi = 0.6
beta = 0.3
t_nu = 5
epi_nu = 5
omega = 0.02

err1 = rnorm(N)
#err1 = rstd(N,mean = 0,sd = 1,nu = epi_nu) # epsilon_t
#err2 = rnorm(N)
err2 = rstd(N,mean = 0,sd = 0.35,nu = t_nu) # eta_t
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
#par(mfrow = c(1,2))
chiplot(gsim,nq = 200,which = 2, main2 = "Chi Bar Plot for GARCH Process")
abline(h = 1,lty = 'dotdash')
#par(mfrow = c(1,1))

# calculate chi_bar
chi.bar(gx,q = 0.999,n = 1)

