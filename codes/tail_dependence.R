library(evd)
library(fGarch)
source("new_sv_funs.R")

set.seed(1)
N = 50000
burnin = 5000
phi = 0.95
t_nu = 5
omega = -0.5

err1 = rnorm(N)
#err1 = rstd(N,mean = 0,sd = 1,nu = epi_nu) # epsilon_t
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
par(mfrow = c(1,2))
chiplot(smdat,nq = 200)
par(mfrow = c(1,1))

# calculate chi_bar
chi.bar(x,q = 0.9995)
