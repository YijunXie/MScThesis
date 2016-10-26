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
