sk_nt = as.matrix(sk_nt)
res = as.numeric(sk_nt[,2:5])
res = matrix(res, ncol = 4)
par(mfrow = c(2,2))
hist(res[,1],col=rgb(1,0,0,0.5), main= expression(beta[0]),
     xlim = c(-1.3,-0.05),
     breaks = 15,ylim = c(0,80),xlab = '')
hist(res[,5],col=rgb(0,0,1,0.5), add = T, breaks = 15)
abline(v = omega)
legend('topleft',c('New Method','stochvol'),
       col = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),lwd = c(5,5))
hist(res[,2],col=rgb(1,0,0,0.5), main = expression(beta[1]),
     xlim = c(0.87,0.99),
     breaks = 15,ylim = c(0,80),xlab = '')
hist(res[,6],col=rgb(0,0,1,0.5), add = T, breaks = 15)
abline(v = phi)
legend('topleft',c('New Method','stochvol'),
       col = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),lwd = c(5,5))
hist(res[,3],col=rgb(1,0,0,0.5), main = expression(delta),
     xlim = c(0.15,0.66),
     breaks = 15,ylim = c(0,80),xlab = '')
hist(res[,7],col=rgb(0,0,1,0.5), add = T, breaks = 10)
abline(v = 0.35)
legend('topright',c('New Method','stochvol'),
       col = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),lwd = c(5,5))
hist(res[,4],col=rgb(1,0,0,0.5), main = expression(nu),
     xlim = c(0,25),
     breaks = 15,ylim = c(0,100),xlab = '')

