

#mu=0.1
sigma_n=0.1
phi=0.98
k1=sqrt(1-phi*phi)/sigma_n


m=seq(25, 250, by=1)
b1=2.5
b2=4
b3=5

c11=2*b1*k1/m
c12=2*b2*k1/m
c13=2*b3*k1/m

#mu=0.5
sigma_n=0.15
phi=0.99
k2=sqrt(1-phi*phi)/sigma_n

c21=2*b1*k2/m
c22=2*b2*k2/m
c23=2*b3*k2/m

par(mfrow=c(1, 2))
plot(m, c11, type='l', ylim=c(0,.85))
lines(m, c12, lty=2)
lines(m, c13, lty=3)

plot(m, c21, type='l', ylim=c(0,.85))
lines(m, c22, lty=2)
lines(m, c23, lty=3)
par(mfrow=c(1, 1))




