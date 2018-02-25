# MNP- clase 7 feb - Prueba-signo
#prueba-signo para hipotesis a una cola
#Rechazar H0 si S=suma de signos + en la diferencia es mayor

mu0=28
sigma=1
n=16
l=5*sigma/sqrt(n) 
# Longitud grande (5 desviaciones estandar)
mu1.U = seq(from= mu0, to= mu0+l, by= 0.005) ; mu1.U
z=(mu0-mu1.U)/sigma ;z
theta1=1-pnorm(z) ;theta1
potenciaSigno=1-pbinom(11, 16, theta1, lower.tail = TRUE) # Equivalente a p(s=>12|theta1)
plot(mu1.U, potenciaSigno)

#Comparacion con la prueba t
alpha=0.05
z.alpha= qnorm(alpha, lower.tail = FALSE) ; z.alpha 
c1=mu0+z.alpha*sigma/sqrt(n) ; c1 # Limite para rechazar H0
zt=(c1-mu1.U)/sigma*sqrt(n) ;zt
potencia.tB=1-pnorm(zt) ; potencia.tB

plot(mu1.U, potencia.tB)
lines(mu1.U, potenciaSigno)
rm(list=ls())
