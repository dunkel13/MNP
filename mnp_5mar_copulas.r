library(MASS)
library(mvtnorm)
library(cubature)
n=100

#Generación de una muestra normal multivariada de media "mu" y varianza y cov
s11=1; s22=1; s12=0.0
cov=cbind(c(s11, s12), c(s12, s22)); #Matriz de var-cov (matriz de correlación)
mu<- c(0,0); # Vector de medias (normal estandar bivariada)


rho=s12/(s11*s22)^0.5 #parámetro de correlación  cov/producto de desvestandar. Da cero porque la covarianza es cero.
b<- rmvnorm(n, mean = mu, sigma = cov, method=c("chol")) # ("eigen", "svd, "chol")
#method método de generación de números aleatorios cholesky
cor(b[,1], b[,2]) #correlación muestral (rho=0)

#Definir una función en R.
#C (cópula Clayton)
#Suponer que ml es un vector normal
m1=mean(b[,1]); m2=mean(b[,2]) #m1 es media de la primer variable, m2 es " de la 2a variable
v1=var(b[,1]); v2=var(b[,2]) #v1=desviación estadar de la primera variable, v2 " 2a
s1=sd(b[,1]); s2=sd(b[,2]) # s1 es desviación estandar de la 1a variable


# Estimación del parámetro de ligadura de la cópula de Clayton

#Estimación del parámetro de max vero
#densidad de una cópula c=C''(u,v) * f(u) f(v); con C cópula
# C(u,v)=(u^-theta) + v^(-theta)-1)^(-1/theta)
# c(u,v)= C''(u,v)=(u^-(tehta) + v^(-theta)-1)^(-1/theta-2) * ((1+theta) * (u*v)^(-theta-1))


fexp<- function(theta){
  u=pnorm(b[,1], mean=m1, sd=s1); v=pnorm(b[,2], mean=m2, sd=s2)
  c=u^(-theta)+v^(-theta)-1 #C mayúscula en apuntes
  if(c<0) c=0 
  c=(-1/theta-2)*log(c)+log(1+theta)-(theta+1)*log(u*v) # c es ln(c); c es c mínuscula en apuntes.
  L=-sum(c) #minimiza el negativo, obteniendo el máximo para theta E.
}
theta=0.5
sol=optim(theta, fexp); sol
thetaE=sol$par; thetaE 

#Función cópula para el cálculo de una probabilidad usando la cópula estimada
#Cálculo de: PR(x1<a1, X2<a2)
copClay=function(a){
  theta=a[3]
  u=pnorm(a[1], mean=m1, sd=s1); v=pnorm(a[2],mean=m2,sd=s2) 
  k=u^(-theta)+ v^(-theta)-1
  if(k<0) k=0
  C=(u^(-theta)+v^(-theta)-1)^(-1/theta)
  return(C)
}
# Cálculo de la probabilidad acumulada en los puntos (0,0) y (1.96,1.96).
#Para la normal estandar la probabilidad en (0,0) es p=0.5^2=0.25 y en 
#(1.96, 1.96) p=0.975^2=0.9506
copClay(c(0,0, thetaE)) 
copClay(c(1.96,1.96, thetaE))

#Cuando no hay independencia de los valores usar:
#Valores no estimados de las probabilidades con base en la normal ultivariada.

pmvnorm(lower=-Inf, upper=c(0,0), mean=mu, sigma=cov)
pmvnorm(lower)
          
 ##Correlación estimada

##Parámetro rho=0
#Rho de Spearman
R=cor(b[,1],b[,2], method="pearson");R

#Tau de Kendall
#El parámetro tau

Tau=function(x){
  F=pmvnorm(lower=-Inf, upper=c(x[1],x[2]), mean=mu, sigma=cov)
  f=dmvnorm(c(x[1],x[2]), mean=mu, sigma=cov)
  I=F*f; return(I)
}
pc=adaptIntegrate(Tau, lowerLimit = c(-10, -10), upperLimit=c(10,10))$integral
tau=4*pc-1;tau

#Estimador no paramétrico del tau
tauk=cor(b[,1], b[,2], method="kendall");tauk

#Estimación del tau de Kendall con base en la cópula estimada
tauCopula=thetaE/(thetaE+2); tauCopula

tauCopulaf=function(x){
  u=x[1]; v=x[2]
  k=u^(-thetaE)+ v^(-thetaE)-1
  if(k<0)k=0
  C=(k)^(-1/thetaE)
  c=(k)^(-1/thetaE-2)* (1+thetaE) * (u*v)^(-thetaE-1)
  I=C*c
  return(I)
}
4*adaptIntegrate(tauCopulaf, lowerLimit = c(0.0,0.0), upperLimit = c(1,1))$integral-1
