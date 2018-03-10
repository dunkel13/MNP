###########################################################################
###########################################################################
###########################################################################

remove(list=ls())

base=read.table("DatoSimulado.txt", header=T)
x=base$x

#### Estimación de parámetros de diferentes distribuciones con base en la muestra


### Weibull
fWei<-function(d){
	a=d[1];b=d[2]
	fw=(a/b)*(x/b)^(a-1)*exp(- (x/b)^a)
      L=-sum(log(fw))
}
d=c(1,1)
sol=optim(d,fWei);sol
a=sol$par[1];b=sol$par[2]

###### Exponencial
fExp<-function(lambda){
	fEx=lambda*exp(- lambda*x)
      L=-sum(log(fEx))
}
solEx=optim(1,fExp);solEx
lambda=solEx$par[1];lambda

###### Normal
fNorm<-function(d){
	mu=d[1];sigma=d[2]
	fN=1/((2*pi)^0.5*sigma)*exp(-0.5*((x-mu)/sigma)^2)
      L=-sum(log(fN))
}
solN=optim(c(1,1),fNorm);solN
mu=solN$par[1];sigma=solN$par[2];

###### Hi cuadrado
fChi<-function(n){
	fc = 1 / (2^(n/2)* gamma(n/2))* x^(n/2-1)* exp(-x/2)
      L=-sum(log(fc))
}
solChi=optim(3,fChi);solChi
N=solChi$par[1]

###### Gamma
#E(X) = a*s ; Var(X) = a*s^2
fGamma<-function(d){
	a=d[1];s=d[2]
	fg= 1/(s^a * gamma(a)) * x^(a-1) * exp(-(x/s))
      L=-sum(log(fg))
}
solGamma=optim(c(1,1),fGamma);solGamma
A=solGamma$par[1];s=solGamma$par[2]

###### Lognormal
#E(X) = exp(µ + 1/2 s^2)  ;  Var(X) = exp(2*µ + s^2)*(exp(s^2) - 1)
fLogN<-function(d){
	mu=d[1];sigma=d[2]
	fln=1/((2*pi)^0.5 * sigma* x)* exp( -((log(x) - mu)^2 / (2 * sigma^2)) )
      L=-sum(log(fln))
}
solLogN=optim(c(0.5,0.5),fLogN);solLogN
muL=solLogN$par[1];sigmaL=solLogN$par[2]


###### Gráfica de las Estimaciones #########

pts=seq(0,3,0.01)								## Definir recorrido de la muestra
fw=(a/b)*(pts/b)^(a-1)*exp(- (pts/b)^a)   			## Weibull
fEx=lambda*exp(- lambda*pts)						## Exponencial
fNor=1/((2*pi)^0.5*sigma)*exp(-0.5*((pts-mu)/sigma)^2)	## Normal
fChi = 1 / (2^(N/2)* gamma(N/2))* pts^(N/2-1)* exp(-pts/2)	## Chi
fGamma= 1/(s^A * gamma(A)) * pts^(A-1) * exp(-(pts/s))	## Gamma
fLogN=1/((2*pi)^0.5 * sigmaL* pts)* exp( -((log(pts) - muL)^2 / (2 * sigmaL^2)) )	## LogNormal

hist(x,prob=TRUE,xlim=c(0,3))
lines(pts,fEx,col=2,lwd=2)
lines(pts,fChi,col=4,lwd=2)

hist(x,prob=TRUE,xlim=c(0,2))
lines(pts,fw,col=3,lwd=2)
lines(pts,fNor,col=1,lwd=2)
lines(pts,fGamma,col=5,lwd=2)
lines(pts,fLogN,col=6,lwd=2)


#### Pruebas KS 
ks.test(x,"pchisq",N)
ks.test(x,"pexp",lambda)

ks.test(x,"pgamma",A,s)
ks.test(x,"plnorm",muL,sigmaL)
ks.test(x,"pweibull",a,b)
ks.test(x,"pnorm",mu,sigma)







###########################################################################
###########################################################################
###########################################################################



## Weibull trasladada

# Leer base
base=read.table("DatoSimulado.txt", header=T)
x=base$x

# Estimar datos con base en una Weibull
### Weibull
fWei<-function(d){
	a=d[1];b=d[2]
	fw=(a/b)*(x/b)^(a-1)*exp(- (x/b)^a)
      L=-sum(log(fw))
}
d=c(1,1)
sol=optim(d,fWei);sol
a=sol$par[1];b=sol$par[2]


### Weibull con parámetro de localización

#Desplazar la distribución
y=x+15

fWei<-function(d){
	a2=d[1];b2=d[2];c2=d[3]
	fw=(a2/b2)*((y-c2)/b2)^(a2-1)*exp(- ((y-c2)/b2)^a2)
      L=-sum(log(fw))
}
d=c(1,1,1)
sol=optim(d,fWei);sol
a2=sol$par[1];b2=sol$par[2];c2=sol$par[3];


### Weibull sin parámetro de localización de la distribución desplazada
fWei<-function(d){
	a3=d[1];b3=d[2]
	fw=(a3/b3)*(y/b3)^(a3-1)*exp(- (y/b3)^a3)
      L=-sum(log(fw))
}
d=c(1,1)
sol=optim(d,fWei);sol
a3=sol$par[1];b3=sol$par[2]

pts=seq(0,3,0.01)
pts2=seq(15,18,0.01)

fw=(a/b)*(pts/b)^(a-1)*exp(- (pts/b)^a)   	
fw2=(a2/b2)*((pts2-c2)/b2)^(a2-1)*exp(- ((pts2-c2)/b2)^a2)
fw3=(a3/b3)*((pts2)/b3)^(a3-1)*exp(- ((pts2)/b3)^a3)

hist(y,prob=TRUE)
lines(pts2,fw,col=2,lwd=2)
lines(pts2,fw2,col=3,lwd=3)
lines(pts2,fw3,col=4,lwd=4)





###########################################################################
###########################################################################
###########################################################################


################################  MIXTURAS ##########################

m=10000
# generar distribuciones
x1=rnorm(m,mean=0,sd=1)
x2=rnorm(m,mean=3/2,sd=1/3)

par(mfrow=c(2,1))
hist(x1,prob=TRUE)
lines(density(x1),col=2)

hist(x2,prob=TRUE)
lines(density(x2),col=2)

# general la mixtura
# f=(c*f1+(1-c)*f2)
b=x1
for (i in 1:length(b)){
   c=runif(1, min=0, max=1)
   if(c >=0.75) b[i]=x2[i]
}

#Graficar las distribuciones
points=seq(-4.5,4.5,0.2)
hist(b,breaks = points,prob=TRUE)
lines(density(b),lwd=3)
lines(density(x2),col=2)
lines(density(x1),col=2)


# Estimar una mixtura (MV) 

fexp<-function(d){
      L=0;mu1=d[1];s1=d[2];mu2=d[3];s2=d[4];c=d[5]
      for(i in 1:length(b)){
          f1=1/((s1*2*pi)^0.5)*exp(-0.5*(b[i]-mu1)^2/s1 )
          f2=1/((s2*2*pi)^0.5)*exp(-0.5*(b[i]-mu2)^2/s2 )
          L=L+log(c*f1+(1-c)*f2)
      }
      L=-L
}

d=c(0,1,1,1,0.7)
sol=optim(d,fexp);sol
pare=sol$par
#nlminb(d,fexp, lower = c(-3,0,-3,0,0), upper = c(3,5,3,5,1))


# Calcular f con los parámetros estimados
mu1=pare[1];s1=pare[2];mu2=pare[3];s2=pare[4];c=pare[5]
pts=seq(-4.5,4.5,0.05)
fmv=rep(0,length(pts))
for(i in 1:length(pts)){
          f1=1/((s1*2*pi)^0.5)*exp(-0.5*(pts[i]-mu1)^2/s1 )
          f2=1/((s2*2*pi)^0.5)*exp(-0.5*(pts[i]-mu2)^2/s2 )
          fmv[i]=(c*f1+(1-c)*f2)
}

par(mfrow=c(2,1))
hist(b,breaks = points,prob=TRUE)
lines(density(b))
plot(pts,fmv)

#Graficar la mixtura
par(mfrow=c(1,1))
hist(b,breaks = points,prob=TRUE)
lines(density(b),col=2)# Estimación No-Paramétrica del R
lines(pts,fmv,lwd=3)   # Estimación MV
