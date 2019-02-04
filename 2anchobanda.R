b=c(20.9,18.2,20,17.3,19.6,13.6,24.9,26.9,23.5,21.8,17,20.4,24.6,22.6,21.2,19.6,14.6,24.4,21.8,18.4,24.8,28.5,11.9,10,25.7,27.2,24.4,30.1,21.6,26,14.6,26.1,22.1,8.4,16.4,19.6,19.6,21.5,20.2,25.2,26.7,22.3,22.9,19.9,16.5,14.1,20.4,16.6,19.1,25.5,16.2,24.7,20,28.4,24.4,15.8,25.6,22.5,17.2,15.8,15.1,16.2,19.9,27.3,22.3,19.3,11.7,14.4,24.5,21.6,12.4,15.9,23.5,22.8,26.6,31,22.2,21.7,25.1,28.8,22.8,21.3,24.5,13.8,14.3,23.6,13.3,28.6,22.9,13.7,15.4,13.1,28.8,11.2,22.3,21.9,11.2,21.2,18.7,15)
n=length(b)
## valor del ancho de banda fijos
library(sm)
## h.Normal

DE=sd(b)
h.N=(1.06*DE)/(n^(1/5));h.N

hnorm(b)
 ## hS:Silverman . R basico
RI=(summary(b)[5]-summary(b)[2])/1.349
A= min(DE,RI)
h.S=(0.9*A)/(n^(1/5));h.S
density(b)$bw

#Sobresuavizado
h.osg <- (1.144*DE)/(n^(1/5));h.osg
hnorm(b)*1.08

##Validación cruzada (Mirarlo c:)
par(mfrow=c(1,2))
h.LSCV <- hcv(b,display="lines",ngrid=32);h.LSCV
sm.density(b,h=h.LSCV)
par(mfrow=c(1,1))

##plug-in (No la pide)
h.PI=hsj(b);h.PI

# Validacion cruzada por maxima verosimilitud
kerNorm= function(xx,xi,hh){
  nn=length(xi)
  zz=(xx-xi)/hh
  kk=(1/((2*pi)^0.5))*exp(-0.5*zz*zz)
  fes= sum(kk)/(hh*nn)
  return(fes)
}

LCV= function(h.mv){
  LCVh=0
  ind=1:n
  for(i in 1:n){
    fi=kerNorm(b[i],b[ind != i], h.mv) #quitar el individuo i TAREA: USAR DENSITY PARA CALCULARLO 
    LCVh=LCVh+log(fi)
  }
  LCVv=-(LCVh)
}

h.MVCV=optimize(LCV,h.S,lower=1.5,upper=4)$minimum;h.MVCV
#Si se quita el criterio de validación, toma el valor del límite inferior


## bootstrap
np=500
points=density(b,n=np)$x;dx=points[2]-points[1];dx
fx=density (b,n=np, kernel="gaussian")$y #(La de la derecha)
m=1000
booth=function(h){
  E=0
  for(j in 1:m){
    muestra=b[sample(1:n,replace=T)]
    fxb=density(muestra,bw=h,n=np,kernel="gaussian")$y
    E= E+ dx*sum((fxb-fx)^2)
  }
  MISE=E/m
}
solu=optimize(booth,1,lower=0.5,upper=4);solu
h.boot<-solu$minimum;h.boot
nmise(sd(b),n,h.boot)

np=2000
points=density(b,n=np)$x;dx=points[2]-points[1];dx
fx=density (b,n=np, kernel="gaussian",bw=h.osg)$y #(La de la derecha)
m=1000
booth=function(h){
  E=0
  for(j in 1:m){
    muestra=b[sample(1:n,replace=T)]
    fxb=density(muestra,bw=h,n=np,kernel="gaussian")$y
    E= E+ dx*sum((fxb-fx)^2)
  }
  MISE=E/m
}
solu=optimize(booth,1,lower=0.5,upper=4);solu
h.boot<-solu$minimum;h.boot
nmise(sd(b),n,h.boot)


#RESUMEN
h.boot #Bootstrap
h.osg #Sobresuav.
h.LSCV #Cruzada LS
h.N #Normal
h.PI #plug-in
h.MVCV #Cruzada MV
h.S #Silverman (peque)

##Grafica con h seleccionado
hist(b,prob="T")
lines(density(b,bw= h.boot, kernel = "gaussian"),col=2,lwd=5)
lines(density(b,bw= h.osg, kernel = "gaussian"),col=1,lwd=3)
lines(density(b,bw= h.LSCV, kernel = "gaussian"),col=3)
lines(density(b,bw= h.N, kernel = "gaussian"),col=4)
lines(density(b,bw= h.PI, kernel = "gaussian"),col=5)
lines(density(b,bw= h.MVCV, kernel = "gaussian"),col=6)
lines(density(b, kernel = "gaussian"),col=9)

# ANcho de banda variable
Den=density(b,kernel="gaussian")
fe=Den$y; points=Den$x;h.S=Den$bw
mg=exp(sum(log(fe))/length(fe));mg
alpha=0.5;
hi=h.S*(fe/mg)^(-alpha);hi
fe2=NULL 
for(i in 1:length(points)){
  z=(points[i]-b)/hi[i]
  k=(1/((2*pi)^0.5))*exp(-0.5*z*z)
  fe2[i]=sum(k)/(hi[i]*n)
}
hist(b, prob=TRUE,main=paste("Histograma"),xlab="edad")
lines(points,fe2,col=2,lwd=4)
lines(points,fe,col=3,lwd=2)
