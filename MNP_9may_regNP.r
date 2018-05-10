### Simulacióm MCO ###
n=100
beta0=10
beta1=1
m=5000 #Num replicas
# consideramos siempre en una simulación las x's como fijas
x<-1:n
bet0=NULL
bet1=NULL
best=NULL
for(j in 1:m ){
  e=rnorm(n, mean=0, sd=2) # errores iid normalmente de varianza constante
  y=beta0+beta1*x+e #muestra
  lm1=lm(y~x) #ejecutamos con el comando lm un MRNL
  bet0[j]=lm1$coefficients[1]
  bet1[j]=lm1$coefficients[2]
}
best[1]<-mean(bet0)
best[2]<-mean(bet1)
est<- data.frame(b0=c(round(best[1],3)), b1=c(round(best[2],3)), row.names="estimaciones");est
cbind(best[1],best[2])
### Regresión por mínimos cuadrados ###
x=c(1,1.2,1.5,2,3,3.7,4,4.5)
y=c(3,3.4,5,2,4.1,5,7,6.5)
plot(x,y);lines(x,y)
w=c(0.33,0.33,0.28,0.06,0,0,0,0)
lm(y~x)
# usando el criterio general de un modelo múltiple sin intercepto donde las variables se han ponderado por sqrt(wi)
wi=sqrt(w)
yw1=wi*y
xw1=wi*x
lm(yw1~wi+xw1-1) #lm(y~x, weights=w)
### Ajuste polinomial ###
x2=x^2; x3=x^3; x4=x^4; x5=x^5; x6=x^6;
lm1=lm(y~x+x2+x3+x4+x5+x6)
X=cbind(rep(1,8), x,x2,x3,x4,x5,x6)
ye=X%*%lm1$coefficients
plot(x,y); lines(x,y); lines(x, ye, col=2, lwd=3)
# clauclo de y estimado cuando x es igual a 5 (predicción)
xn=5
xf=c(1,xn,xn^2,xn^3,xn^4,xn^5,xn^6)
Xf=rbind(X,xf) # x futuro, linea de info con la cual se estima el vallor de y cuando x es igual a 5
yf=Xf%*%lm1$coefficients #y estimado
plot(x,y, xlim=c(1,5.5), ylim=c(-15,8)); lines(x,y)
lines(c(x,xn), yf, col=2, lwd=3)
### Modelo no paramétrico por polinomios locales ###
# Con la siguiente instrucción se obtiene el estimador de Nardaya-Watson, este no admite polinomios luego sus suavizamiento no es el mejor
hs=1
xs=seq(min(x),5,0.05)
ys=ksmooth(x,y,kernel="normal", bandwidth=hs, x.points=xs)
ys$x #puntos de evaluación
ys$y # estimaciones
plot(x,y, xlim=c(1,5)); lines(x,y)
lines(ys$x, ys$y, col=2)
# Ksmooth es suavizador con base en el estimador Local mean(Naradaya) para cada punto un modelo de regresión polinómico centrado en x
### pendientes ###
a=12; b=11
ys$x[b]
(ys$y[b]-ys$y[a])/(ys$x[b]-ys$x[a])
a=36; b=35
(ys$y[b]-ys$y[a])/(ys$x[b]-ys$x[a])
# ancho de banda para que coincida sm.regression con el siguiente código
h=0.3
#h=0.4
#h=0.5
points=seq(min(x), 5,0.05);
a=rep(0, length(points))
b=rep(0, length(points))
for(i in 1:length(points)){
  z=(x-points[i])/h
  k=(1/((2*pi)^0.5))*exp(-0.5*z*z)
  wi=k/sum(k)
  xt=x-points[i] #estimador no paramétrico, polinomio local grado 1
  xt2=(x-points[i])^2
  xt3=(x-points[i])^3
  lm1=lm(y~xt, weights=wi) 
  # lm1=lm(y~xt+xt2,weights=wi)
  # lm1=lm(y~xt+xt2+xt3,weights=wi)
  a[i]=lm1$coefficients[1]
}
points #points evaluados
a # estimaciones
plot(x,y, xlim=c(1,5));
lines(points, a, lwd=3) #linear local
lines(ys$x, ys$y, col=2) # estimador Nadaraya-Watson
library(sm)
h.sm=h
ye=sm.regression(x,y, h=h.sm, model="none", eval.points=points, xlim=c(1,5))
lines(points, a, lwd=3)
lines(x,y)
# ye coincide con la estimación de polinomio local de grado 1
rm(list=ls())
