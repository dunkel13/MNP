###############################################
#### ESTIMADOR KERNELL DE LA DENSIDAD #########
###############################################
#### Generación de la muestra
b=c(20.9,18.2,20,17.3,19.6,13.6,24.9,26.9,23.5,21.8,17,20.4,24.6,22.6,21.2,19.6,14.6,24.4,21.8,18.4,24.8,28.5,11.9,10,25.7,27.2,24.4,30.1,21.6,26,14.6,26.1,22.1,8.4,16.4,19.6,19.6,21.5,20.2,25.2,26.7,22.3,22.9,19.9,16.5,14.1,20.4,16.6,19.1,25.5,16.2,24.7,20,28.4,24.4,15.8,25.6,22.5,17.2,15.8,15.1,16.2,19.9,27.3,22.3,19.3,11.7,14.4,24.5,21.6,12.4,15.9,23.5,22.8,26.6,31,22.2,21.7,25.1,28.8,22.8,21.3,24.5,13.8,14.3,23.6,13.3,28.6,22.9,13.7,15.4,13.1,28.8,11.2,22.3,21.9,11.2,21.2,18.7,15)
n=length(b)

### Procedimiento del R básico
hist(b, prob=TRUE)
al=density(b,kernel = "epanechnikov")
# estimación de la linea del gráfico, hay varias opciones en el kernell
# f_{k}(x)= \frac{1}{n}\sum{i=1}_^n k \frac{x_{i}-x}{h}
# unif k(u) = \frac{1}{2} I_[-1,1](u)
# epanech k(u) = \frac{3}{4} (1-u^2) I_[-1,1](u)
lines(al)
density(b,bw=0.16)$x #valores sobre los que se calcula la densidad
density(b,bw=0.16)$y # valores de la densidad
# valores de f estimado
max(b)

## código detallado
## valor del ancho de banda
DE=sd(b)
RI=(summary(b)[5]-summary(b)[2])/1.35
# el 1.35 sale de minimizar el error cuadratico medio del estimador \hat{f_{k}} para h
minimo=min(DE, RI)
h=(0.9*minimo)/(n^(1/5));h
al$bw
# Selección de puntos donde se estimará la densidad 
points0=c(5,8,10,12,14,16,18,20,22,24,26,28,30,32,40)
# denisdad con el R básico evaluada en los puntos seleccionados
hist(b, breaks=points0, prob=TRUE, main=paste("Histograma "), xlab="x")
# breaks marca de clase, prob= TRUE frecuencias relativas
lines(density(b,bw=0.5, kernel="epanechnikov"))
# ancho de banda pequeño no suaviza, ancho de banda grande suaviza mucho

# otros puntos9
points1=seq(min(b), max(b), 0.01); points1
#anchos de banda muy pequeñitos 
hist(b, breaks=points1, prob=TRUE, main=paste("Hist "), xlab="x")
lines(density(b, kernell="epanechnikov"))

# Ahora si pasemos al código :v

h=2.5
## ciclo para el cálculo de la densidad. Uniforme
points=c(7.5,12.5,17.5,22.5, 27.5, 32.5)
fe=rep(0,length(points))
for (i in 1:length(points)){
  # z es (x_{i}-x)/h que es igual a u
  z=(points[i]-b)/h
  k=rep(0.5, length(z)); for(j in 1:length(z)){
                                if(z[j]<(-1) | z[j]>=1) k[j]=0
  }
  fe[i]=sum(k)/(h*n)
}
fe
hist(b, prob=TRUE, main=paste("Histograma "), xlab="edad")$density
lines(points,fe, col=2, lwd=4)

## ciclo para el cálculo de la densidad. Epanech
# points=seq(min(b)-h/2, max(b)+h/2, by=0.05 );points
h=2 #h=4
# h ancho de banda
points= seq(min(b)-2*h, max(b)+2*h, length.out = 50);points 
fe=rep(0,length(points))
for (i in 1:length(points)){
  # z es (x_{i}-x)/h que es igual a u
  z=(points[i]-b)/h
  k=(1/((2*pi)^0.5))*exp(-0.5*z*z)
  fe[i]=sum(k)/(h*n)
}
fe
hist(b, prob=TRUE, main=paste("Histograma "), xlab="edad")$density
lines(points,fe, col=2, lwd=4)
lines(density(b, bw=2, kernel="gaussian"))
# Ciclo para el cálculo de la densidad. Uniforme
h=2.5 #h=4
# h ancho de banda
points= seq(min(b)-2*h, max(b)+2*h, length.out = 1000);points 
fe=rep(0,length(points))
for (i in 1:length(points)){
  # z es (x_{i}-x)/h que es igual a u
  z=(points[i]-b)/h
  k=rep(0.5, length(z)); for(j in 1:length(z)){
    if(z[j]<(-1) | z[j]>=1) k[j]=0
  }
  fe[i]=sum(k)/(h*n)
}
fe
hist(b, prob=TRUE, main=paste("Histograma "), xlab="edad")$density
lines(points,fe, col=2, lwd=4)
lines(density(b, bw=2, kernel="gaussian"))
# Ciclo para el cálculo de la densidad. Epanechnikov
h=4.5 #h=4
# h ancho de banda
#points= seq(min(b)-2*h, max(b)+2*h, length.out = 1000);points 
fe=rep(0,length(points))
for (i in 1:length(points)){
  # z es (x_{i}-x)/h que es igual a u
  z=(points[i]-b)/h
  k= 3/4*(1-z*z) 
  for(j in 1:length(z)){
    if(z[j]<(-1) | z[j]>=1) k[j]=0
  }
  fe[i]=sum(k)/(h*n)
}
fe
hist(b, prob=TRUE, main=paste("Histograma "), xlab="edad")$density
lines(points,fe, col=2, lwd=4)
lines(density(b, bw=2, kernel="epanechnikov"))

hist(b, prob=TRUE, main=paste("Histograma "), xlab="edad")$density
lines(density(b, bw=2, kernel="epanechnikov"), col=5)
lines(density(b, bw=2, kernel="gaussian"), col=4)
lines(density(b, bw=2, kernel="biweight"), col=3)
lines(density(b, bw=2, kernel="cosine"), col=2)
lines(density(b, bw=2, kernel="optcosine"), col=1)

############################
########## Mixtura #########
############################
## Generar la muestra, estimador MV y estimador kernell, dan lo mismo
# para MV debe empezar a asumir dist
m=10000
x1=rnorm(m, mean=0, sd=1)
x2=rnorm(m, mean=3/2, sd=1/3)

hist(x1, prob=TRUE)
lines(density(x1), col=2)
hist(x2, prob=TRUE)
lines(density(x2), col=2)

b=x1
for(i in 1: length(b)){
  c=runif(1, min=0, max=1)
  if(c>=0.75) b[i]=x2[i]
}
hist(x2, prob=TRUE)
lines(density(x2), col=2)
lines(density(x1), col=2)
