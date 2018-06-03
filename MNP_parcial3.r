### MNP parcial 3 ###

############################################################
### punto 1
############################################################
mes=1:24
ipc=c(1.56,    1.22,    0.47,    0.03,    0.29,    0.35,    0.17,    0.91,    2.21,    1.7,    0.94,    0.78,    0.48,    0.28,    0.31,    0.5,    0.33,    0.35,    0.48,    0.53,    1.29,    2.3,    1.71,    1)
Invext=c(15086.5,    16308.4,    16272.7,    16097.1,    16265.7,    16385.5,    16869.2,    16885,    17129,    17701.9,    18550.5,    18644.7,    18871.3,    19050.1,    19218.8,    19352.4,    19669.9,    20196,    20562.8,    20859.4,    20960.1,    21045.4,    19877.2,    19961.4)
Import=c(2193.2 , 2242.6, 1725.2, 2159.2, 1899, 2015.4, 1646.5, 1736.5, 1412.6, 1297.5, 2330.9, 1343.3, 1398.5, 1552.1, 1260.5, 1273.4, 1177.1, 1074, 1173.4, 1189.3, 1184.9, 1374.5, 1273, 2664.6)



############################################################
### punto 2 
############################################################
cit<-read.table("cittarium.txt", header=TRUE, dec=",")
attach(cit)
altura<-cit$ALTO
longitud<-cit$LARGOLC
library(sm)
#H0: La asociación entre altura y longitud es lineal
#H1: La asociación entre altura y longitud no es lineal
regnp=sm.regression(altura,longitud, model="linear",test=TRUE, col=2, lwd=2)
regnp$p # P-valor de la prueba 
# como el p-valor de la prueba es cero se rechaza H0 lineal.

estacion1<-subset(cit, ESTACION=="1")
estacion2<-subset(cit, ESTACION=="2")
E12<-rbind(estacion1, estacion2)
#H0: la curva de regresión (largo vs ancho) es la misma
#H1: la curva de regresión es diferente
ancova1 <-sm.ancova(E12$LARGOLC, E12$ANCHOCC, E12$ESTACION , model="equal")
# se usa sm.ancova con el argumento model="equal" para probar el sistema hipótesis
#el p-valor de la prueba es 0.0244, luego a un nivel de significacnia del 5% se rechaza H0
# se concluye que la relación largo y ancho de las especies  difieren entre las estaciones uno y dos 

############################################################
### punto 3
############################################################
qxTmor<-read.table("qxTmor.txt", header=TRUE, dec=",")
x<-qxTmor$edad
y<- qxTmor$qx
plot(x,qx, xlim=c(20,100)) #lines(edad,qx)
t1=31;t2=60
tm1=(x-t1);tm2=(x-t2)
tm1[t1 > x] <- 0;tm2[t2 > x] <- 0
tm12=tm1^2 ; tm22 = tm2^2
tm13=tm1^3 ; tm23 = tm2^3
x2=x^2;x3=x^3
xs=seq(min(x), 100, 0.1);#xs=seq(min(x), 5, 0.1);
xs2=xs^2 ; xs3 = xs^3
tm1s=(xs-t1);tm2s=(xs-t2)
tm1s[t1 > xs] <- 0;tm2s[t2 > xs] <- 0
tm1s2 =tm1s^2 ; tm2s2 =tm2s^2 
tm1s3 =tm1s^3 ; tm2s3 =tm2s^3 
lmq = lm( y ~ x + x2 + tm12 + tm22);lmq
b0=lmq$coefficients[1];b1=lmq$coefficients[2];b2=lmq$coefficients[3];b3=lmq$coefficients[4];b4=lmq$coefficients[5]
ypw = b0 + b1 * xs + b2*xs2 + b3 * tm1s2 + b4 * tm2s2
plot(x,y, xlim=c(0,100), ylim=c(0, max(ypw)));
lines(xs,ypw)
# Se usa una regresión por spline cuadrático, en la cual los nodos son las edades 31 y 60
# En el gráfico se ve que la proyeción de qx va a umentando hasta una edad x=100
