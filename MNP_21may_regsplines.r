##############################################
######### Regresión kernel #################
##############################################

############ Datos #########
mes=1:24
ipc=c(1.56,    1.22,    0.47,    0.03,    0.29,    0.35,    0.17,    0.91,    2.21,    1.7,    0.94,  0.78,    0.48,    0.28,    0.31,    0.5,    0.33,    0.35,    0.48,    0.53,    1.29,    2.3,    1.71,    1)
Invext=c(15086.5,    16308.4,    16272.7,    16097.1,    16265.7,    16385.5,    16869.2,    16885,    17129,    17701.9,    18550.5,    18644.7,    18871.3,    19050.1,    19218.8,    19352.4,    19669.9,    20196,    20562.8,    20859.4,    20960.1,    21045.4,    19877.2,    19961.4)


##############  Reg. simple  #################
#### ipc con mes
y=ipc;x=mes;h=2
points=seq(min(mes), max(mes), 0.1)
a=rep(0,length(points))
for (i in 1:length(points)){
  z=(x-points[i])/h;       k=(1/((2*pi)^0.5))*exp(-0.5*z*z)
  wi=k/sum(k);             xt=x-points[i]
  lm1=lm(y~xt,weights=wi); a[i]=lm1$coefficients[1]
}
plot(mes,ipc);lines(points,a)
# Verificación procedimiento sm
library(sm)
reg1=sm.regression(mes,ipc,h=2)


################### Regresión Bivariada ##################

h1=2.4;h2=800;
pts1=seq(min(mes), max(mes), 0.1)
pts2=seq(min(Invext), max(Invext), length.out=length(pts1))

a=rep(0,length(pts1))
for (i in 1:length(pts1)){
  z=(pts1[i]-mes)/h1
  k1=(1/((2*pi)^0.5))*exp(-0.5*z*z)
  z=(pts2[i]-Invext)/h2
  k2=(1/((2*pi)^0.5))*exp(-0.5*z*z)
  wi=(k1*k2)/sum(k1*k2)
  xmes=mes-pts1[i];xmes2=(mes-pts1[i])^2
  xInvext=Invext-pts2[i];xInvext2=(Invext-pts2[i])^2
  xMesInvest=(mes-pts1[i])*(Invext-pts2[i])
  lm1=lm(ipc~xmes+xInvext,weights=wi)
  #lm1=lm(ipc~xmes+xInvext+xMesInvest,weights=wi)
  #lm1=lm(ipc~xmes+xInvext+xmes2+xInvext2,weights=wi)
  #lm1=lm(ipc~xmes+xInvext+xmes2+xInvext2+xMesInvest,weights=wi)
  a[i]=lm1$coefficients[1]
}
plot(mes,ipc,ylim=c(0,2.4));lines(pts1,a)

#Verificación paquete sm
ye=sm.regression(cbind(mes,Invext),ipc,h=c(h1,h2),eval.points=cbind(pts1,pts2),kernel="gaussian")
plot(mes,ipc)
lines(pts1,diag(ye$estimate))


###################################
##### Regresión por splines  ######
###################################


x=c(1,    1.2,    1.5,    2,    3,    3.7,    4,    4.5)
y=c(3,    3.4,    5,    2,    4.1,    5,    7,    6.5)
plot(x,y);lines(x,y)
t1=2;t2=4
tm1=(x-t1);tm2=(x-t2)
tm1[t1 > x] <- 0;tm2[t2 > x] <- 0
tm12=tm1^2 ; tm22 = tm2^2
tm13=tm1^3 ; tm23 = tm2^3
x2=x^2;x3=x^3

# Definir secuencias para graficar ecuaciones estimadas
xs=seq(min(x), max(x), 0.1);#xs=seq(min(x), 5, 0.1);
xs2=xs^2 ; xs3 = xs^3
tm1s=(xs-t1);tm2s=(xs-t2)
tm1s[t1 > xs] <- 0;tm2s[t2 > xs] <- 0
tm1s2 =tm1s^2 ; tm2s2 =tm2s^2 
tm1s3 =tm1s^3 ; tm2s3 =tm2s^3 

# Piece-wise regression - un nodo
lm.1 = lm( y ~ x + tm1)
b0=lm.1$coefficients[1];b1=lm.1$coefficients[2];b2=lm.1$coefficients[3];
ypw = b0 + b1 * xs + b2 * tm1s
plot(x,y);lines(xs,ypw)

# Piece-wise regression - dos nodos
lm.2 = lm( y ~ x + tm1 + tm2);lm.2
b0=lm.2$coefficients[1];b1=lm.2$coefficients[2];b2=lm.2$coefficients[3];b3=lm.2$coefficients[4]
ypw = b0 + b1 * xs + b2 * tm1s + b3 * tm2s
plot(x,y);lines(xs,ypw)

# Spline cuadrático
lm.3 = lm( y ~ x + x2 + tm12 + tm22);lm.3
b0=lm.3$coefficients[1];b1=lm.3$coefficients[2];b2=lm.3$coefficients[3];b3=lm.3$coefficients[4];b4=lm.3$coefficients[5]
ypw = b0 + b1 * xs + b2*xs2 + b3 * tm1s2 + b4 * tm2s2
plot(x,y,ylim=c(0,7.5));lines(xs,ypw)

# Spline Cúbico
lm.4 = lm( y ~ x + x2  + x3 + tm13 + tm23);lm.4
b0=lm.4$coefficients[1];b1=lm.4$coefficients[2];b2=lm.4$coefficients[3];b3=lm.4$coefficients[4];b4=lm.4$coefficients[5];b5=lm.4$coefficients[6]
ypw = b0 + b1 * xs + b2 * xs2 + b3*xs3 + b4 * tm1s3 + b5 * tm2s3
plot(x,y);lines(xs,ypw)
