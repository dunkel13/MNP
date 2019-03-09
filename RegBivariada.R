rm.list=ls()

#----------------------
# Lectura
#----------------------

# Datos referido a la vivienda en 506 barrios de Boston en 1978
# AGE: Proporción de unidades construidas antes de 1940 ocupadas por el propietario
# lstat: Menor estatus de la población
# rm: Número medio de habitaciones por vivienda

datos<-read.table("lstat.txt",head=TRUE)
head(datos)
attach(datos)

pairs(~lstat+rm+age, data=datos, col=2)



#----------------------
# Regresión con sm
#----------------------

library(sm)
covariables<-cbind(age,lstat)
h1=10;h2=3;h<-c(h1,h2)
modelo<-sm.regression(covariables,rm,h=h,design.mat=NA,model="none")


#----------------------
#con codigo
#----------------------

y<-rm
x1<-age
x2<-lstat

points1=seq(min(x1), max(x1), 5);points1
points2=seq(min(x2), max(x2), 2);points2


a=matrix(0,length(points1),length(points2))
for(i in 1:length(points1)){
     z1=(x1-points1[i])/h1
     k1=(1/((2*pi)^0.5))*exp(-0.5*z1*z1)   
     xt1=x1-points1[i]
     #xt12=(x1-points1[i])^2
     #xt13=(x1-points1[i])^3

  for(j in 1:length(points2)){
     z2=(x2-points2[j])/h2
     k2=(1/((2*pi)^0.5))*exp(-0.5*z2*z2)   
     xt2=x2-points2[j]
     #xt22=(x2-points2[j])^2
     #xt23=(x2-points2[j])^3
     wi=k1*k2
     lm1=lm(y~xt1+xt2+xt1*xt2,weights=wi) #ajusta un polinomio grado 1
     a[i,j]=lm1$coefficients[1]
     }
}
a

persp(points1,points2,a, phi=10, theta=-60, xlim=c(0,100), ylim=c(4,40),ticktype="detailed")

plot(lstat,rm)
lines(points2,diag(a),col=2,lwd=2)
