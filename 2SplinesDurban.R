remove(list=ls())

library(splines)
# Simular valores de X y Y
n = 200
x = seq(0,1,length=n)
y = sin(3*pi*x) + 0.5*rnorm(n)
plot(x,y,pch=1,bty="l",col=12)
lines(x,sin(3*pi*x),col=6)

#####################################
## Polinomios truncados de grado p ##
#####################################
tpoly=function(x,t,p){
    B=NULL
    for(i in 1:length(t)){
        B=cbind(B,(x-t[i])^p * (x>t[i]))
    }
    return(B)
}
K=5 # # de nodos
knots=seq(0,1,length=(K+2))[-c(1,K+2)]
#j=1:k+1;#knots=(j-1)/K;

B0=tpoly(x,knots,0)
B1=tpoly(x,knots,1)
B2=tpoly(x,knots,2)
B3=tpoly(x,knots,3)

# Gráfica de los polinomios truncados
par(mfrow=c(2,2))
plot(x,B0[,1],type="n",ylim=c(0,1),ylab="")
title("Polinomios truncados de grado 0",cex.main=1)
for(i in 1:K){lines(x,B0[,i],col=i,lty=i)}
plot(x,B1[,1],type="n",ylim=c(0,1),ylab="")
title("Polinomios truncados de grado 1",cex.main=1)
for(i in 1:K){lines(x,B1[,i],col=i,lty=i)}
plot(x,B2[,1],type="n",ylim=c(0,1),ylab="")
title("Polinomios truncados de grado 2",cex.main=1)
for(i in 1:K){lines(x,B2[,i],col=i,lty=i)}
plot(x,B3[,1],type="n",ylim=c(0,1),ylab="")
title("Polinomios truncados de grado 3",cex.main=1)
for(i in 1:K){lines(x,B3[,i],col=i,lty=i)}
par(mfrow=c(1,1))

# Gráfica de los polinomios truncados grado 3
B3=tpoly(x,knots,3)
B3=cbind(rep(1,n),x,x^2,x^3,B3)
plot(x,B3[,1],type="n",ylim=c(0,1),ylab="")
title("Polinomios truncados de grado 3",cex.main=1)
for(i in 1:(K+4)){lines(x,B3[,i],col=i,lty=i)}

# Regresión spline cúbico
plot(x,y,pch=".",cex=3)
B3=tpoly(x,knots,3);B3=cbind(rep(1,n),x,x^2,x^3,B3)
at.1=solve(t(B3)%*%B3)%*%t(B3)%*%y
lines(x,B3%*%at.1,col=6,lwd=2)
lines(x,sin(3*pi*x),col=4)



################
## B Splines  ##
################

bspline = function(x, xl, xr, ndx, bdeg){
    dx = (xr-xl)/ndx
    knots = seq(xl-bdeg*dx, xr+bdeg*dx, by=dx)
    B = spline.des(knots,x,bdeg+1,0*x)$design
    B
}


xl=-0.0000001
xr=1.0000001
bdeg=3
ndx=8
B=bspline(x,xl,xr,ndx,bdeg)
BB=bspline(x,xl,xr,ndx=9,1)
par(mfrow=c(2,2))
plot(x,BB[,6],type="l",col=12,ylab="",ylim=c(0,1))
points(x[c(89,112,134)],BB[c(89,112,134),6],col=12,pch=19)
points(x[c(89,112,134)],rep(0,3),col=12,pch=19)
plot(x,BB[,1],type="n",ylim=c(0,1),ylab="")
for(i in 1:10){lines(x,BB[,i],col=i,lty=i)}
plot(x,B[,6],type="l",col=12,ylab="",ylim=c(0,1))
points(x[c(51,76,100,125,150)],B[c(51,76,100,125,150),6],col=12,pch=19)
points(x[c(51,76,100,125,150)],rep(0,5),col=12,pch=19)
plot(x,B[,1],type="n",ylim=c(0,1),ylab="")
for(i in 1:10){lines(x,B[,i],col=i,lty=i)}
par(mfrow=c(1,1))
