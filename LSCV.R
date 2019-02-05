b=c(20.9,18.2,20,17.3,19.6,13.6,24.9,26.9,23.5,21.8,17,20.4,24.6,22.6,21.2,19.6,14.6,24.4,21.8,18.4,24.8,28.5,11.9,10,25.7,27.2,24.4,30.1,21.6,26,14.6,26.1,22.1,8.4,16.4,19.6,19.6,21.5,20.2,25.2,26.7,22.3,22.9,19.9,16.5,14.1,20.4,16.6,19.1,25.5,16.2,24.7,20,28.4,24.4,15.8,25.6,22.5,17.2,15.8,15.1,16.2,19.9,27.3,22.3,19.3,11.7,14.4,24.5,21.6,12.4,15.9,23.5,22.8,26.6,31,22.2,21.7,25.1,28.8,22.8,21.3,24.5,13.8,14.3,23.6,13.3,28.6,22.9,13.7,15.4,13.1,28.8,11.2,22.3,21.9,11.2,21.2,18.7,15)
n=length(b)

kerNorm= function(xx,xi,hh){
  nn=length(xi)
  zz=(xx-xi)/hh
  kk=(1/((2*pi)^0.5))*exp(-0.5*zz*zz)
  fes= sum(kk)/(hh*nn)
  return(fes)
}

LSCV= function(h.LSCV){
  ##Se calcula la integral de f estimado al cuadrado
  points=density(b)$x;dx=points[2]-points[1];dx # determino el valor de dx
  I=0
  for(i in 1:n){
    fxb=kerNorm(b[i],b, h.LSCV)
    I= I+ dx*sum((fxb)^2)  
  } 
  ##Se calcula la sumatoria de f estimado sin la i-esima observación
  LSCVh=0
  ind=1:n
  for(i in 1:n){
    fi=kerNorm(b[i],b[ind != i], h.LSCV) 
    LSCVh=LSCVh+fi
  }
  LSCVh=I-(2/n)*(LSCVh)
  return(LSCVh)
}
sol=optimize(LSCV,2,lower=0,upper=4)
h.LSCV<-sol$minimum;h.LSCV #1.675497
library(sm)
h.LSCVs <- hcv(b,ngrid=32);h.LSCVs #2.286215
