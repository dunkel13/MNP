library(sm)
n=100
x1=rnorm(n, mean=0, sd=1)
x2=rnorm(n, mean=3/2, sd=1/3)

h.x=hnorm(x1)
h.y=hnorm(x2)
kerNorm= function(xx,xi,hh){
  nn=length(xi)
  zz=(xx-xi)/hh
  kk=(1/((2*pi)^0.5))*exp(-0.5*zz*zz)
  fes= sum(kk)/(hh*nn)
  return(fes)
}

f.z=function(x1,x2){
  z<-matrix(NA, nr=100, nc=100)
  sum.doble<-matrix(NA, nr=100, nc=100)
  f.est<-NULL
  for(i in 1:100){
    sum.d=NULL
    f.est[i]=kerNorm(x1[i],x1, h.x)
    for(j in 1:100){
      sum.d[j]=f.est[i]*kerNorm(x2[j],x2, h.y)
    }
    sum.doble[i,]<-sum.d
  }
  z<-(1/(n*h.x*h.y))*sum.doble
  return(z)
}
x<-sort(x1)
y<-sort(x2)
persp(x,y,f.z(x,y), phi = 30, expand = 0.5, col = "lightblue")
