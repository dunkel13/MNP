################################  MIXTURAS ##########################
# generar distribuciones
{
  m=10000
  x1=rnorm(m,mean=0,sd=1)
  x2=rnorm(m,mean=3/2,sd=1/3)
}
# genera la mixtura
{
  # f=(c*f1+(1-c)*f2)
  b=x1
  for (i in 1:length(b)){
    c=runif(1, min=0, max=1)
    if(c >=0.75) b[i]=x2[i]
  } 
}
# Estimar una mixtura (MV) 
{
  fexp<-function(d){
    L=0;mu1=d[1];s1=d[2];mu2=d[3];s2=d[4];c=d[5]
    for(i in 1:length(b)){
      f1=1/((s1*2*pi)^0.5)*exp(-0.5*(b[i]-mu1)^2/s1 )
      f2=1/((s2*2*pi)^0.5)*exp(-0.5*(b[i]-mu2)^2/s2 )
      L=L+log(c*f1+(1-c)*f2)
    }
    L=-L
  }
  d=c(0,1,1,1,0.7) # valor inicial
  sol=optim(d,fexp);sol
  pare=sol$par
  #nlminb(d,fexp, lower = c(-3,0,-3,0,0), upper = c(3,5,3,5,1))
}
# Calcular f con los parámetros estimados
{
  mu1=pare[1];s1=pare[2];mu2=pare[3];s2=pare[4];c=pare[5]
  pts=seq(-4.5,4.5,0.05)
  fmv=rep(0,length(pts))
  for(i in 1:length(pts)){
    f1=1/((s1*2*pi)^0.5)*exp(-0.5*(pts[i]-mu1)^2/s1 )
    f2=1/((s2*2*pi)^0.5)*exp(-0.5*(pts[i]-mu2)^2/s2 )
    fmv[i]=(c*f1+(1-c)*f2)
  }
}
#Graficar la mixtura
{
  points=seq(-4.5,4.5,0.2)
  par(mfrow=c(2,1))
  hist(b,breaks = points,prob=TRUE)
  lines(density(b))
  plot(pts,fmv)
  # gráficos sobrepuestos
  par(mfrow=c(1,1))
  hist(b,breaks = points,prob=TRUE)
  lines(density(b),col=2)# Estimación No-Paramétrica del R
  lines(pts,fmv,lwd=3)   # Estimación MV
}
