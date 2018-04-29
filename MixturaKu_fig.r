### Mixtura normales estimador kernel de densidad uniforme ###
x<-seq(from=-3, to=3, by=0.01)
fdp<- function(x){ y=0.75*dnorm(x,0,1)+0.25*dnorm(x,1.5,1/3)}
curve(fdp(x), from=-3, to=3)

x1<-rnorm(100, mean=0, sd=1)
x2<-rnorm(100, mean=1.5, sd=1/3)
x<-sample(c(x1,x2),100)
dx<-density(x, kernel="rectangular")
lines(dx)
