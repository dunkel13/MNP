### MNP - Método de validación cruzada ###

x<- rnorm(100, mean=20, sd=1)
n=length(x)
h=2

# Función para el kernel gaussiano   
kernelf<- function(mi,h){
  points= seq(min(mi)-2*h, max(mi)+2*h, length.out = 99);points 
  fe=rep(0,length(points))
  for (i in 1:length(points)){
    z=(points[i]-mi)/h
    k=(1/((2*pi)^0.5))*exp(-0.5*z*z)
    fe[i]=sum(k)/(h*(n-1))
  }
  return(fe)
}

# función para hallar el valor de h que hace máxima LCV o la productoria de las f. de densidad estimadas. 
m=NULL

f<- function(h){
  fest=rep(0,length(x))
  for(i in 1:length(x)){
    #asigna a m la estimación kernel eliminando de la muestra el i-ésimo elemento
    m[i]<-x[-i]
    #calcula la estimación de la densidad para cada muestra sin el i-ésimo elemento 
    fest[i]=kernelf(m[i],h)
    print(fest)
  }
  #verosimilitud de la muestra por validación cruzada
  LCV= prod(fest)
  return(LCV)
}
# se halla el máximo con optimise

optimise(f, c(0,4), maximum = TRUE)
rm(list=ls())
