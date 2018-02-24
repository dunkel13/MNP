# solución planteada para 10000 muestras de tamaño n
n<-10000
x<-rnorm(10000, mean=20, sd=5)
# espacio muestral: una población normal de parámetros mu=20 y sigma:25
m<-matrix(nrow=n, ncol=50, NA)
for (k in 1:10000){
    m[k,]<-sample(x,50)
    }
View(m)
#calculamos la media muestral de cada muestra, se la asignamos a una componente de vector de medias vm
vm<-matrix(nrow=10000, ncol=1, NA)
for(i in 1:10000){ 
    vm[i,]<-mean(m[i,])
    } 
View(vm)
# prueba de hipótesis H0: mu=20 vs H1:mu>20
alpha<-0.05
mu=20
sigma=25
Z<-matrix(nrow=10000, ncol=1, NA)
for(j in 1:10000) {   
    Z[j,]<-((vm[j,]-20)/sqrt(sigma))
    } 
zalpha<- qnorm(alpha, mean=0, sd=1, lower.tail = FALSE)
zalpha
# Se quiere contar el número de muestras en las que Z>zalpha, es decir, en las que se rechaza H0.
c<-matrix(nrow=10000, ncol=1, NA)
for (i in 1:10000){
  if(Z[i,]>=zalpha){
    c[i,]<-1
  } else c[i,]<-0
}
sum(c)
