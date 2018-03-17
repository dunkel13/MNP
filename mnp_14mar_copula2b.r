library(MASS);library(mvtnorm);library(cubature)
rm(list=ls())

########################################################################
####################### ANALISIS DE CORRELACIÓN  #######################
########################################################################

# Generación de una muestra normal multivariada: media "mu" y varianzas y cov "Cov"

Sigma11 = 1; Sigma22 = 1; Sigma12 = 0.8 ; Mu1=0 ; Mu2 = 0	# Fijar Parámetros
Cov=cbind(c(Sigma11, Sigma12),c(Sigma12, Sigma22));		# Matriz de var-cov
Mu<-c(Mu1,Mu2); 								# Vector de medias
Lx1=5*Sigma11;Lx2=5*Sigma22						# para definir rango de X1 y X2

n=100
b<-rmvnorm(n, mean = Mu, sigma = Cov, method=c("chol"))

#######################################################################
##### Estimación del parámetro de ligadura de la cópula de Clayton  ###
#######################################################################

### A) Estimación del parámetro de ligadura por Máxima Verosimilitud
# Densidad de una cópula c = C''(u,v) * f(u) f(v) ; con C la cópula
# C(u,v)=( u^(-theta) + v^(-theta) -1 )^(-1/theta)
# c(u,v) = C''(u,v) = ( u^(-theta) + v^(-theta) -1 )^(-1/theta - 2) * (1+theta) * (u*v)^(-theta-1)

# Estimación de marginales asumiendo normalidad
m1=mean(b[,1]);m2=mean(b[,2])
s1=sd(b[,1]);s2=sd(b[,2])

# loglikelihood de la copula
fexp<-function(theta){
	u=pnorm(b[,1],mean=m1,sd=s1);v=pnorm(b[,2],mean=m2,sd=s2)
	# ud=dnorm(b[,1],mean=m1,sd=s1);vd=dnorm(b[,2],mean=m2,sd=s2) #densidades
	k= u^(-theta) + v^(-theta) -1
	k[k<0] <- 0
	c= (-1/theta - 2)*log( k ) + log(1+theta) -(theta+1)*log(u*v)
	#c=c + log(ud) + log(vd) # log de la densidad. Genra el mismo estimador pues theta no depende de las densidades pues ud y vd.
	L=-sum(c)
}
# valor de arranque
theta=0.5
# Uso de optimizadores
thetaE=optim(theta,fexp)$par;thetaE
thetaE=nlminb(theta,fexp,lower = 0, upper = 3)$par;thetaE
thetaE=optimize(fexp, 0.2, lower = 0 , upper = 3)$minimum;thetaE


### B) Uso de la cópula "C" para cálculo de probabilidades: Pr( X1<a1 , X2<a2 )
# Función cópula para el cálculo de una probabilidad
CopClay=function(a,theta){
	u = pnorm(a[1],mean=m1,sd=s1); v = pnorm(a[2],mean=m2,sd=s2)
	k=u^(-theta) + v ^(-theta) -1
	if(k<0) k=0
	C=( u^(-theta) + v ^(-theta) -1 )^(-1/theta)
	return(C)
}

# Cálcuo de la probabilidad acumulada en los puntos (0,0) y (1.96 , 1.96).
# Para la normal estándard la probabilidad en (0.0) es p=0.5^2= 0.25 y en 
# (1.96 , 1.96) p = 0.975^2 = 0.9506
CopClay( c(0,0) , thetaE )
CopClay( c(1.96 , 1.96) , thetaE  )

# Valores de las probabilidades con base en la normal multivariada
pmvnorm(lower=-Inf,upper=c(0,0), mean = Mu, sigma = Cov)[1]
pmvnorm(lower=-Inf,upper=c(1.96,1.96), mean = Mu, sigma = Cov)[1]

# Función para el cálculo de la densidad de la cópula "c"
densiCop=function(x){
		u=pnorm(x[1],mean=m1,sd=s1);v=pnorm(x[2],mean=m2,sd=s2)
		ud=dnorm(x[1],mean=m1,sd=s1);vd=dnorm(x[2],mean=m2,sd=s2)
		k= u^(-thetaE) + v^(-thetaE) -1
		if(k<0) k=0
		c= ( k )^(-1/thetaE - 2) * (1+thetaE) * (u*v)^(-thetaE-1)
		c=c*ud*vd
		return(c)
}
# Evaluación de una probabilidad acumulada usando la cópula
adaptIntegrate(densiCop, lowerLimit = c(-Lx1,-Lx2), upperLimit = c(1.96, 1.96))$integral


#######################################################
#############  Correlación Estimada ###################
#######################################################

### A) Correlación de Pearson
# A.1) El Parámetro Rho = 0
Rho = Sigma12 / (Sigma11*Sigma22)^0.5;Rho

# A.2) Estimación de Rho
R=cor(b[,1],b[,2],method = "pearson" );R


### B) El Tau de Kendall
# B.1) El parámetro tau: 4 * E [ F(x1 ,x2) ] - 1
Tau=function(x){
	F=pmvnorm(lower=-Inf,upper=c(x[1],x[2]), mean = Mu, sigma = Cov)
	f=dmvnorm( c(x[1],x[2]) ,  mean = Mu , sigma = Cov)
	I=F*f;return(I)
}
EFx1x2=adaptIntegrate(Tau, lowerLimit = c(-Lx1, -Lx2), upperLimit = c(Lx1, Lx2))$integral
tau=4*EFx1x2-1;tau

# Si la distribución es normal
tau= 2 / pi * asin(Rho); tau

# B.2) Estimación del tau de Kendall
# B.2.1) Estimación no paramétrica
tauk=cor(b[,1],b[,2],method = "kendall" );tauk

# B.2.2) Estimador basado en una Cópula
# B.2.2.1) Con base en la relación entre "tau" y el parámetro de ligadura "theta"
tauCopula=(thetaE)/(thetaE+2);tauCopula

# B.2.2.2) Con base en la cópula estimada. Dará igual al anterior caso (B.2.2.1)
tauCopulaf=function(x){
	u=x[1];v=x[2]
	k = u^(-thetaE) + v^(-thetaE) -1
	if(k<0) k=0
	C=( k )^(-1/thetaE)
	c= ( k )^(-1/thetaE - 2) * (1+thetaE) * (u*v)^(-thetaE-1)
	I= C * c
	return(I)
}
tauCopula=4*adaptIntegrate(tauCopulaf, lowerLimit = c(0.0, 0.0), upperLimit = c(1, 1))$integral-1;tauCopula

# B.3) Uso del tau de Kendall para estimar el parámetro de ligadura de una Cópula
thetaE2 = 2*tauk/(1-tauk); thetaE2
thetaE;thetaE2
# Probabilidades acumuladas en los puntos (0,0) y (1.96 , 1.96):
# Basada en el estimador máximo verosímil
CopClay( c(0,0) , thetaE ) ; CopClay( c(1.96 , 1.96) ,  thetaE )
# Basada en el estimador obtenido del tau
CopClay( c(0,0) , thetaE2) ; CopClay( c(1.96 , 1.96) , thetaE2 )


### C) El Rho de Spearman
# C.1) El parámetro
# C.1.1) Una forma: 12 * E[ F(x1)*F(x2) ] - 3


# C.1.2) Otra forma. 12* II F(x1,x2) f(x1) f(x2)dx1 dx2 - 3


# En la normal: 
RhoS = asin(Rho/2) * 6 / pi;RhoS

# En general RhoS = cor (F(x) , F(y))



# C.2) Estimación del Rho de Spearman
# C.2.1) Estimación no paramétrica
Rs=cor(b[,1],b[,2],method = "spearman" );Rs

# C.2.2) Estimador basado en una Cópula
# C.2.2.1) Una forma: 12 * E[ U*V ] - 3 = 12 II u*v * c dudv - 3


# C.2.2.2) Otra forma. 12* II C(u,v) du dv - 3



### D) Sigma de Schweizer y Wolff
# D.1) El parámetro SigmaSW. II | F(x1,x2) - F(x1)*F(x2) | dx dy


# B.2) Estimación de la SigmaSW.  II | C(u,v) - u*v | du dv


#############################################################################
########## Simulación para analizar propiedades de los estimadores ##########
#############################################################################


m=10000

# Preparación vectores de estimación no paramétrica
RS=NULL;taukS=NULL;RsS=NULL

# Preparación vectores de estimación cópula
tauCopulaS=NULL;RsCopulaS=NULL;SigSWES=NULL
thetaES=NULL

for(i in 1:m){

	b<-rmvnorm(n, mean = Mu, sigma = Cov, method=c("chol"))#("eigen", "svd", "chol")

	taukS[i]=cor(b[,1],b[,2],method = "kendall")
	RsS[i]=cor(b[,1],b[,2],method = "spearman")
	RS[i]=cor(b[,1],b[,2],method = "pearson")
	
	# Estimación copula

	
	# Estimación tau kendall basado en copulas
	
	# Estimación rho Spearman basado en copulas
	
	# Sigma de Schweizer y Wolff

}

# Medias y desviaciones estándar de los vectores simulados
