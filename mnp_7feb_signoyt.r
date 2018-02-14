#noparamétricos - analisis de potencia (prueba t y signo)
#Asignación de parámetros, tamaño de muestra, rango bajo H1 y alpha.

mu0=28
sigma=1
n=16
l=5*sigma/sqrt(n) # Longitud grande (5 desviaciones estandar)
alpha=0.05

# 1) Prueba-t: A cola derecha, asumiendo H1:mu > mu0 (Tipo B).
# Potencia = pr(media > mu0 + z.alpha / sqrt(n) | mu1)

mu1.L = seq(from= mu0-l, to= mu0, by= 0.005) ; mu1.L
mu1.U = seq(from= mu0, to= mu0+l, by= 0.005) ; mu1.U
mu1= c(mu1.L, mu1.U) ; mu1
z.alpha= qnorm(alpha, lower.tail = FALSE) ; z.alpha 
c1=mu0+z.alpha*sigma/sqrt(n) ; c1 # Limite para rechazar H0
z=(c1-mu1.U)/sigma*sqrt(n) ;z
potencia.tB=1-pnorm(z) ; potencia.tB
plot(mu1.U, potencia.tB) 


  
# Tarea H1:mu < mu0 y H1:mu != mu0.

# 2) Prueba-t: A dos colas.

z.alpha2=-qnorm(alpha/2)
c1=mu0-z.alpha2*sigma/sqrt(n)
c2=mu0+z.alpha2*sigma/sqrt(n)
z1=(c1-mu1)/sigma*sqrt(n)
# algunos valores de z1 son positivos, lo cual no tiene sentido porque c1-mu1 <0 dado que min(mu1)=26.75 y c1=26.04
z2=(c2-mu1)/sigma*sqrt(n)

potencia=pnorm(z1)+(1-pnorm(z2))
plot(mu1, potencia)
lines(mu1, rep(0.05, length(mu1))) 
# No sabemos que hace axis(...)
axis(2, at= seq (0.01, 0.05))

  
# 3) Prueba-signo
#prueba-signo para hipotesis a una cola
#Rechazar H0 si S=suma de signos + en la diferencia es mayor
  
z=mu0-mu1.U/sigma ;z
theta1=1-pnorm(z) ;theta1
potenciaSigno=1-pbinom(11, 16, theta1, lower.tail = TRUE) # Equivalente a p(s=>12|theta1)
plot(mu1.U, potenciaSigno)

#Comparacion con la prueba t
plot(mu1.U, potencia.tB)
lines(mu1.U, potenciaSigno)

# En este escenario la desviacion estandar es baja, mirar en otros casos. Mirar prueba del signo de Wilcoxon.
# Trabajo: Potencia bajo H1:mu < mu0 y H1:mu != mu0, prueba del signo de Wilcoxon.
