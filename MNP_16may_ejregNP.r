#-------------------------------------
# Regresión y ANCOVA no paramétrico
#-------------------------------------

rm.list=ls()

#--------------------------
# Librerías
#--------------------------

library(sm)
library(car)

#--------------------------
# Lectura base de datos
#--------------------------

Ramiro=read.table("Ramirob.txt", head=TRUE, dec=",")
head(Ramiro)
attach(Ramiro)
summary(Ramiro)
Epoca=as.factor(Epoca)
is.factor(Epoca)

#--------------------------
# Graficos bivariados
#----------------------------

pairs(~Talla+Peso+FB, data=Ramiro, col=2, pch=16)
scatterplotMatrix(~Talla+Peso+FB, data=Ramiro)


#-----------------------------------
# Regresión No Paramétrica
#-----------------------------------


regnp=sm.regression(Talla,Peso, model="none", col=2, lwd=2,se=T)

regnp2=sm.regression(Talla,FB, model="none", col=2, lwd=2,se=T)

regnp3=sm.regression(Peso, FB, model="none", col=2, lwd=2,se=T)


# No efecto

regnp=sm.regression(Talla,Peso, model="no effect", col=2, lwd=2)
regnp2=sm.regression(Talla,FB, model="no effect", col=2, lwd=2)
regnp3=sm.regression(Peso, FB, model="no effect", col=2, lwd=2)

# Efecto lineal
regnp=sm.regression(Talla,Peso, model="linear", col=2, lwd=2)

# ANCOVA

ancova1 <-sm.ancova(Talla,Peso, Sitio, model="equal") #PRUEBA GLOBAL, como pvalor=0.7245m No R.H0
#suponiendo que rechazamos H0, hacemos las comparaciones por pares
Guajira=subset(Ramiro,Sitio=="G")
Isla=subset(Ramiro,Sitio=="I")
Flores=subset(Ramiro,Sitio=="F")

GI=rbind(Guajira, Isla)
#H0: las curvas de regresión son iguales (model="equal")
ancova2 <-sm.ancova(GI$Talla,GI$Peso, GI$Sitio, model="equal") #está comparando dos Sitios Guajira e Isla, y el p-valor=0.4829 ent. NO R.H0
GF=rbind(Guajira, Flores)
ancova3 <-sm.ancova(GF$Talla,GF$Peso, GF$Sitio, model="equal")
IF=rbind(Isla, Flores)
ancova4 <-sm.ancova(IF$Talla,IF$Peso, IF$Sitio, model="equal")

ancova5=sm.ancova(GI$Talla,GI$FB, GI$Sitio, model="equal")
ancova6=sm.ancova(GF$Talla,GF$FB, GF$Sitio, model="equal")
ancova7=sm.ancova(IF$Talla,IF$FB, IF$Sitio, model="equal")
