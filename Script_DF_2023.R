#####################################################################
####### SCRIPT para el TP ESTIMACION de DIVERSIDAD FUNCIONAL ########
################# ECOLOGIA DE COMUNIDADES - 2023 ####################
#####################################################################

#####PRIMERA PARTE########

rm(list=ls())
ls()

## EStablecemos el directorio de trabajo, donde estan los archivos
setwd("F:/Juli/Documents/Eco Comunidades")


## Carga de la informacion base
## Abri el archivo de Excel con todos los datos 
## y generá los archivos txt de Sitios y Rasgos, luego cargá los datos en R
Rasgos_Especies <- read.table("Rasgos.txt", header=T, row.names = 1)
Sitios <- read.table("Sitios.txt", header=T)

head(Rasgos_Especies) #para chequeo
head(Sitios) #para chequeo

# Aclaración nomenclatura de los sitios
# Sitio 1: Humedal
# Sitio 2: Agroecosistema 
# Sitio 3: Forestacion madura de Eucalipto


# Cargamos los paquetes necesarios... ¡y crucemos los dedos!
library(FD) #para cargar el paquete estadistico que calcula el indice FD de Petchey y Gaston 2002 (basado en el árbol funcional)
library(cluster)

Rasgos_Dis <-gowdis(Rasgos_Especies) # generá la matriz de distancia entre especies basada en los estados para los rasgos
# El método de Gower permite calcular distancias basandose en variables cuantitativas y cualitativas 

Rasgos_Arbol <- hclust(Rasgos_Dis, method='ward.D2') # Construimos el árbol para las especies de todos los sitios por el método de Ward

windows(10,10) 
plot(Rasgos_Arbol, hang = -1) 

## Si FD es la suma de las ramas que unen a las especies para un sitio:
## Analizá el árbol. 
## Suponé una comunidad 1 que presenta las siguientes especies: Tu_am, Tu_ru, My_mo, Ze_au, Co_pi, No_cy, Mo_ba y Zo_ca
##Rta: FD=

## ¿Como se modifica FD si por un cambio ambiental se pierden en los ensambles Tu_am, Co_pi y Mo_ba?
##Rta: No hay cambios significativos porque hay alta redundancia de esos caracteres

## ¿Y si solo se pierde No_cy?
##Rta:disminuye mucho

## Ahora calcule FD usando VEGAN con la funcion treedive
## ¿en cual ambiente espera que FD sea menor? Pensar en términos de filtros ambientales y redundancia funcional

?treedive() 
FD <- treedive(Sitios, Rasgos_Arbol, match.force = TRUE)  # Calculamos el índice de diversidad funcional FD (P&G 2002) a partir del árbol
FD
windows(12,4)
par(mfrow=c(1,3))
barplot(FD, ylim=c(0,max(FD)), ylab="FD")

## ¿Como se relaciona la DF con la riqueza taxonómica? 
##Rta.: 

?specnumber

Sobs  <- specnumber(Sitios)       # Riqueza observada
Sobs
barplot(Sobs, ylim=c(0,max(Sobs)), ylab="xxxxx")# poné nombre al eje y
plot(Sobs,FD, ylab="xxxx", xlab = "xxxx")# ponele nombre a los ejes

#Interpretar los resultados


###SEGUNDA PARTE#####

###EXTRAER INFO DE UN SOLO SITIO (Forestal, el 3)
#extraemos los nombres de las especies para el tercer sitio
spp.s3<- names(Sitios[3, Sitios[3,] > 0]) 
#llamar a spp.s3 para ver si lo hizo bien

#extraemos los rasgos para esas especies
rasgos.s3<-Rasgos_Especies[spp.s3,]    
#mirar qué hizo la sentencia anterior

mat.dist3<- vegdist (rasgos.s3, method="gower") #distancia entre especies de ese sitio
dendro.s3<-hclust(mat.dist3, method= "ward.D2") #construccion del arbol para el sitio elegido
windows(4,4) 
plot(dendro.s3, hang = -1)


FD3<- treedive(Sitios[3, Sitios[3,] > 0], dendro.s3, match.force = F)
FD3


##¿Por qué el FD que obtuvimos antes para el sitio forestal no es igual a FD3?
FD
FD3

####TERCERA PARTE##### 

### FD para un solo rasgo (Insectivoro)para todos los sitios.
## ¿Se animan? 
# Una ayudita...


#1. extraer el rasgo de interes para esas especies
#
INSEC<-Rasgos_Especies[XXXXXXX,xxxx]


#3. Calcular la distancia de Gower entre especies en base a si son insectivoras o no
#4. Calcular el dendrograma
#5. Graficarlo

####CUARTA PARTE######### 

####NOS INDEPENDIZAMOS DEL NUMERO DE ESPECIES################

# necesitamos cargar un paquete que calcula otros índices de diversidad funcional
?FD
install.packages("FD")

library(FD)

?fdisp ##esta función va a calcular el indices de diversidad funcional de Laliberté y Legendre (2010). 
      ##Es clave leer el argumento

## Vamos a pedir que calcule FDis que sirve para este tipo de variables y no está inflado por la cantidad de especies
Sitios.Mx<- as.matrix(Sitios)
FDis <- fdisp(Rasgos_Dis, Sitios.Mx)
FDis$FDis

windows(8,4)
par(mfrow=c(1,2))
barplot(FD, ylim=c(0,max(FD)), ylab="FD")
barplot(FDis$FDis, ylim=c(0,max(FDis$FDis)), ylab = "FDis")


#¿Conclusiones?


#####QUINTA PARTE#### 

################################################################################################ Beta funcional

#Volvamos a los índices basado en árboles

install.packages("picante")
library(picante) #esta libreria trabaja con objetos (arboles) tipo phylo (una extension distinta a hclust).
# Los arboles obtenidos con hclus deben transformarse

Arbol.phylo <- as.phylo(Rasgos_Arbol) 

?phylosor #chusmear que hace esta funcion

FSor  <- phylosor(Sitios, Arbol.phylo) # tambien se usa en arboles filogeneticos
FSor
DistFunc<- 1-FSor # ¿por que hacemos esto?
DistFunc

Hclust.FSor <- hclust(DistFunc, method='single') # ¿que elementos agrupamos aqui?
windows()
plot(Hclust.FSor) # visualicemos...

# ¿cuales son los ambientes mas parecidos funcionalmente? ¿Y los mas diferentes?

####SEXTA PARTE#####

#####ACOPLAMIENTO DE CAMBIOS EN LAS DIVERSIDADES TAXONOMICA Y FUNCIONAL
### Beta funcional vs taxonómica

DistTax <- vegdist(Sitios, method="bray", binary=TRUE)#Sorensen taxonomico
DistTax
windows()
plot(DistFunc, DistTax, xlim=c(0,1), ylim=c(0,1))
abline(0,1)

# ¿Hay tendencia a que el cambio funcional este acoplado al cambio taxonomico?


#####SEPTIMA PARTE- TAREA EN CASA########

############################################################################
### otros alfa funcionales con el paquete FD
#####################################################################################
#sacar un # al inicio de las sentencias que siguen para correrlo en sus casas

library(FD) #si ya lo tenían cargado esta sentencia no hay que correrla

#?FD
#?dbFD 
##esta función va a calcular varios indices de diversidad funcional. 
##Es clave leer el argumento
?dbFD
DF_Varios <- dbFD(Rasgos_Dis, Sitios, calc.FGR =T)

## Si en el argunmento agregamos calc.FGR= T le estamos pidiendo que forme grupos a posteriori
## al ejecutarse y llegar a los grupos funcionales va a preguntar como los queremos definir,
## por altura (h) o por cantidad (g). Prueben. Elegir una altura de corte (¡clase de clusters!)... 
## o cantidad de grupos (¿sabemos cuantos habrá?)

#Tira un warning, no se preocupen que la función corre igual.

DF_Varios
## leer el argumento de dbFD para entender lo que proporciona la salida

windows(9,3)
par(mfrow=c(1,3))
barplot(DF_Varios$FRic, type="h", ylab = "FRic")
barplot(DF_Varios$FEve, type="h", ylab = "FEve")
barplot(DF_Varios$FDiv, type="h", ylab = "FDiv")


###Comparemos FDis con RaoQ

windows(6,3)
par(mfrow=c(1,2))
barplot(DF_Varios$FDis, ylab = "FDis")
barplot(DF_Varios$RaoQ, ylab = "RaoQ")
