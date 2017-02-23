# directorio de trabajo
setwd("~/2_Otros_Proyectos_Conabio_2/Mario_Carlos_Errores")

# paquetes
library(randomGLM)
library(randomForest)
library(plyr)
library(raster)

# cargamos datos
datos<- read.csv(file="dat.csv",head=TRUE,sep=",")
datos <- as.data.frame(dat)

# muestreo para validacion cruzada manual hecha por julian
dobleces <- function(datos,dobleces=10,seed=NULL){
  if (!is.null(seed)) set.seed(seed)
  
  # inicializamos los dobleces
  conjuntos <- list()
  # tamaño de las muestras
  tam <- floor(nrow(datos)/dobleces)

  # muestreos
  index <- 1:nrow(datos)
  for (i in 1:(dobleces-1)){
    seleccion <- sample(index,tam)
    conjuntos[[i]] <- datos[seleccion, ]
    index <- index[! index %in% seleccion]
  }
  conjuntos[[dobleces]]<-datos[index,]
  return(conjuntos)
}

# cabecera
head(datos)
nombres <- names(datos)

# separamos la base en 10 conjuntos usando la funcion anterior
# todo los próximos resultados por lo tanto están basados en
# 10-fold cross validation
conjuntos <- dobleces(datos,10)

# inicializamos resultados finales
resultadosfinales <- matrix(0,nrow(datos),3)
colnames(resultadosfinales)<-c("predicciones","limsup","liminf")

for (f in 1:10){

# train - test a lo largo de los 10-folds
train <- ldply(conjuntos[-f], data.frame)
names(train) <- nombres
test <- ldply(conjuntos[f], data.frame)
names(test) <- nombres

# numero de modelos en los que estará basado el
# randomGLM
nummodels <- 100

# modelo
RGLM = randomGLM(train[,13:18], train[,11], classify=FALSE, keepModels=TRUE, nThreads=1,nBags=nummodels,
                 maxInteractionOrder=1)
# models in each bag
modelbag <- RGLM$models

# inicializar la matriz que guardara las prediccions individuales
preds <- matrix(0,nrow(test),nummodels)
limsup <- matrix(0,nrow(test),nummodels)
liminf <- matrix(0,nrow(test),nummodels)

# pruebas intervalos de confianza de la prediccion
#prueba1 <- predict(modelbag[[1]],test,se.fit=TRUE)
#prueba2 <- predict(modelbag[[1]], test, interval="predict")

# predicciones individuales
for (i in 1:100){
  prediccion <- predict(modelbag[[i]],test,interval="predict")
  preds[,i]<- prediccion[,1]
  limsup[,i]<- prediccion[,3]
  liminf[,i]<- prediccion[,2]
  remove(prediccion)
}

# predicciones usuales
resultadosfinales[(1+(f-1)*nrow(test)):(f*nrow(test)),1] <- rowMeans(preds)

# promedio limites superiores
resultadosfinales[(1+(f-1)*nrow(test)):(f*nrow(test)),2] <- max(limsup)

# promedio limites inferiores
resultadosfinales[(1+(f-1)*nrow(test)):(f*nrow(test)),3] <- min(liminf)

}

head(resultadosfinales)

# datos originales en el orden de la validacion cruzada
datos <- ldply(conjuntos, data.frame)
names(datos) <- nombres

# correlacion observado vs predicho 
cor(datos[,11],resultadosfinales[,1])

# errores absolutos
errores <- abs(datos[,11]-resultadosfinales[,1])

# correlacion entre el tamaño de los intervalos y los errores
cor(errores,resultadosfinales[,2]-resultadosfinales[,3]) #muy baja

# la relacion entre los dos es muy baja pero positiva, asi que de cierta manera
# puedes asumir que son intervalos de predicción que incluso podrias graficar:

# base completa de la zona para hacer el mapa una vez que se tenga
# el modelo
covar<-read.table ("COV1.txt", sep="\t", header=T)

# modelo con datos completos
RGLM = randomGLM(datos[,13:18], datos[,11], classify=FALSE, keepModels=TRUE, nThreads=1,nBags=nummodels,
                 maxInteractionOrder=1)
# models in each bag
modelbag <- RGLM$models

# inicializar la matriz que guardara las prediccions individuales
preds <- matrix(0,nrow(covar),nummodels)
limsup <- matrix(0,nrow(covar),nummodels)
liminf <- matrix(0,nrow(covar),nummodels)


# predicciones individuales
for (i in 1:100){
  prediccion <- predict(modelbag[[i]],covar,interval="predict")
  preds[,i]<- prediccion[,1]
  limsup[,i]<- prediccion[,3]
  liminf[,i]<- prediccion[,2]
  remove(prediccion)
  
}

# predicciones usuales
predicciones <- rowMeans(preds)
limitesuperior <- apply(limsup, 1, max)
limiteinferior <- apply(liminf, 1, min)

# juntamos la info
datos <- data.frame(predicciones,desviacionesestandar=limitesuperior-limiteinferior,limitesuperior,
                    limiteinferior,X=covar$X,Y=covar$Y)
# clean memory
remove(covar)
remove(preds)
remove(limsup)
remove(liminf)

# varianza de las predicciones
scatter.fill(x=datos$Y,y=datos$X,z=datos$desviacionesestandar,nlevels=6,
             pch=".",main="incertitumbre de las pred")
# limite superior
scatter.fill(x=datos$Y,y=datos$X,z=datos$limitesuperior,nlevels=6,
             pch=".",main="limite superior")
# limite inferior
scatter.fill(x=datos$Y,y=datos$X,z=datos$limiteinferior,nlevels=6,
             pch=".",main="limite inferior")
