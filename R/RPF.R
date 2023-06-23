#' Regresión Polinomial
#'
#' Realiza una regresión polinomial a partir del grado indicado
#'
#' @param x (vector) Datos de la variable independiente (explicativa)
#' @param y (vector) Datos de la variable dependiente (respuesta)
#' @param grado (valor) Grado al que se desea realizar el polinomio
#' @return Una lista con los valores de Y ajustada, los coeficientes de regresión,el error estandar, la R^2 y la varianza
#' @export
#'
#' @examples
#'\dontrun
#' #Cargando la libreria
#'library(ReP)
#'
#'Datos para el modelo
#'X<-(1, 3, 5, 7, 9)
#'Y<-(4, 8, 12, 16, 20)
#'salida <-  regpol(x=X,y=Y, grado=3)
#'
regpol<- function(x,y, grado)
{
  x<-x
  y<-y
  grado<-grado
  #Crear matriz para las potencias de x
  #Crear una matriz vacia
  matrizX<- matrix(NA, nrow = grado + 1, ncol = grado + 1)
  # Llenar la matriz
  for (i in 0:grado) {
    for (j in 0:grado) {
      matrizX[i + 1, j + 1] <- sum(x^(i+j))
    }
  }
  #----------------------------------------------------------------------
  #Crear una matriz para los valores de y
  # Crear una matriz vacía para almacenar los valores
  matrizY <- matrix(NA, nrow = grado + 1, ncol = 1)
  # Calcular las sumatorias
  matrizY[1, 1] <- sum(y)
  for (i in 1:grado) {
    matrizY[i + 1, 1] <- sum(x^i * y)
  }
  #-----------------------------------------------------------------------
  #Obtener los coeficientes para el polinomio
  matrizInv<- solve(matrizX)#Matriz inversa
  #calcular coeficientes
  coefB<- matrizInv%*%matrizY
  #-----------------------------------------------------------------------
  #Calcular el error
  n <- length(x)
  matrizD <- matrix(NA, nrow = length(x), ncol = grado + 1)
  for (i in 0:grado) {
    matrizD[,i + 1] <- x^i
  }
  PredY <- matrizD %*% coefB
  sr<- sum((y-PredY)^2)
  error <- sqrt(sr/(n-(grado+1)))
  #-----------------------------------------------------------------------------
  #Rcuadrada
  tss<- sum((y-mean(y))^2)
  sigma<-sum((y-PredY)^2)/(n-grado+1)
  r2<- 1-(sqrt(sigma)/sqrt(tss))
  #varianza
  varianza<- sum((x-mean(x))^2)/n
  #------------------------------------------------------------------------------
  #r, y observado, ajustado, error, varianza
  #error estandar de los coeficientes
  invxx <- solve(t(matrizD) %*% matrizD)#Matriz inversa x´x
  matrizCov <- invxx * error^2 #Matriz de covarianzas
  erroresBetas <- sqrt(diag(matrizCov)) #Diagonal de la matriz
  #--------------------------------------------------------------------------
  #Significancia de las betas
  valorT<- coefB/erroresBetas
  gradosL<-n-(grado+1)
  valorP<- 2 * (1 - pt(abs(valorT), gradosL))
  #Crear lista de regreso
  df <- data.frame(Betas = coefB,
                   Error_estand = erroresBetas,
                   Vt = valorT,
                   vp = valorP)
  listaFinal<- list(df,"Error estandar residual"=error, "r2"=r2,
                    "Varianza"=varianza)
  return(listaFinal)
}
