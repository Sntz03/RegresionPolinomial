#' Predicciones de la regresión polinomial
#'
#' Realiza predicciones a partir de un modelo dado
#'
#' @param x (vector) Datos de la variable independiente (explicativa)
#' @param betas (vector) Coeficientes del modelo polinomial
#' @return Una matriz con los valores de Y calculados
#' @export
#'
#' @examples
#'\dontrun
#'#Directorio de trabajo
#'ruta<- "~/ruta/datos.csv"
#'(datos<- read.csv(ruta))
#'
#'#Cargando libreria
#'library(ReP)
#'
#'b<- regpol(x=datos$x, y=datos$y, 4)
#'beta<-b[[1]][1]
#'nuevaX<-c(5,6,7,8,9,10)
#'
#'prediccion(nuevaX, beta)


prediccion <- function(x, betas) {
  n <- length(x)
  matrizB <- as.matrix(betas)
  grado <- length(matrizB) - 1
  y <- matrix(NA, nrow = n, ncol = 1)
  for (i in 1:n) {
    y[i] <- 0
    for (j in 0:grado) {
      y[i] <- y[i] + (matrizB[j+1, 1] * x[i]^j)
    }
  }
  return(y)
}
