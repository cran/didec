#' distance function based on T
#'
#' @param X a data frame for vector X
#' @param Y a data frame for vector Y
#' @param estim.method An optional character string specifying a method for estimating the directed dependence coefficient.
#' @param mutual use mutual perfect dependence or not
#'
#' @return a value for distance between two vectors
#'
#' @keywords internal
dist.Tq <- function(X, Y, 
                    estim.method = c("copula"),
                    mutual = FALSE) {
  
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)
  
  if (nrow(X) != nrow(Y)) {
    stop("The numbers of rows of X and Y should be equal.")
  }
  
  if (!mutual) {
    D <- (1 - didec(X, Y, estim.method = estim.method)) * (1 - didec(Y, X, estim.method = estim.method))
  }
  
  if (mutual) {
    D <- 1 - mean(didec(X, Y, estim.method = estim.method), didec(Y, X, estim.method = estim.method))
  }
  
  if (D < 0) {
    D <- 0
  }
  if (D > 1) {
    D <- 1
  }
  return(D)
}

#' distance function based on multivariate concordance measures
#'
#' @param X a data frame for vector X
#' @param Y a data frame for vector Y
#' @param method kendall / footrule
#' 
#' @return a value for distance between two vectors
#'
#' @keywords internal
dist.concor.M <- function(X, Y, method = c("footrule")) {
  
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)
  
  if (nrow(X) != nrow(Y)) {
    stop("The numbers of rows of X and Y should be equal.")
  }
  
  df <- data.frame(X, Y)
  D <- 1 - concor.M(df, method = method)
  return(D)
}

#' Distance Matrix Computation using distance function based on T^q
#'
#' @param X a data frame for vector X
#' @param mutual use type B function (mutual perfect dependence) or not
#' @param estim.method An optional character string specifying a method for estimating the directed dependence coefficient.
#' @return an object of class "dist"
#'
#' @importFrom stats as.dist
#'
#' @keywords internal
dist.mat.T <- function(X,
                       estim.method = c("copula"), 
                       mutual = FALSE) {
  df <- as.data.frame(X)
  dX <- length(df)
  cn <- colnames(df)
  dist <- matrix(0, dX, dX)
  colnames(dist) <- rownames(dist) <- cn
  for (i in c(1L:(dX - 1))) {
    for (j in c((i + 1):dX)) {
      dist[i, j] <- dist[j, i] <- dist.Tq(df[, i], df[, j], estim.method = estim.method, mutual = mutual)
    }
  }
  dist <- as.dist(dist)
  return(dist)
}

#' Distance Matrix Computation using distance function based on multivariate concordance measures
#'
#' @param X a data frame for vector X
#' @param method kendall / footrule
#'
#' @return an object of class "dist"
#'
#' @importFrom stats as.dist
#'
#' @keywords internal
dist.mat.concor <- function(X, method = c("footrule")) {

  df <- as.data.frame(X)
  dX <- length(df)
  cn <- colnames(X)
  dist <- matrix(0, dX, dX)
  colnames(dist) <- rownames(dist) <- cn
  for (i in c(1:(dX - 1))) {
    for (j in c((i + 1):dX)) {
      dist[i, j] <- dist[j, i] <- dist.concor.M(df[, i], df[, j], method = method)
    }
  }
  dist <- as.dist(dist)
  return(dist)
}
