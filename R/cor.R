#' Estimate for T^q(Y|X) based on function Codec
#'
#' @param X a data frame for input vector X
#' @param Y a data frame for output vector Y
#'
#' @return the value of T^q(Y|X)
#'
#' @importFrom FOCI codec
#'
#' @keywords internal
Codec.Tq <- function(X, Y) {
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)
  dY <- length(Y)
  if (nrow(X) != nrow(Y)) {
    stop("Objects of different size.")
  }
  ZW <- numeric()
  weight <- numeric()
  weight[1] <- 0
  ZW[1] <- as.numeric(codec(Y[, 1], X))
  if (dY > 1) {
    for (i in 2L:dY) {
      weight[i] <- as.numeric(codec(Y[, i], Y[, 1L:(i - 1)]))
      ZW[i] <- as.numeric(codec(Y[, i], data.frame(X, Y[, 1L:(i - 1)])))
    }
  }
  return(CodecTq = 1 - (dY - sum(ZW)) / (dY - sum(weight)))
}


#' Estimate for T^q_bar(Y|X) based on function Codec & a sample of all / all increasing / all decreasing permutations
#'
#' @param X a data frame for input vector X
#' @param Y a data frame for output vector Y
#' @param method permuatation methods: sample / increasing / decreasing / full
#'
#' @return the value of T^q_bar(Y|X)
#'
#' @importFrom gtools permutations
#'
#' @keywords internal
Codec.Tq.Perm <- function(X, Y, method = c("sample")) {
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)
  if (nrow(X) != nrow(Y)) {
    stop("Objects of different size")
  }
  dY <- length(Y)
  if (dY > 1) {
    if (method == "sample") {
      perm <- permutations(n = dY, r = dY, v = 1L:dY)
      perm <- perm[sample(1L:factorial(dY), size = dY, replace = FALSE), ]
    }
    if (method == "increasing") {
      perm <- matrix(1L:dY, dY, dY + 1, byrow = T)[, 1L:dY]
    }
    if (method == "decreasing") {
      perm <- matrix(dY:1L, dY, dY + 1, byrow = T)[, 1L:dY]
    }
    if (method == "full") {
      perm <- permutations(n = dY, r = dY, v = 1L:dY)
    }
    Results <- numeric()
    for (l in seq_len(nrow(perm))) {
      Y <- Y[, perm[l, ]]
      Results[l] <- Codec.Tq(X, Y)
    }
    Result <- mean(Results)
  }
  if (dY == 1) {
    Result <- Codec.Tq(X, Y)
  }
  return(Result)
}


#' Estimation for multivariate concordance measures
#'
#' @param X a data frame for vector X
#' @param method kendall / footrule
#'
#' @return a value of the estimator for the multivariate concordance measures
#'
#' @keywords internal
concor.M <- function(X, method = c("footrule")) {
  Idx.Method <- c("kendall", "footrule")
  if (!method %in% Idx.Method) {
    stop("'method' should be one of 'kendall','footrule'.")
  }
  X <- as.data.frame(X)
  n <- nrow(X)
  d <- length(X)

  if (method == "kendall") {
    S <- c()
    xn <- c(1L:n)
    for (i in xn) {
      for (j in xn[xn != i]) {
        L <- c()
        for (k in 1L:d) {
          if (X[i, k] <= X[j, k]) {
            L <- c(L, TRUE)
          }
        }
        if (sum(L) == d) {
          S <- c(S, 1)
        }
      }
    }
    t <- ((2^d * sum(S)) / (n * (n - 1)) - 1) / (2^(d - 1) - 1)
    return(t)
  }

  if (method == "footrule") {
    R <- X
    for (j in 1L:d) {
      R[, j] <- rank(X[, j], ties.method = "random")
    }
    L <- c()
    for (i in 1L:n) {
      Li <- max(R[i, ]) - min(R[i, ])
      L <- c(L, Li)
    }
    phi <- 1 - ((d + 1) / (d - 1)) * (sum(L) / (n^2 - 1))
    return(phi)
  }
}
