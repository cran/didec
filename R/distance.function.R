#' distance function based on T^q
#'
#' @param X a data frame for vector X
#' @param Y a data frame for vector Y
#' @param perm permuted version or not
#' @param perm.method permutation methods: sample / increasing / decreasing / full
#' @param mutual use mutual perfect dependence or not
#'
#' @return a value for distance between two vectors
#'
#' @keywords internal
dist.Tq <- function(X, Y,
                    perm = TRUE, perm.method = c("decreasing"),
                    mutual = FALSE) {
  Idx.Perm <- c("sample", "increasing", "decreasing", "full")
  if (!perm.method %in% Idx.Perm) {
    stop("'perm.method' should be one of 'sample', 'increasing', 'decreasing', 'full'.")
  }

  X <- as.data.frame(X)
  Y <- as.data.frame(Y)
  if (nrow(X) != nrow(Y)) {
    stop("The numbers of rows of X and Y should be equal.")
  }

  if (!mutual) {
    D <- (1 - didec(X, Y, perm = perm, perm.method = perm.method)) * (1 - didec(Y, X, perm = perm, perm.method = perm.method))
  }

  if (mutual) {
    D <- 1 - mean(didec(X, Y, perm = perm, perm.method = perm.method), didec(Y, X, perm = perm, perm.method = perm.method))
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
  Idx <- c("kendall", "footrule")
  if (!method %in% Idx) {
    stop("'method' should be one of 'kendall','footrule'.")
  }
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)
  if (nrow(X) != nrow(Y)) {
    stop("Objects of different size")
  }
  df <- data.frame(X, Y)
  D <- 1 - concor.M(df, method = method)
  return(D)
}

#' Distance Matrix Computation using distance function based on T^q
#'
#' @param X a data frame for vector X
#' @param mutual use type B function (mutual perfect dependence) or not
#'
#' @return an object of class "dist"
#'
#' @importFrom stats as.dist
#'
#' @keywords internal
dist.mat.T <- function(X, mutual = FALSE) {
  df <- as.data.frame(X)
  dX <- length(df)
  cn <- colnames(df)
  dist <- matrix(0, dX, dX)
  colnames(dist) <- rownames(dist) <- cn
  for (i in c(1L:(dX - 1))) {
    for (j in c((i + 1):dX)) {
      dist[i, j] <- dist[j, i] <- dist.Tq(df[, i], df[, j], perm = FALSE, mutual = mutual)
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
#' @importFrom factoextra get_dist
#' @importFrom copBasic footCOP
#' @importFrom stats as.dist
#'
#' @keywords internal
dist.mat.concor <- function(X, method = c("footrule")) {
  Idx.Method <- c("kendall", "footrule")
  if (!method %in% Idx.Method) {
    stop("'method' should be one of 'kendall','footrule'.")
  }
  df <- as.data.frame(X)
  if (method == "kendall") {
    dist <- get_dist(t(X), method = "kendall")
  }
  if (method == "footrule") {
    dX <- length(df)
    cn <- colnames(X)
    dist <- matrix(0, dX, dX)
    colnames(dist) <- rownames(dist) <- cn
    for (i in c(1:(dX - 1))) {
      for (j in c((i + 1):dX)) {
        dist[i, j] <- dist[j, i] <- 1 - footCOP(para = df[, c(i, j)], as.sample = TRUE)
      }
    }
    dist <- as.dist(dist)
  }
  return(dist)
}
