#' Computes the directed dependence coefficient.
#'
#' @description
#' The directed dependence coefficient (\code{didec}) estimates the degree of functional dependence of a random vector Y on a random vector X, based on an i.i.d. sample of (X,Y).
#'
#' @param X A numeric matrix or data.frame/data.table. Contains the predictor vector X.
#' @param Y A numeric matrix or data.frame/data.table. Contains the response vector Y.
#' @param trans A logical. If \code{TRUE} the inputs of X are standardized (transformed) before \code{didec} is computed.
#' @param trans.method An optional character string specifying the data standardization method. This must be one of the strings \code{"standardization"} (default), \code{"rank"} or \code{"rescaling"}. \code{"standardization"} centers and scales each predictor to zero mean and unit variance (z-score). \code{"rank"} uses the rank of values instead of the values themselves. \code{"rescaling"} rescales each predictor to \eqn{[0,1]} (min–max normalization).
#' @param estim.method An optional character string specifying a method for estimating the directed dependence coefficient. This must be one of the strings \code{"codec"} or \code{"copula"} (default).
#' @param perm A logical. If \code{TRUE} a version of \code{didec} is computed that takes into account the permutations (specified by \code{perm.method}) of the response variables.
#' @param perm.method An optional character string specifying a method for permuting the response variables. This must be one of the strings \code{"sample"}, \code{"increasing"}, \code{"decreasing"} (default) or \code{"full"}. The version \code{"full"} is invariant under permutations of the response variables.
#'
#' @details
#' The directed dependence coefficient (didec) is an extension of Azadkia & Chatterjee's measure of functional dependence (Azadkia & Chatterjee, 2021) to a vector of response variables introduced in (Ansari & Fuchs, 2025).
#' \code{estim.method} specifies two methods for estimating the directed dependence coefficient. \code{"codec"} uses the function \code{codec} which estimates Azadkia & Chatterjee’s measure of functional dependence and is provided in the R package \code{FOCI}. \code{"copula"} estimates the directed dependence coefficient based on a dimension reduction principle; see (Fuchs 2024). The value returned by \code{didec} may be positive or negative. In the asymptotic limit, however, it is guaranteed to lie between \eqn{0} and \eqn{1}.
#'
#' By definition, \code{didec} is invariant under permutations of the variables within the predictor vector X. Invariance under permutations within the \eqn{q}-dimensional response vector Y is achieved by computing the arithmetic mean over all possible permutations.
#' In addition to the option \code{"full"} of running all \eqn{q!} permutations of \eqn{(1, ..., q)}, less computationally intensive options are also available: a random selection of \eqn{q} permutations \code{"sample"}, cyclic permutations such as \eqn{(1,2,...,q)}, \eqn{(2,...,q,1)} either \code{"increasing"} or \code{"decreasing"}.
#' Note that when the number of variables \eqn{q} is large, choosing \code{"full"} may result in long computation times.
#'
#' @return The degree of functional dependence of the random vector Y on the random vector X.
#'
#' @export
#'
#' @author Yuping Wang, Sebastian Fuchs, Jonathan Ansari
#'
#' @references
#' J. Ansari, S. Fuchs, A direct extension of Azadkia & Chatterjee's rank correlation to multi-response vectors, Available at \url{https://arxiv.org/abs/2212.01621}, 2025.
#' 
#' M. Azadkia, S. Chatterjee, A simple measure of conditional dependence, Ann. Stat. 49 (6), 2021.
#' 
#' S. Fuchs, Quantifying directed dependence via dimension reduction, J. Multivariate Anal. 201, Article ID 105266, 2024.
didec <- function(X, Y,
                  trans = FALSE, trans.method = c("standardization"),
                  estim.method = c("copula"),
                  perm = FALSE, perm.method = c("decreasing")) {
  
  X <- as.data.frame(X); Y <- as.data.frame(Y)
  p <- ncol(X); q <- ncol(Y)
  
  if (nrow(X) != nrow(Y)) {
    stop("The number of rows of X and Y should be equal.")
  }
  
  # remove NAs
  if (!(sum(is.na(X))==0 & sum(is.na(Y))==0)) {
    ok <- complete.cases(X,Y)
    if (sum(ok) < 2) {
      stop("Number of rows with no NAs should be bigger than 1.")
    }
    X <- X[ok,];Y <- Y[ok,]
  }
  
  # Y cannot be constant
  for (i in 1L:q){
    if (!(length(unique(Y[,i])) > 1)) {
      stop(paste0("Column ",i," of Y is constant."))
    }
  }
  
  for (i in 1L:p){
    if (!(length(unique(X[,i])) > 1)) {
      warning(paste0("Column ",i," of X is constant."))
    }
  }
  
  Idx.Estim <- c("copula", "codec")
  if (!estim.method %in% Idx.Estim) {
    stop("'estim.method' should be one of 'copula', 'codec'.")
  }
  
  if (trans) {
    
    # X cannot be constant
    for (i in 1L:p){
      if (!(length(unique(X[,i])) > 1)) {
        stop(paste0("Column ",i," of X is constant."))
      }
    }
    
    Idx.Trans <- c("standardization", "rank", "rescaling")
    if (!trans.method %in% Idx.Trans) {
      stop("'trans.method' should be one of 'standardization', 'rank', 'rescaling'.")
    }
    
    if (trans.method == "standardization") {
      for(i in 1L:p){
        X[,i] <- (X[,i] - mean(X[,i])) / sd(X[,i])
      }
    }
    if (trans.method == "rescaling") {
      for(i in 1L:p){
        X[,i] <- (X[,i] - min(X[,i])) / ( max(X[,i]) - min(X[,i]))
      }
    }
    if (trans.method == "rank") {
      for(i in 1L:p){
        X[,i] <- rank(X[,i], ties.method = "random")
      }
    }
  }

  if (estim.method == "copula") {

    if (!perm) {
      DIDEC <- Copula.Tq(X, Y)
    }
    if (perm) {
      
      Idx.Perm <- c("sample", "increasing", "decreasing", "full")
      if (!perm.method %in% Idx.Perm) {
        stop("'perm.method' should be one of 'sample', 'increasing', 'decreasing', 'full'.")
      }
  
      DIDEC <- Copula.Tq.Perm(X, Y, method = perm.method)
    }
  }

  if (estim.method == "codec") {

    if (!perm) {
      DIDEC <- Codec.Tq(X, Y)
    }
    if (perm) {
      
      Idx.Perm <- c("sample", "increasing", "decreasing", "full")
      if (!perm.method %in% Idx.Perm) {
        stop("'perm.method' should be one of 'sample', 'increasing', 'decreasing', 'full'.")
      }
      
      DIDEC <- Codec.Tq.Perm(X, Y, method = perm.method)
    }
  }

  return(DIDEC)
}
