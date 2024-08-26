#' Computes the directed dependence coefficient.
#'
#' @description
#' The directed dependence coefficient (\code{didec}) estimates the degree of directed dependence of a random vector Y on a random vector X, based on an i.i.d. sample of (X,Y).
#'
#' @param X A numeric matrix or data.frame/data.table. Contains the predictor vector X.
#' @param Y A numeric matrix or data.frame/data.table. Contains the response vector Y.
#' @param perm A logical. If \code{True} a version of \code{didec} is computed that takes into account the permutations (specified by \code{perm.method}) of the response variables.
#' @param perm.method An optional character string specifying a method for permuting the response variables. This must be one of the strings \code{"sample"}, \code{"increasing"}, \code{"decreasing"} (default) or \code{"full"}. The version \code{"full"} is invariant with respect to permutations of the response variables.
#'
#' @details
#' The directed dependence coefficient (didec) is an extension of Azadkia & Chatterjee's measure of directed dependence (Azadkia & Chatterjee, 2021) to a vector of response variables introduced in (Ansari & Fuchs, 2023).
#' Its calculation is based on the function \code{codec} which estimates Azadkia & Chatterjee’s measure of directed dependence and is provided in the R package \code{FOCI}.
#'
#' By definition, \code{didec} is invariant with respect to permutations of the variables within the predictor vector X. Invariance with respect to permutations within the response vector Y is achieved by computing the arithmetic mean over all possible (or chosen) permutations.
#' In addition to the option \code{"full"} of running all \eqn{q!} permutations of \eqn{(1, ..., q)}, less computationally intensive options are also available (here, \eqn{q} denotes the number of response variables): a random selection of \eqn{q} permutations \code{"sample"}, cyclic permutations such as \eqn{(1,2,...,q)}, \eqn{(2,...,q,1)} either \code{"increasing"} or \code{"decreasing"}.
#' Note that when the number of variables \eqn{q} is large, choosing \code{"full"} may result in long computation times.
#'
#' @return The degree of directed dependence of the random vector Y on the random vector X.
#'
#' @export
#'
#' @author Yuping Wang, Sebastian Fuchs, Jonathan Ansari
#'
#' @references
#' M. Azadkia, S. Chatterjee, A simple measure of conditional dependence, Ann. Stat. 49 (6), 2021.
#'
#' J. Ansari, S. Fuchs, A simple extension of Azadkia & Chatterjee's rank correlation to multi-response vectors, Available at https://arxiv.org/abs/2212.01621, 2024.
didec <- function(X, Y,
                  perm = FALSE, perm.method = c("decreasing")) {
  Idx.Perm <- c("sample", "increasing", "decreasing", "full")
  if (!perm.method %in% Idx.Perm) {
    stop("'perm.method' should be one of 'sample', 'increasing', 'decreasing', 'full'.")
  }

  X <- as.data.frame(X)
  Y <- as.data.frame(Y)
  if (nrow(X) != nrow(Y)) {
    stop("The number of rows of X and Y should be equal. ")
  }


  if (!perm) {
    DIDEC <- Codec.Tq(X, Y)
  }
  if (perm) {
    DIDEC <- Codec.Tq.Perm(X, Y, method = perm.method)
  }

  return(DIDEC)
}
