#' Multivariate feature ordering by conditional independence.
#'
#' @description
#' A variable selection algorithm based on the directed dependence coefficient (\code{\link{didec}}).
#'
#' @param X A numeric matrix or data.frame/data.table. Contains the predictor vector X.
#' @param Y A numeric matrix or data.frame/data.table. Contains the response vector Y.
#' @param pre.selected An integer vector for indexing pre-selected predictor variables from X.
#' @param perm A logical. If \code{True} a version of \code{didec} is computed that takes into account the permutations (specified by \code{perm.method}) of the response variables.
#' @param perm.method An optional character string specifying a method in \code{didec} for permuting the response variables. This must be one of the strings \code{"sample"}, \code{"increasing"}, \code{"decreasing"} (default) or \code{"full"}. The version \code{"full"} is invariant with respect to permutations of the response variables.
#' @param autostop A logical. If \code{True} the algorithm stops at the first non-increasing value of \code{didec}.
#'
#' @details
#' \code{mfoci} is a forward feature selection algorithm for multiple-outcome data that employs the directed dependence coefficient (\code{didec}) at each step.
#' \code{mfoci} is proved to be consistent in the sense that the subset of predictor variables selected via \code{mfoci} is sufficient with high probability.
#'
#' If \code{autostop == TRUE} the algorithm stops at the first non-increasing value of \code{didec}, thereby selecting a subset of variables.
#' Otherwise, all predictor variables are ordered according to their predictive strength measured by \code{didec}.
#'
#' @return A data.frame listing the selected variables.
#'
#' @export
#'
#' @author Sebastian Fuchs, Jonathan Ansari, Yuping Wang
#'
#' @importFrom rlang is_empty
#'
#' @references
#' J. Ansari, S. Fuchs, A simple extension of Azadkia & Chatterjee's rank correlation to multi-response vectors, Available at https://arxiv.org/abs/2212.01621, 2024.
#'
#' @examples
#' library(didec)
#' data("bioclimatic")
#' X <- bioclimatic[, c(9:12)]
#' Y <- bioclimatic[, c(1,8)]
#' mfoci(X, Y, pre.selected = c(1, 3))
mfoci <- function(X, Y, pre.selected = NULL,
                  perm = FALSE, perm.method = c("decreasing"),
                  autostop = TRUE) {
  Idx.Perm <- c("sample", "increasing", "decreasing", "full")
  if (!perm.method %in% Idx.Perm) {
    stop("'perm.method' should be one of 'sample', 'increasing', 'decreasing', 'full'.")
  }

  X <- as.data.frame(X)
  Y <- as.data.frame(Y)
  dX <- length(X)
  Var_dX <- colnames(X)
  if (nrow(X) != nrow(Y)) {
    stop("The rows of X and Y should be equal.")
  }

  if (sum(duplicated(pre.selected)) != 0) {
    warning("Remove duplicate items from 'pre.selected'.")
    pre.selected <- pre.selected[-which(duplicated(pre.selected))]
  }

  if (!all(pre.selected %in% 1L:dX)) {
    warning("The items in 'pre.selected' should belong to [1, n], where n is the number of X's features, i.e. the number of columns of X.")
    pre.selected <- pre.selected[-which(!pre.selected %in% 1L:dX)]
    if (is_empty(pre.selected)) {
      pre.selected <- NULL
    }
  }

  num_features <- dX - length(pre.selected)
  l <- min(dX, num_features)
  Var_pre <- Var_dX[pre.selected]
  Var_l <- Var_dX[which(!Var_dX %in% Var_pre)]
  Var_add <- character()

  H <- data.frame(Feature = rep(0, l + 1), number = rep(0, l + 1))
  colnames(H)[2] <- "T"
  Z <- as.data.frame(X[, Var_pre])
  X <- as.data.frame(X[, Var_l])
  colnames(Z) <- Var_pre
  colnames(X) <- Var_l
  H[1, 1] <- "pre.selected"
  if (length(Var_pre) != 0) {
    H[1, 2] <- didec(Z, Y, perm = perm, perm.method = perm.method)
  } else {
    H[1, 2] <- 0
  }

  if (l != 0) {
    st <- H[1, 2]
    for (l2 in 2L:(l + 1)) {
      h1 <- data.frame(feat = Var_l, Tq = 0)
      if (length(Z) > 0) {
        for (l1 in seq_along(Var_l)) {
          h1[l1, 2] <- didec(cbind(Z, X[, Var_l[l1]]), Y, perm = perm, perm.method = perm.method)
        }
        H[l2, ] <- h1[which.max(h1$Tq), ]
        Var_add <- H[l2, 1]
        Z <- cbind(Z, X[, Var_add])
        colnames(Z)[length(colnames(Z))] <- Var_add
      } else {
        for (l1 in seq_along(Var_l)) {
          h1[l1, 2] <- didec(X[, Var_l[l1]], Y, perm = perm, perm.method = perm.method)
        }
        H[l2, ] <- h1[which.max(h1$Tq), ]
        Var_add <- H[l2, 1]
        Z <- data.frame(X[, Var_add])
        colnames(Z) <- Var_add
      }
      if (length(Var_l) > 0) Var_l <- Var_l[-which(Var_l == H[l2, 1])]
      X <- as.data.frame(X[, Var_l])
      colnames(X) <- Var_l
      if (autostop && H[l2, 2] <= st) {
        H <- H[1:(l2 - 1), ]
        break
      }
      st <- H[l2, 2]
    }
  }

  rownames(H) <- 0L:(nrow(H) - 1)
  H <- H[, c(2, 1)]
  return(H)
}
