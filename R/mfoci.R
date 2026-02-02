#' Multivariate feature ordering by conditional independence.
#'
#' @description
#' A variable selection algorithm based on the directed dependence coefficient (\code{\link{didec}}).
#' 
#' @param X A numeric matrix or data.frame/data.table. Contains the predictor vector X.
#' @param Y A numeric matrix or data.frame/data.table. Contains the response vector Y.
#' @param trans A logical. If \code{TRUE} the inputs of X are standardized (transformed) before the variable selection.
#' @param trans.method An optional character string specifying a method for data standardization. This must be one of the strings \code{"standardization"} (default), \code{"rank"} or \code{"rescaling"}.
#' @param estim.method An optional character string specifying a method for estimating the directed dependence coefficient \code{didec}. This must be one of the strings \code{"codec"} or \code{"copula"} (default).
#' @param perm A logical. If \code{TRUE} a version of \code{didec} that takes into account the permutations of the response variables is used in the variable selection algorithm.
#' @param perm.method An optional character string specifying a method for permuting the response variables. This must be one of the strings \code{"sample"}, \code{"increasing"}, \code{"decreasing"} (default) or \code{"full"}. 
#' @param pre.selected An integer vector for indexing pre-selected components from predictor X.
#' @param select.method An optional character string specifying a feature selection method. This must be one of the strings \code{"forward"} (default) or \code{"subset"}.
#' @param autostop A logical. If \code{True} (default) the forward feature selection algorithm stops at the first non-increasing value of \code{didec}.
#' @param max.num An integer for limiting the maximal number of selected variables if \code{select.method == "subset"}.
#'
#' @details
#' \code{mfoci} involves a forward feature selection algorithm for multiple-outcome data that employs the directed dependence coefficient (\code{didec}) at each step.
#'
#' If \code{autostop == TRUE} the algorithm stops at the first non-increasing value of \code{didec}, thereby selecting a subset of variables.
#' Otherwise, all predictor variables are ranked according to their predictive strength measured by \code{didec}.
#' 
#' In addition to the forward feature selection algorithm, this function also provides a best subset selection, which can be accomplished by \code{select.method == "subset"}. 
#' This method selects features by calculating the directed dependence coefficient of all possible feature combinations. 
#' Note that the features selected by this method are not ordered.
#'
#' @return A list containing:
#' \describe{\item{features}{A vector listing all features in X;}
#' \item{pre.selected.features}{A vector listing the pre.selected features in X if \code{pre.selected != NULL};}
#' \item{selected.features}{A data.frame listing the selected and ranked variables and the corresponding values of the directed dependence coefficient if \code{select.method == "forward"}; A vector listing the selected features if \code{select.method == "subset"};}
#' \item{valueT}{The values of the directed dependence coefficient if \code{select.method == "subset"}.}
#' }
#'
#' @export
#'
#' @author Sebastian Fuchs, Jonathan Ansari, Yuping Wang
#'
#' @importFrom rlang is_empty
#'
#' @references
#' J. Ansari, S. Fuchs, A direct extension of Azadkia & Chatterjee's rank correlation to multi-response vectors, Available at \url{https://arxiv.org/abs/2212.01621}, 2025.
#'
#' @examples
#' library(didec)
#' df <- as.data.frame(bioclimatic)
#' X <- df[, c(9:12)]
#' Y <- df[, c(1,8)]
#' mfoci(X, Y, pre.selected = c(1, 3))
mfoci <- function(X, Y,
                  trans = FALSE, trans.method = c("standardization"), 
                  estim.method = c("copula"), 
                  perm = FALSE, perm.method = c("decreasing"), 
                  pre.selected = NULL,
                  select.method = c("forward"),
                  autostop = TRUE,
                  max.num = NULL) {
  
  X <- as.data.frame(X); Y <- as.data.frame(Y)
  Var_dX <- colnames(X); Var_dY <- colnames(Y)
  dX <- length(X)
  
  if (nrow(X) != nrow(Y)) {
    stop("The rows of X and Y should be equal.")
  }

  # Y cannot be constant
  for (i in 1L:ncol(Y)){
    if (!(length(unique(Y[,i])) > 1)) {
      stop(paste0("Column ",i," of Y is constant."))
    }
  }
  
  Idx.Select <- c("forward", "subset")
  if (!select.method %in% Idx.Select) {
    stop("'select.method' should be one of 'forward', 'subset'.")
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

  mylist <- list()

  mylist$features <- Var_dX
  # mylist$response <- Var_dY

  num_features <- dX - length(pre.selected)
  l <- min(dX, num_features)
  Var_pre <- Var_dX[pre.selected]
  Var_l <- Var_dX[which(!Var_dX %in% Var_pre)]
  Var_add <- character()

  if (length(pre.selected)!=0) {
    mylist$pre.selected.features <- Var_pre
  }
  
  if (select.method == "subset"){
    
    if (length(Var_pre) != 0) {
      Z <- as.data.frame(X[, Var_pre]); colnames(Z) <- Var_pre
      X <- as.data.frame(X[, Var_l]); colnames(X) <- Var_l
      
      st <- didec(Z, Y, estim.method = estim.method, perm = perm, perm.method = perm.method, trans = trans, trans.method = trans.method)
      sl <- which(!Var_dX %in% Var_pre)
      
      ps <- powerset(sl)
      
      ####
      
      if (is.null(max.num)) {
        max.num <- dX-length(Var_pre)
        warning("'max.num' defaults to the number of variables in the predictor vector X excluding the pre-selected variabls.")
        # if (max.num > 5) {
        #   warning("select.method 'subset' might take significant computational time.")
        # }
      }
      
      if (!is.null(max.num)) {
        if ((!is.numeric(max.num)) || max.num%%1!=0) {
          stop("'max.num' should be an integer.")
        }
        if (max.num > (dX-length(Var_pre)) || max.num < 1L) {
          stop("'max.num' should be in [1,n], where n denotes the number of variables in the predictor vector X excluding the pre-selected variabls.")
        }
        
        for (i in length(ps):1L) {
          if (length(unlist(ps[i]))> max.num) {
            ps<- ps[-i]
          }
        }
        
        if (max.num > 5) {
          warning("select.method 'subset' might take significant computational time.")
        }
      }
      
      ####
      
      h1 <- data.frame(ind.ps = c(1L:length(ps)), Tq = 0)

      for (l1 in seq_along(ps)) {
        h1[l1, 2] <- didec(cbind(Z, X[, Var_dX[unlist(ps[l1])]]), Y, estim.method = estim.method, perm = perm, perm.method = perm.method, trans = trans, trans.method = trans.method)
      }
      
      ind.Feat.set <- which.max(h1$Tq)
      
      if(h1[ind.Feat.set, 2] > st){
        mylist$selected.features <- Var_dX[unlist(ps[ind.Feat.set])]
        
        if (perm == FALSE) {mylist$valueT <- h1[ind.Feat.set, 2] }
        if (perm == TRUE) {mylist$valueT <- h1[ind.Feat.set, 2] }
      }
      else{
        mylist$selected.features <- NULL 
        if (perm == FALSE) {mylist$valueT <- st}
        if (perm == TRUE) {mylist$valueT <- st}
      }
    }
    
    if (length(Var_pre) == 0) {
      X <- as.data.frame(X[, Var_l])
      colnames(X) <- Var_l
      ps <- powerset(1L:length(Var_l))
      
      ####
      
      if (is.null(max.num)) {
        max.num <- dX
        warning("'max.num' defaults to the number of variables in the predictor vector X.")
        
        # if (max.num > 5) {
        #   warning("select.method 'subset' might take significant computational time.")
        # }
      }
      
      if (!is.null(max.num)) {
        if ((!is.numeric(max.num)) || max.num%%1!=0) {
          stop("'max.num' should be an integer.")
        }
        if (max.num > dX || max.num < 1L) {
          stop("'max.num' should be in [1,n], where n denotes the number of variables in the predictor vector X.")
        }

        for (i in length(ps):1L) {
          if (length(unlist(ps[i]))> max.num) {
            ps<- ps[-i]
          }
        }
        
        if (max.num > 5) {
          warning("select.method 'subset' might take significant computational time.")
        }
      }
      
      ####
      
      h1 <- data.frame(ind.ps = c(1L:length(ps)), Tq = 0)
      for (l1 in seq_along(ps)) {
        h1[l1, 2] <- didec(X[, Var_dX[unlist(ps[l1])]], Y, estim.method = estim.method, perm = perm, perm.method = perm.method, trans = trans, trans.method = trans.method)
      }
      ind.Feat.set <- which.max(h1$Tq)
      mylist$selected.features <- Var_dX[unlist(ps[ind.Feat.set])]
      if (perm == FALSE) {mylist$valueT <- h1[ind.Feat.set, 2] }
      if (perm == TRUE) {mylist$valueT <- h1[ind.Feat.set, 2] }
    }
  }
  
  if (select.method == "forward") {
    if (perm == FALSE) {colnam <- "valueT"}
    if (perm == TRUE) {colnam <- "valueT"}
    
    H <- data.frame(Feature = rep(0, l + 1), number = rep(0, l + 1))
    colnames(H)[2] <- colnam
    Z <- as.data.frame(X[, Var_pre])
    X <- as.data.frame(X[, Var_l])
    colnames(Z) <- Var_pre
    colnames(X) <- Var_l
    H[1, 1] <- "pre.selected"
    if (length(Var_pre) != 0) {
      H[1, 2] <- didec(Z, Y, estim.method = estim.method, perm = perm, perm.method = perm.method, trans = trans, trans.method = trans.method)
    } else {
      H[1, 2] <- 0
    }
    
    if (l != 0) {
      st <- H[1, 2]
      for (l2 in 2L:(l + 1)) {
        h1 <- data.frame(feat = Var_l, Tq = 0)
        if (length(Z) > 0) {
          for (l1 in seq_along(Var_l)) {
            h1[l1, 2] <- didec(cbind(Z, X[, Var_l[l1]]), Y, estim.method = estim.method, perm = perm, perm.method = perm.method, trans = trans, trans.method = trans.method)
          }
          H[l2, ] <- h1[which.max(h1$Tq), ]
          Var_add <- H[l2, 1]
          Z <- cbind(Z, X[, Var_add])
          colnames(Z)[length(colnames(Z))] <- Var_add
        } else {
          for (l1 in seq_along(Var_l)) {
            h1[l1, 2] <- didec(X[, Var_l[l1]], Y, estim.method = estim.method, perm = perm, perm.method = perm.method, trans = trans, trans.method = trans.method)
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
    # H <- H[, c(2, 1)]
    
    mylist$selected.features <- H
  }

  
  return(mylist)
}