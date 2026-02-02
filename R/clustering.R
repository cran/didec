#' Cluster a set of variables using distance function based on predictive measure
#'
#' @param X a data frame for a set of variables X
#' @param estim.method An optional character string specifying a method for estimating the directed dependence coefficient.
#' @param mutual type B function or not
#'
#' @return a list for hierarchical clustering result
#'
#' @keywords internal
clust.Tq <- function(X,
                     estim.method = c("copula"),
                     mutual = FALSE) {

  df <- as.data.frame(X)
  dX <- length(df)
  cn <- colnames(df)

  # variable classes needed to be clustered
  class <- list()
  for (i in 1L:dX) {
    class <- append(class, cn[i], after = length(class))
  }

  S <- 0 # record clustering step
  nam <- c() # record the classes before every clustering step
  L <- list() # record the final clustering result

  while (length(class) > 1) {
    S <- S + 1

    #### classes before clustering
    aln0 <- c()
    for (i in seq_along(class)) {
      aln0 <- c(aln0, paste(unlist(class[i]), collapse = ","))
    }
    aln <- paste(aln0, collapse = ";")
    nam <- c(nam, aln)

    #### distance matrix for distances between classes
    dist <- matrix(0, length(class), length(class))
    for (i in c(1:(length(class) - 1))) {
      for (j in c((i + 1):length(class))) {
        dist[i, j] <- dist.Tq(df[, unlist(class[i])], df[, unlist(class[j])], estim.method = estim.method, mutual = mutual)
      }
      for (j in c(1:i)) {
        dist[i, j] <- NA
      }
    }
    dist[length(class), ] <- NA

    index <- which(dist == min(dist, na.rm = TRUE), arr.ind = TRUE) # index for two classes having minimum distance
    new <- list(c(unlist(class[index[1, 1]]), unlist(class[index[1, 2]]))) # two classes having minimum distance

    #### the classes clustered and their distance in this clustering step
    L[S] <- list(c(paste(unlist(class[index[1, 1]]), collapse = ","), paste(unlist(class[index[1, 2]]), collapse = ","), min(dist, na.rm = TRUE)))

    #### update the classes needed to be clustered
    class <- class[-c(index[1, 1], index[1, 2])]
    class <- append(class, list(new))
  }
  names(L) <- nam
  return(L)
}

#' Clustering a set of variables using distance function based on multivariate concordance measures
#'
#' @param X a data frame for vector X
#' @param method kendall / footrule
#'
#' @return a list for hierarchical clustering result
#'
#' @keywords internal
clust.concor.M <- function(X, method = c("footrule")) {

  df <- as.data.frame(X)
  dX <- length(df)
  cn <- colnames(df)

  # variable classes needed to be clustered
  class <- list()
  for (i in c(1:dX)) {
    class <- append(class, cn[i], after = length(class))
  }

  S <- 0 # record clustering step
  nam <- c() # record the classes before every clustering step
  L <- list() # record the final clustering result

  while (length(class) > 1) {
    S <- S + 1

    #### classes before clustering
    aln0 <- c()
    for (i in seq_along(class)) {
      aln0 <- c(aln0, paste(unlist(class[i]), collapse = ","))
    }
    aln <- paste(aln0, collapse = ";")
    nam <- c(nam, aln)

    #### distance matrix for distances between classes
    dist <- matrix(0, length(class), length(class))
    for (i in c(1:(length(class) - 1))) {
      for (j in c((i + 1):length(class))) {
        dist[i, j] <- dist.concor.M(df[, unlist(class[i])], df[, unlist(class[j])], method = method)
      }
      for (j in c(1:i)) {
        dist[i, j] <- NA
      }
    }
    dist[length(class), ] <- NA

    index <- which(dist == min(dist, na.rm = TRUE), arr.ind = TRUE) # index for two classes having minimum distance
    new <- list(c(unlist(class[index[1, 1]]), unlist(class[index[1, 2]]))) # two classes having minimum distance

    #### the classes clustered and their distance in this clustering step
    L[S] <- list(c(paste(unlist(class[index[1, 1]]), collapse = ","), paste(unlist(class[index[1, 2]]), collapse = ","), min(dist, na.rm = TRUE)))

    #### update the classes needed to be clustered
    class <- class[-c(index[1, 1], index[1, 2])]
    class <- append(class, list(new))
  }

  names(L) <- nam
  return(L)
}

#' Read a dendrogram from a list for hierarchical clustering result
#'
#' @param clust a list for hierarchical clustering result
#' @param step whether using clustering step as y axis or not
#'
#' @return an object of class "dendrogram"
#'
#' @importFrom phylogram read.dendrogram
#' @importFrom stats as.dendrogram
#'
#' @keywords internal
dendrogram <- function(clust, step = TRUE) {
  s <- length(clust)
  pre <- c()
  m <- c()
  f <- c()

  if (step == TRUE) {
    for (i in c(1:s)) {
      n <- unlist(clust[i])[1:2]

      if ((n[1] %in% pre) && (n[2] %in% pre)) {
        m1 <- which(f == n[1])
        m2 <- which(f == n[2])
        new <- paste("(", m[m1], ":", i - m1, ",", m[m2], ":", i - m2, ")", sep = "")
        m <- c(m, new)
        f <- c(f, paste(n[1], n[2], sep = ","))
        pre <- c(pre, paste(n[1], n[2], sep = ","))
      }

      if ((!n[1] %in% pre) && (n[2] %in% pre)) {
        m2 <- which(f == n[2])
        new <- paste("(", n[1], ":", i, ",", m[m2], ":", i - m2, ")", sep = "")
        m <- c(m, new)
        f <- c(f, paste(n[1], n[2], sep = ","))
        pre <- c(pre, n[1], paste(n[1], n[2], sep = ","))
      }

      if (n[1] %in% pre && !n[2] %in% pre) {
        m1 <- which(f == n[1])
        new <- paste("(", m[m1], ":", i - m1, ",", n[2], ":", i, ")", sep = "")
        m <- c(m, new)
        f <- c(f, paste(n[1], n[2], sep = ","))
        pre <- c(pre, n[2], paste(n[1], n[2], sep = ","))
      }

      if (!n[1] %in% pre && !n[2] %in% pre) {
        new <- paste("(", n[1], ":", i, ",", n[2], ":", i, ")", sep = "")
        m <- c(m, new)
        f <- c(f, paste(n[1], n[2], sep = ","))
        pre <- c(pre, n[1], n[2], paste(n[1], n[2], sep = ","))
      }
    }
  }

  if (step == FALSE) {
    for (i in c(1:s)) {
      n <- unlist(clust[i])[1:3]

      if (n[1] %in% pre && n[2] %in% pre) {
        m1 <- which(f == n[1])
        m2 <- which(f == n[2])
        new <- paste("(", m[m1], ":", as.numeric(n[3]) - as.numeric(unlist(clust[[length(clust)]][m1])[3]), ",", m[m2], ":", as.numeric(n[3]) - as.numeric(unlist(clust[[length(clust)]][m2])[3]), ")", sep = "")
        m <- c(m, new)
        f <- c(f, paste(n[1], n[2], sep = ","))
        pre <- c(pre, paste(n[1], n[2], sep = ","))
      }

      if (!n[1] %in% pre && n[2] %in% pre) {
        m2 <- which(f == n[2])
        new <- paste("(", n[1], ":", n[3], ",", m[m2], ":", as.numeric(n[3]) - as.numeric(unlist(clust[[length(clust)]][m2])[3]), ")", sep = "")
        m <- c(m, new)
        f <- c(f, paste(n[1], n[2], sep = ","))
        pre <- c(pre, n[1], paste(n[1], n[2], sep = ","))
      }

      if (n[1] %in% pre && !n[2] %in% pre) {
        m1 <- which(f == n[1])
        new <- paste("(", m[m1], ":", as.numeric(n[3]) - as.numeric(unlist(clust[[length(clust)]][m1])[3]), ",", n[2], ":", n[3], ")", sep = "")
        m <- c(m, new)
        f <- c(f, paste(n[1], n[2], sep = ","))
        pre <- c(pre, n[2], paste(n[1], n[2], sep = ","))
      }

      if (!n[1] %in% pre && !n[2] %in% pre) {
        new <- paste("(", n[1], ":", n[3], ",", n[2], ":", n[3], ")", sep = "")
        m <- c(m, new)
        f <- c(f, paste(n[1], n[2], sep = ","))
        pre <- c(pre, n[1], n[2], paste(n[1], n[2], sep = ","))
      }
    }
  }

  text <- paste(m[length(m)], ";", sep = "")
  dend <- as.dendrogram(read.dendrogram(text = text))
  return(dend)
}
