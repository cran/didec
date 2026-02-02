#' Diameter of a class of variables based on different distance function
#'
#' @param X a data frame for a set of variables X
#' @param dist.func PD / MPD / kendall / footrule
#' @param estim.method An optional character string specifying a method for estimating the directed dependence coefficient.
#'
#' @return a value
#'
#' @keywords internal
diam <- function(X, dist.func = "PD", estim.method = c("copula")) {
  
  df_X <- as.data.frame(X)
  n_X <- length(df_X)
  if (n_X == 1) {
    diam <- 1
  }
  if (n_X != 1) {
    if (dist.func == "PD") {
      diam <- min(1 - dist.mat.T(df_X, estim.method = estim.method))
    }
    if (dist.func == "MPD") {
      diam <- min(1 - dist.mat.T(df_X, estim.method = estim.method, mutual = TRUE))
    }
    if (dist.func == "kendall") {
      diam <- min(1 - dist.mat.concor(df_X, method = "kendall"))
    }
    if (dist.func == "footrule") {
      diam <- min(1 - dist.mat.concor(df_X, method = "footrule"))
    }
  }
  return(diam)
}


#' Split of two classes of variables based on different distance function
#'
#' @param X a data frame for a set of variables X
#' @param Y a data frame for a set of variables Y
#' @param dist.func PD / MPD / kendall / footrule
#' @param estim.method An optional character string specifying a method for estimating the directed dependence coefficient.
#'
#' @return a value
#'
#' @keywords internal
split <- function(X, Y, dist.func = "PD", estim.method = c("copula")) {
  
  df_X <- as.data.frame(X)
  df_Y <- as.data.frame(Y)
  n_X <- length(df_X)
  n_Y <- length(df_Y)
  if (n_Y == 0) {
    sim <- NA
  }
  if (n_Y != 0) {
    sim <- c()
    for (i in 1L:n_X) {
      df <- as.data.frame(cbind(df_X[, i], df_Y))
      if (dist.func == "PD") {
        sim <- c(sim, max(1 - dist.mat.T(df, estim.method = estim.method)[1L:n_Y]))
      }
      if (dist.func == "MPD") {
        sim <- c(sim, max(1 - dist.mat.T(df, estim.method = estim.method, mutual = TRUE)[1L:n_Y]))
      }
      if (dist.func == "kendall") {
        sim <- c(sim, max(1 - dist.mat.concor(df, method = "kendall")[1L:n_Y]))
      }
      if (dist.func == "footrule") {
        sim <- c(sim, max(1 - dist.mat.concor(df, method = "footrule")[1L:n_Y]))
      }
    }
  }
  return(max(sim))
}


#' Average diameter & Maximum split of every partition of a given dendrogram
#'
#' @param X a data frame for a set of variables X
#' @param dend a dendrogramm
#' @param dist.func PD / MPD / kendall / footrule
#' @param estim.method An optional character string specifying a method for estimating the directed dependence coefficient.
#'
#' @return a data frame
#'
#' @importFrom stats as.hclust
#' @import dendextend
#'
#' @keywords internal
Adiam.Msplit <- function(X, dend = dend, dist.func = "PD", estim.method = c("copula")) {

  if (!inherits(dend, "dendrogram")) {
    stop("'dend' should be a 'dendrogram'.")
  }
  df <- as.data.frame(X)
  var <- names(df)
  for (i in seq(var)) {
    if (grepl(" ", var[i])) {
      var[i] <- unlist(strsplit(var[i], " "))[1]
    }
  }
  names(df) <- var
  # all partitions of the given dendrogram
  all_partition <- list()
  for (i in 2L:length(df)) {
    part <- list()
    for (j in 1L:i) {
      part[j] <- list(which(var %in% names(which(cutree(as.hclust(dend), k = i) == j))))
    }
    all_partition[i - 1] <- list(part)
  }
  # Adiam and Msplit for every partition
  partition <- c()
  ad <- c()
  ms <- c()
  dis <- c()
  for (i in seq_along(all_partition)) {
    clust <- all_partition[i][[1]]
    nclust <- length(clust)
    d <- c()
    s <- c() # diamter and split for k. class in i. partition
    for (k in 1L:nclust) {
      c <- clust[[k]]
      cX <- df[, c]
      cY <- df[, -c]
      d <- c(d, diam(cX, dist.func = dist.func, estim.method = estim.method))
      s <- c(s, split(cX, cY, dist.func = dist.func, estim.method = estim.method))
    }
    ad_partition <- mean(d)
    ms_partition <- max(s)
    dis_partition <- sqrt((1 - ad_partition)^2 + ms_partition^2)
    partition <- c(partition, nclust)
    ad <- c(ad, ad_partition)
    ms <- c(ms, ms_partition)
    dis <- c(dis, dis_partition)
  }
  tradeoff <- data.frame(partition = partition, ad = ad, ms = ms, dis = dis)
  return(tradeoff)
}


#' Silhouette value for the i. variable given variable partition
#'
#' @param i the index of the variable
#' @param df a data frame for all variables
#' @param partition a partition
#' @param dist.func PD / MPD / kendall / footrule
#' @param estim.method An optional character string specifying a method for estimating the directed dependence coefficient.
#'
#' @return a value for Silhouette
#'
#' @keywords internal
Silhouette <- function(i, df, partition, dist.func = "PD", estim.method = c("copula")) {

  a_i <- function(obj, cls, dist.func = dist.func, estim.method = estim.method) {
    if (!identical(obj, cls[, 1])) {
      stop("Error")
    }
    if (dist.func == "PD") {
      a <- mean(dist.mat.T(cls, estim.method = estim.method)[1L:(length(cls) - 1)])
    }
    if (dist.func == "MPD") {
      a <- mean(dist.mat.T(cls, estim.method = estim.method, mutual = TRUE)[1L:(length(cls) - 1)])
    }
    if (dist.func == "kendall") {
      a <- mean(dist.mat.concor(cls, method = "kendall")[1L:(length(cls) - 1)])
    }
    if (dist.func == "footrule") {
      a <- mean(dist.mat.concor(cls, method = "footrule")[1L:(length(cls) - 1)])
    }
    return(a)
  }

  b <- c()
  for (j in seq_along(partition)) {
    if (i %in% partition[[j]]) {
      if (length(partition[[j]]) == 1) {
        return(0)
      }
      ai <- a_i(df[, i], df[, c(i, partition[[j]][partition[[j]] != i])], dist.func = dist.func, estim.method = estim.method)
    } else {
      b <- c(b, a_i(df[, i], df[, c(i, partition[[j]])], dist.func = dist.func, estim.method = estim.method))
    }
  }
  bi <- min(b)
  return((bi - ai) / max(ai, bi))
}


#' Silhouette coefficients given a dendrogram
#'
#' @param X a data frame for a set of variables X
#' @param dend a dendrogramm
#' @param dist.func PD / MPD / kendall / footrule
#' @param estim.method An optional character string specifying a method for estimating the directed dependence coefficient.
#'
#' @return a data frame
#'
#' @importFrom stats as.hclust
#' @import dendextend
#'
#' @keywords internal
Silhouette.coefficient <- function(X, dend, dist.func = "PD", estim.method = c("copula")) {
  
  if (!inherits(dend, "dendrogram")) {
    stop("'dend' should be a 'dendrogram'.")
  }
  df <- as.data.frame(X)
  var <- names(df)
  for (i in seq(var)) {
    if (grepl(" ", var[i])) {
      var[i] <- unlist(strsplit(var[i], " "))[1]
    }
  }
  names(df) <- var
  # all partitions of the given dendrogram
  all_partition <- list()
  for (i in 2L:length(df)) {
    part <- list()
    for (j in 1L:i) {
      part[j] <- list(which(var %in% names(which(cutree(as.hclust(dend), k = i) == j))))
    }
    all_partition[i - 1] <- list(part)
  }
  # Silhouette coefficient for every partition
  num <- c()
  ASW <- c()
  for (i in seq_along(all_partition)) {
    partition <- all_partition[i][[1]]
    s <- c()
    for (j in seq_along(X)) {
      s <- c(s, Silhouette(j, X, partition, dist.func = dist.func, estim.method = estim.method))
    }
    num <- c(num, length(partition))
    ASW <- c(ASW, mean(s))
  }
  Silhouette_Index <- data.frame(partition = num, ASW = ASW)
  return(Silhouette_Index)
}
