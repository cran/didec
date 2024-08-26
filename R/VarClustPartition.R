#' Hierarchical variable clustering.
#'
#' @description
#' \code{VarClustPartition} is a hierarchical variable clustering algorithm based on the directed dependence coefficient (\code{\link{didec}}) or a concordance measure (Kendall tau \eqn{\tau} or Spearman's footrule) according to a pre-selected number of clusters or an optimality criterion (Adiam&Msplit or Silhouette coefficient).
#'
#' @param X A numeric matrix or data.frame/data.table. Contains the variables to be clustered.
#' @param dist.method An optional character string computing a distance function for clustering. This must be one of the strings \code{"PD"} (default), \code{"MPD"}, \code{"kendall"} or \code{"footrule"}.
#' @param linkage A logical. If \code{TRUE} a linkage method is used.
#' @param link.method An optional character string selecting a linkage method. This must be one of the strings \code{"complete"} (default), \code{"average"} or \code{"single"}.
#' @param part.method An optional character string selecting a partitioning method. This must be one of the strings \code{"optimal"} (default) or \code{"selected"}.
#' @param criterion An optional character string selecting a criterion for the optimal partition, if \code{part.method = "optimal"}. This must be one of the strings \code{"Adiam&Msplit"} (default) or \code{"Silhouette"}.
#' @param num.cluster An integer value for the selected number of clusters, if \code{part.method = "selected"}.
#' @param plot A logical. If \code{TRUE} a dendrogram is plotted with colored branches according to the corresponding partitioning method.
#'
#' @details
#' \code{VarClustPartition} performs a hierarchical variable clustering based on the directed dependence coefficient (\code{\link{didec}}) and provides a partition of the set of variables.
#'
#' If \code{dist.method =="PD"} or \code{dist.method =="MPD"}, the clustering is performed using \code{didec} either as a directed ("PD") or as a symmetric ("MPD") dependence coefficient.
#' If \code{dist.method =="kendall"} or \code{dist.method =="footrule"}, clustering is performed using either multivariate Kendall's tau ("kendall") or multivariate Spearman's footrule ("footrule").
#'
#' Instead of using one of the above-mentioned four multivariate measures for the clustering, the option \code{linkage == TRUE} enables the use of bivariate linkage methods,
#' including complete linkage (\code{link.method == "complete"}), average linkage (\code{link.method == "average"}) and single linkage (\code{link.method == "single"}).
#' Note that the multivariate distance methods are computationally demanding because higher-dimensional dependencies are included in the calculation, in contrast to linkage methods which only incorporate pairwise dependencies.
#'
#' A pre-selected number of clusters \code{num.cluster} can be realized with the option \code{part.method == "selected"}.
#' Otherwise (\code{part.method == "optimal"}), the number of clusters is determined by maximizing the intra-cluster similarity (similarity within the same cluster) and minimizing the inter-cluster similarity (similarity among the clusters). Two optimality criteria are available:
#'
#' \code{"Adiam&Msplit"}: \emph{Adiam} measures the intra-cluster similarity and \emph{Msplit} measures the inter-cluster similarity.
#'
#' \code{"Silhouette"}: A mixed coefficient incorporating the intra-cluster similarity and the inter-cluster similarity. The optimal number of clusters corresponds to the maximum Silhouette coefficient.
#'
#' @return A list containing a dendrogram without colored branches (\bold{dendrogram}), an integer value determining the number of clusters after partitioning (\bold{num.cluster}), and a list containing the clusters after partitioning (\bold{clusters}).
#'
#' @export
#'
#' @author Yuping Wang, Sebastian Fuchs
#'
#' @importFrom cowplot plot_grid
#' @importFrom grDevices recordPlot
#' @importFrom stats as.hclust
#' @importFrom stats as.dendrogram
#' @importFrom stats hclust
#' @import dendextend
#'
#' @references
#' S. Fuchs, Y. Wang, Hierarchical variable clustering based on the predictive strength between random vectors, Int. J. Approx. Reason. 170, Article ID 109185, 2024.
#'
#' P. Hansen, B. Jaumard, Cluster analysis and mathematical programming, Math. Program. 79 (1) 191–215, 1997.
#'
#' L. Kaufman, Finding Groups in Data, John Wiley & Sons, 1990.
#'
#' @examples
#' library(didec)
#' n  <- 50
#' X1 <- rnorm(n,0,1)
#' X2 <- X1
#' X3 <- rnorm(n,0,1)
#' X4 <- X3 + X2
#' X  <- data.frame(X1=X1,X2=X2,X3=X3,X4=X4)
#' vcp <- VarClustPartition(X,
#'                             dist.method = c("PD"),
#'                             part.method = c("optimal"),
#'                             criterion   = c("Silhouette"),
#'                             plot        = TRUE)
#' vcp$clusters
#' \donttest{
#' data("bioclimatic")
#' X   <- bioclimatic[c(2:4,9)]
#' vcp1 <- VarClustPartition(X,
#'                           linkage     = TRUE,
#'                           link.method = c("complete"),
#'                           dist.method = "PD",
#'                           part.method = "optimal",
#'                           criterion   = "Silhouette",
#'                           plot        = TRUE)
#' vcp1$clusters
#' vcp2 <- VarClustPartition(X,
#'                           linkage     = TRUE,
#'                           link.method = c("complete"),
#'                           dist.method = "footrule",
#'                           part.method = "optimal",
#'                           criterion   = "Adiam&Msplit",
#'                           plot        = TRUE)
#' vcp2$clusters
#' }
VarClustPartition <- function(X, dist.method = c("PD"),
                              linkage = FALSE, link.method = c("complete"),
                              part.method = c("optimal"), criterion = c("Adiam&Msplit"),
                              num.cluster = NULL, plot = FALSE) {
  X <- as.data.frame(X); n_X <- length(X)

  if (n_X < 2L) {
    stop("'X' should have more than one columns")
  }

  var <- names(X)
  for (i in seq(var)) {
    if (grepl(" ", var[i])) {
      var[i] <- unlist(strsplit(var[i], " "))[1]
    }
  }
  names(X) <- var

  Idx.dist <- c("PD", "MPD", "kendall", "footrule")
  if (!dist.method %in% Idx.dist) {
    stop("'dist.method' should be one of 'PD', 'MPD', 'kendall','footrule'.")
  }

  if (linkage) {
    Idx.link <- c("complete", "average", "single")
    if (!link.method %in% Idx.link) {
      stop("'link.method' should be one of 'complete', 'avergae', 'single'.")
    }
  }

  Idx.part <- c("optimal", "selected")
  if (!part.method %in% Idx.part) {
    stop("'part.method' should be one of 'optimal', 'selected'.")
  }

  if (part.method == "optimal") {
    Idx.crit <- c("Adiam&Msplit", "Silhouette")
    if (!criterion %in% Idx.crit) {
      stop("'criterion' should be one of 'Adiam&Msplit', 'Silhouette'.")
    }
  }

  if (part.method == "selected") {
    if (is.null(num.cluster)) {
      num.cluster <- 1L
      warning("'num.cluster' defaults to 1.")
    }
    if (num.cluster%%1!=0) {
      stop("'num.cluster' should be an integer.")
    }
    if (num.cluster > n_X || num.cluster < 1L) {
      stop("'num.cluster' should be in [1,n], where n denotes the number of columns of X.")
    }
  }

  mylist <- list()

  if (linkage == FALSE) {
    if (dist.method == "PD"){
      clust <- clust.Tq(X)
    }
    if (dist.method == "MPD"){
      clust <- clust.Tq(X, mutual = TRUE)
    }
    if (dist.method %in% c("kendall", "footrule")){
      clust <- clust.concor.M(X, method = dist.method)
    }
    dend <- dendrogram(clust)
    mylist$dendrogram <- dend
  }

  if (linkage == TRUE) {
    if (dist.method == "PD"){
      distmat <- dist.mat.T(X)
    }
    if (dist.method == "MPD"){
      distmat <- dist.mat.T(X, mutual = TRUE)
    }
    if (dist.method %in% c("kendall", "footrule")){
      distmat <- dist.mat.concor(X, method = dist.method)
    }
    dend <- as.dendrogram(hclust(distmat, method = link.method))
    mylist$dendrogram <- dend
  }

  if (part.method == "selected") {
    mylist$num.cluster <- num.cluster
  }

  if (part.method == "optimal") {
    if (criterion == "Adiam&Msplit") {
      trf <- Adiam.Msplit(X, dend, dist.func = dist.method)
      num.cluster <- trf$partition[which(trf$dis == min(trf$dis))]
    }
    if (criterion == "Silhouette") {
      trf <- Silhouette.coefficient(X, dend, dist.func = dist.method)
      num.cluster <- trf$partition[which(trf$ASW == max(trf$ASW))]
    }
    mylist$num.cluster <- num.cluster
  }

  cls <- list()
  for (j in 1L:num.cluster) {
    cls[j] <- list(var[which(var %in% names(which(cutree(as.hclust(dend), k = num.cluster) == j)))])
  }
  mylist$clusters <- cls

  if (plot == TRUE) {
    if (linkage == FALSE) {
      ylab = "clustering step"
    }
    if (linkage == TRUE) {
      ylab = "dissimilarity"
    }
    plot_dendrogram(dend = mylist$dendrogram, num.cluster = mylist$num.cluster, linkage = linkage, ylab = ylab)
    p1 <- recordPlot()
    if (part.method == "selected") {
      p <- plot_grid(p1, rel_heights = c(1, .6))
    }
    if (part.method == "optimal") {
      if (criterion == "Adiam&Msplit") {
        p2 <- plot_Adiam.Msplit(trf)
      }
      if (criterion == "Silhouette") {
        p2 <- plot_Silhouette.coefficient(trf)
      }
      p <- plot_grid(p1, p2, rel_heights = c(1, .6))
    }
    mylist$plot <- p
  }
  return(mylist)
}
