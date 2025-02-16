% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VarClustPartition.R
\name{VarClustPartition}
\alias{VarClustPartition}
\title{Hierarchical variable clustering.}
\usage{
VarClustPartition(
  X,
  dist.method = c("PD"),
  linkage = FALSE,
  link.method = c("complete"),
  part.method = c("optimal"),
  criterion = c("Adiam&Msplit"),
  num.cluster = NULL,
  plot = FALSE
)
}
\arguments{
\item{X}{A numeric matrix or data.frame/data.table. Contains the variables to be clustered.}

\item{dist.method}{An optional character string computing a distance function for clustering. This must be one of the strings \code{"PD"} (default), \code{"MPD"}, \code{"kendall"} or \code{"footrule"}.}

\item{linkage}{A logical. If \code{TRUE} a linkage method is used.}

\item{link.method}{An optional character string selecting a linkage method. This must be one of the strings \code{"complete"} (default), \code{"average"} or \code{"single"}.}

\item{part.method}{An optional character string selecting a partitioning method. This must be one of the strings \code{"optimal"} (default) or \code{"selected"}.}

\item{criterion}{An optional character string selecting a criterion for the optimal partition, if \code{part.method = "optimal"}. This must be one of the strings \code{"Adiam&Msplit"} (default) or \code{"Silhouette"}.}

\item{num.cluster}{An integer value for the selected number of clusters, if \code{part.method = "selected"}.}

\item{plot}{A logical. If \code{TRUE} a dendrogram is plotted with colored branches according to the corresponding partitioning method.}
}
\value{
A list containing a dendrogram without colored branches (\bold{dendrogram}), an integer value determining the number of clusters after partitioning (\bold{num.cluster}), and a list containing the clusters after partitioning (\bold{clusters}).
}
\description{
\code{VarClustPartition} is a hierarchical variable clustering algorithm based on the directed dependence coefficient (\code{\link{didec}}) or a concordance measure (Kendall tau \eqn{\tau} or Spearman's footrule) according to a pre-selected number of clusters or an optimality criterion (Adiam&Msplit or Silhouette coefficient).
}
\details{
\code{VarClustPartition} performs a hierarchical variable clustering based on the directed dependence coefficient (\code{\link{didec}}) and provides a partition of the set of variables.

If \code{dist.method =="PD"} or \code{dist.method =="MPD"}, the clustering is performed using \code{didec} either as a directed ("PD") or as a symmetric ("MPD") dependence coefficient.
If \code{dist.method =="kendall"} or \code{dist.method =="footrule"}, clustering is performed using either multivariate Kendall's tau ("kendall") or multivariate Spearman's footrule ("footrule").

Instead of using one of the above-mentioned four multivariate measures for the clustering, the option \code{linkage == TRUE} enables the use of bivariate linkage methods,
including complete linkage (\code{link.method == "complete"}), average linkage (\code{link.method == "average"}) and single linkage (\code{link.method == "single"}).
Note that the multivariate distance methods are computationally demanding because higher-dimensional dependencies are included in the calculation, in contrast to linkage methods which only incorporate pairwise dependencies.

A pre-selected number of clusters \code{num.cluster} can be realized with the option \code{part.method == "selected"}.
Otherwise (\code{part.method == "optimal"}), the number of clusters is determined by maximizing the intra-cluster similarity (similarity within the same cluster) and minimizing the inter-cluster similarity (similarity among the clusters). Two optimality criteria are available:

\code{"Adiam&Msplit"}: \emph{Adiam} measures the intra-cluster similarity and \emph{Msplit} measures the inter-cluster similarity.

\code{"Silhouette"}: A mixed coefficient incorporating the intra-cluster similarity and the inter-cluster similarity. The optimal number of clusters corresponds to the maximum Silhouette coefficient.
}
\examples{
library(didec)
n  <- 50
X1 <- rnorm(n,0,1)
X2 <- X1
X3 <- rnorm(n,0,1)
X4 <- X3 + X2
X  <- data.frame(X1=X1,X2=X2,X3=X3,X4=X4)
vcp <- VarClustPartition(X,
                            dist.method = c("PD"),
                            part.method = c("optimal"),
                            criterion   = c("Silhouette"),
                            plot        = TRUE)
vcp$clusters
\donttest{
data("bioclimatic")
X   <- bioclimatic[c(2:4,9)]
vcp1 <- VarClustPartition(X,
                          linkage     = TRUE,
                          link.method = c("complete"),
                          dist.method = "PD",
                          part.method = "optimal",
                          criterion   = "Silhouette",
                          plot        = TRUE)
vcp1$clusters
vcp2 <- VarClustPartition(X,
                          linkage     = TRUE,
                          link.method = c("complete"),
                          dist.method = "footrule",
                          part.method = "optimal",
                          criterion   = "Adiam&Msplit",
                          plot        = TRUE)
vcp2$clusters
}
}
\references{
S. Fuchs, Y. Wang, Hierarchical variable clustering based on the predictive strength between random vectors, Int. J. Approx. Reason. 170, Article ID 109185, 2024.

P. Hansen, B. Jaumard, Cluster analysis and mathematical programming, Math. Program. 79 (1) 191–215, 1997.

L. Kaufman, Finding Groups in Data, John Wiley & Sons, 1990.
}
\author{
Yuping Wang, Sebastian Fuchs
}
