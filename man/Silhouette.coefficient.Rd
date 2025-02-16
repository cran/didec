% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimal.partition.R
\name{Silhouette.coefficient}
\alias{Silhouette.coefficient}
\title{Silhouette coefficients given a dendrogram}
\usage{
Silhouette.coefficient(X, dend, dist.func = "PD")
}
\arguments{
\item{X}{a data frame for a set of variables X}

\item{dend}{a dendrogramm}

\item{dist.func}{PD / MPD / kendall / footrule}
}
\value{
a data frame
}
\description{
Silhouette coefficients given a dendrogram
}
\keyword{internal}
