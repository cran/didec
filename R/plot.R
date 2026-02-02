#' plotting a dendrogram with colored branches
#'
#' @param dend a dendrogram
#' @param num.cluster the number of colored branches
#' @param linkage logical; if 'True', the linkage method is used
#' @param ylab a string
#' @param cex.lab a value
#' @param cex.axis a value
#'
#' @return plot
#'
#' @import dendextend
#' @importFrom graphics par
#' @importFrom graphics axis
#'
#' @keywords internal
plot_dendrogram <- function(dend, num.cluster = num.cluster, linkage = FALSE,
                            ylab = ylab, cex.lab = 0.6, cex.axis = 0.6) {
  if (!inherits(dend, "dendrogram")) {
    stop("'dend' should be a 'dendrogram'.")
  }
  n_X <- length(unlist(dend))
  if (!is.null(num.cluster)) {
    if (num.cluster < 1) {
      num.cluster <- 1
      warning("'num.cluster' has to be positive! Automatically set to 1.")
    }
    if (num.cluster > n_X) {
      num.cluster <- n_X
      warning(paste0("'num.cluster' cannot exceed the number of variables. Automatically set to ", n_X))
    }
    dend <- color_branches(dend, k = num.cluster)
  }
  labels_cex(dend) <- 0.6
  plot(dend, cex.lab = cex.lab, cex.axis = cex.axis, ylab = ylab)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if (linkage == FALSE) {
    if (n_X <= 4) {
      par(cex.axis = cex.axis)
      plot(dend, cex.lab = cex.lab, ylab = ylab, yaxt = "n")
      axis(2, 0L:(n_X - 1))
    }
  } else {
    plot(dend, cex.lab = cex.lab, cex.axis = cex.axis, ylab = ylab)
  }
}

#' Plot the trade-off of Adiam and Msplit
#'
#' @param tradeoff a data frame
#' @param main main title
#' @param sub sub title
#'
#' @return ggplot
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#'
#' @keywords internal
plot_Adiam.Msplit <- function(tradeoff,
                              main = NULL, sub = NULL) {
  tradeoff <- data.frame(partition = tradeoff[, 1], ad = tradeoff[, 2], ms = tradeoff[, 3], dis = tradeoff[, 4])
  tradeoff[, 2] <- 1 - tradeoff[, 2]
  n_best_parti <- tradeoff$partition[which(tradeoff$dis == min(tradeoff$dis))]
  x_best_parti <- tradeoff$ad[which(tradeoff$dis == min(tradeoff$dis))] + 0.03
  y_best_parti <- tradeoff$ms[which(tradeoff$dis == min(tradeoff$dis))] + 0.02
  p <- ggplot(tradeoff[, 1:2], aes(tradeoff$ad, tradeoff$ms)) +
    geom_point(pch = 16) +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(
      x = "1-average diameter",
      y = "maximum split",
      title = main,
      subtitle = sub
    ) +
    annotate("segment",xend = tradeoff$ad, yend = tradeoff$ms, x = 0, y = 0, linetype = "dotted")+ #geom_segment(aes(xend = tradeoff$ad, yend = tradeoff$ms), x = 0, y = 0, linetype = "dotted") +
    annotate("segment",xend = tradeoff$ad[which(tradeoff$dis == min(tradeoff$dis))], yend = tradeoff$ms[which(tradeoff$dis == min(tradeoff$dis))], x = 0, y = 0, linetype = "dotted", colour = "#990000")+ # geom_segment(aes(xend = tradeoff$ad[which(tradeoff$dis == min(tradeoff$dis))], yend = tradeoff$ms[which(tradeoff$dis == min(tradeoff$dis))]), x = 0, y = 0, linetype = "dotted", colour = "#990000") +
    annotate("text",x_best_parti, y_best_parti,label = n_best_parti, colour = "#990000", size = 3)+ #geom_text(aes(x_best_parti, y_best_parti), label = n_best_parti, colour = "#990000", size = 3) +
    theme(
      axis.title.x = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7)
    )
  return(p)
}

#' Plot the Silhouette coefficient
#'
#' @param Silhouette_Index a data frame of Silhouette coefficient
#' @param main main title
#' @param sub sub title
#'
#' @return ggplot
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 scale_x_continuous
#'
#' @keywords internal
plot_Silhouette.coefficient <- function(Silhouette_Index,
                                        main = NULL,
                                        sub = NULL) {
  Silhouette_Index <- data.frame(partition = Silhouette_Index[, 1], ASW = Silhouette_Index[, 2])
  partition <- Silhouette_Index$partition
  n_var <- partition[length(partition)]
  ASW <- Silhouette_Index$ASW
  n_best_parti <- partition[which(ASW == max(ASW))]
  p <- ggplot(Silhouette_Index[, 1:2], aes(partition, ASW)) +
    geom_point(pch = 16) +
    ylim(0, 1) +
    labs(
      x = "number of clusters",
      y = "average silhouette width",
      title = main,
      subtitle = sub
    ) +
    geom_vline(aes(xintercept = n_best_parti), colour = "#990000", linetype = "dotted") +
    theme(
      axis.title.x = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7)
    ) +
    scale_x_continuous(breaks=2L:n_var)
  return(p)
}
