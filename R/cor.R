#' Markov product estimate from single (q=1) endogenuous and (p>=1) exogenous variables based on dimension reduction
#'
#' @param X a data frame for input vector X
#' @param y a data frame for output vector Y
#'
#' @return a value
#'
#' @importFrom RANN nn2
#' @importFrom stats sd
#'
#' @keywords internal
MPhi <- function(X,y) {

  df <- data.frame(X,y)
  p <- ncol(df) - 1; n <- nrow(df)
  
  # R_i and R_N(i)
  df$rk_y <- rank(df[,p+1],ties.method="random")/(n+1)
  df$rk_y_nn <- rank(df$y[nn2(df[,1L:p],df[,1L:p],k=2)$nn.idx[,2]],ties.method = "random")/(n+1)
  df <- data.frame(y1=df$rk_y, y2=df$rk_y_nn)
  return(df)
}


#' Estimate for \eqn{\xi}(Y,X) based on dimension reduction principle
#'
#' @param X a data frame for input vector X
#' @param y a data frame for output vector Y
#'
#' @return a value
#'
#' @importFrom stats complete.cases
#'
#' @keywords internal
CopulaCorr <- function(X,y){

  X <- as.data.frame(X); y <- as.data.frame(y)

  MP <- MPhi(X,y)
  n <- nrow(MP)
  MP <- MP*(n+1)

  # Normalization
  D <- data.frame(X,y)
  A <- as.data.frame(table(D[,length(D)]))
  z1 <- rank(sort(sample(x=as.numeric(A$Var1), sum(as.numeric(A$Freq)), prob=A$Freq, replace = TRUE)),ties.method="random")/(n+1)
  z2 <- rank(sort(sample(x=as.numeric(A$Var1), sum(as.numeric(A$Freq)), prob=A$Freq, replace = TRUE)),ties.method="random")/(n+1)
  Mp <- MPhi(z1,z2)
  n <- nrow(Mp)
  Mp <- Mp*(n+1)

  TT <- 1-3/(n^2-1)*sum(abs(MP$y1-MP$y2)) + 3/(n^2-1)*((sum(MP$y2)+sum(MP$y1)-n*(n+1)))
  TTmax <- 1-3/(n^2-1)*sum(abs(Mp$y1-Mp$y2)) + 3/(n^2-1)*((sum(Mp$y2)+sum(Mp$y1)-n*(n+1)))
  return(TT/TTmax)
}


#' Estimate for T(Y,X) based on dimension reduction principle
#'
#' @param X a data frame for input vector X
#' @param Y a data frame for output vector Y
#'
#' @return a value
#'
#' @keywords internal
Copula.Tq <- function(X,Y){

  dY <- length(Y)
  ZW <- numeric()
  weight <- numeric()

  weight[1] <- 0
  ZW[1] <- as.numeric(CopulaCorr(X,Y[,1]))
  if(dY>1){
    for(i in 2L:dY){
      weight[i] <- as.numeric(CopulaCorr(Y[,1L:(i-1)],Y[,i]))
      ZW[i] <- as.numeric(CopulaCorr(data.frame(X,Y[,1L:(i-1)]), Y[,i]))
    }
  }
  return(1-(dY-sum(ZW))/(dY-sum(weight)))
}



#' Estimate for T_bar(Y,X) based on dimension reduction principle
#'
#' @param X a data frame for input vector X
#' @param Y a data frame for output vector Y
#' @param method permuatation methods: sample / increasing / decreasing / full
#'
#' @return a value
#'
#' @importFrom gtools permutations
#'
#' @keywords internal
Copula.Tq.Perm <- function(X,Y,method=c("sample")){

  dY <- length(Y)

  if(dY>1){
    if(method=="increasing"){
      perm <- matrix(1L:dY,dY,dY+1,byrow=T)[,1L:dY]
    }
    if(method=="decreasing"){
      perm <- matrix(dY:1L,dY,dY+1,byrow=T)[,1L:dY]
    }
    if(method=="sample"){
      perm <- permutations(n = dY, r = dY, v = 1L:dY)
      perm <- perm[sample(1L:factorial(dY), size = dY, replace = FALSE),]
    }
    if (method == "full") {
      perm <- permutations(n = dY, r = dY, v = 1L:dY)
    }

    Result <- numeric()
    for (l in 1L:nrow(perm)){
      Y <- Y[,perm[l,]]
      Result[l] <- Copula.Tq(X,Y)
    }
    Result <- mean(Result)
  }
  if(dY==1){
    Result <- Copula.Tq(X,Y)
  }
  return(Result)
}


#' Estimate for \eqn{\xi}(Y,X) using codec function
#'
#' @param X a data frame for input vector X
#' @param y a data frame for output vector Y
#'
#' @return a value
#'
#' @importFrom stats complete.cases
#' @importFrom FOCI codec
#'
#' @keywords internal
CodecCorr <- function(X,y){

  X <- as.data.frame(X); y <- as.data.frame(y)

  df <- data.frame(X,y)
  p <- ncol(df) - 1

  return(as.numeric(codec(df[,p+1], df[,1L:p])))
}

#' Estimate for T(Y,X) based on function codec
#'
#' @param X a data frame for input vector X
#' @param Y a data frame for output vector Y
#'
#' @return a value 
#'
#' @keywords internal
Codec.Tq <- function(X, Y) {

  dY <- length(Y)
  ZW <- numeric()
  weight <- numeric()

  weight[1] <- 0
  ZW[1] <- as.numeric(CodecCorr(X,Y[,1]))
  if(dY>1){
    for(i in 2L:dY){
      weight[i] <- as.numeric(CodecCorr(Y[,1L:(i-1)],Y[,i]))
      ZW[i] <- as.numeric(CodecCorr(data.frame(X,Y[,1L:(i-1)]), Y[,i]))
    }
  }
  return(CodecTq = 1 - (dY - sum(ZW)) / (dY - sum(weight)))
}


#' Estimate for T_bar(Y,X) based on function Codec & a sample of all / all increasing / all decreasing permutations
#'
#' @param X a data frame for input vector X
#' @param Y a data frame for output vector Y
#' @param method permuatation methods: sample / increasing / decreasing / full
#'
#' @return a value
#'
#' @importFrom gtools permutations
#'
#' @keywords internal
Codec.Tq.Perm <- function(X, Y, method = c("sample")) {

  dY <- length(Y)

  if (dY > 1) {

    if (method == "sample") {
      perm <- permutations(n = dY, r = dY, v = 1L:dY)
      perm <- perm[sample(1L:factorial(dY), size = dY, replace = FALSE), ]
    }
    if (method == "increasing") {
      perm <- matrix(1L:dY, dY, dY + 1, byrow = T)[, 1L:dY]
    }
    if (method == "decreasing") {
      perm <- matrix(dY:1L, dY, dY + 1, byrow = T)[, 1L:dY]
    }
    if (method == "full") {
      perm <- permutations(n = dY, r = dY, v = 1L:dY)
    }

    Results <- numeric()
    for (l in seq_len(nrow(perm))) {
      Y <- Y[, perm[l, ]]
      Results[l] <- Codec.Tq(X, Y)
    }
    Result <- mean(Results)
  }
  if (dY == 1) {
    Result <- Codec.Tq(X, Y)
  }
  return(Result)
}


#' Estimation for multivariate concordance measures
#'
#' @param X a data frame for vector X
#' @param method kendall / footrule
#'
#' @return a value of the estimator for the multivariate concordance measures
#' 
#' @importFrom pcaPP cor.fk
#' @importFrom copBasic footCOP
#'
#' @keywords internal
concor.M <- function(X, method = c("footrule")) {
  
  X <- as.data.frame(X)
  n <- nrow(X)
  d <- length(X)

  if (method == "kendall") {
    if (d == 2) {
      t <- cor.fk(X[,1], X[,2])
      return(t)
    }
    
    if (d > 2) {
      S <- c()
      xn <- c(1L:n)
      for (i in xn) {
        for (j in xn[xn != i]) {
          L <- c()
          for (k in 1L:d) {
            if (X[i, k] <= X[j, k]) {
              L <- c(L, TRUE)
            }
          }
          if (sum(L) == d) {
            S <- c(S, 1)
          }
        }
      }
      t <- ((2^d * sum(S)) / (n * (n - 1)) - 1) / (2^(d - 1) - 1)
      return(t)
    }
  }

  if (method == "footrule") {
    if (d == 2) {
      phi <- footCOP(para = X, as.sample = TRUE)
      return(phi)
    }
    if (d>2) {
      R <- X
      for (j in 1L:d) {
        R[, j] <- rank(X[, j], ties.method = "random")
      }
      L <- c()
      for (i in 1L:n) {
        Li <- max(R[i, ]) - min(R[i, ])
        L <- c(L, Li)
      }
      phi <- 1 - ((d + 1) / (d - 1)) * (sum(L) / (n^2 - 1))
      return(phi)
    }
  }
}
