#' Title Generate Data for two outcomes.
#'
#' @param xm A mean matrix with element (i,j) giving the mean of outcome i and time j.
#' @param xs A mean matrix with element (i,j) giving the s.d. of outcome i and time j.
#' @param corln Correlation between the outcomes
#' @param Sigma1 Time Dependence matrix for Outcome 1
#' @param Sigma2 Time Dependence matrix for Outcome 0
#'
#' @return Longitudinal Data with two outcomes
#' @export
#'
#' @examples

gen.fun.dep = function(xm, xs, corln = 0.5, Sigma1 = Sigma1, Sigma2 = Sigma2){

  u1 = c(0, mvrnorm(1,rep(0,nrow(Sigma1)),Sigma1))
  u2 = c(0, mvrnorm(1,rep(0,nrow(Sigma2)),Sigma2))
  ar = matrix(0, 2, 7)
  for(i in 2:7){
    sigma = matrix(0, 2, 2)
    sigma[1,1] = xs[1,i]^2
    sigma[2,2] = xs[2,i]^2
    sigma[1,2] = corln*xs[1,i]*xs[2,i]
    sigma[2,1] = corln*xs[1,i]*xs[2,i]
    # for(j in 1:7){
    ar[,i] = MASS::mvrnorm(1, xm[,i], sigma) + c(u1[i],u2[i])
    # }
  }
  ar = t(apply(ar, 1, cumsum))
}

#' Title Discretize Continuous Data into discrete bins
#'
#' @param X Longitudinal Data with continuous outcomes
#'
#' @return Longitudinal Data with discretized outcomes
#' @export
#'
#' @examples
discretize = function(X){
  X1 = X
  for(k in 1:dim(X)[2]){
    for(i in 2:dim(X)[3]){
      x1 = X[,k,i]
      x1 = x1%>%
        cut(breaks = c(-Inf, xm[k,i] - 3*xs[k,i], xm[k,i] - xs[k,i],
                       xm[k,i] + xs[k,i], xm[k,i] + 3*xs[k,i], Inf),
            labels = 1:5) %>%
        `dim<-`(dim(x1))
      x1 = as.numeric(x1)-1
      X1[,k,i] = x1
    }
  }
  return(X1)
}


