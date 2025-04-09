#' Title
#'
#' @param X_c Longitudinal Data for placebo of dimension (nx, K, T), where nx is the number of
#' patients in the placebo gorup, K is the number of outcomes, and T is the number of time points
#' @param Y_c Longitudinal Data for Treatment of dimension (ny, K, T), where ny is the number of
#' patients in the treatment gorup, K is the number of outcomes, and T is the number of time points
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{T.stat}: Test Statistic
#'   \item \code{T.sd}: Standard Error..
#'   \item \code{p.value}: p-value
#' }

#' @export
#'
#' @examples

uniUstat = function(X_c,Y_c){
  m = dim(X_c)[1]
  n = dim(Y_c)[1]
  K = dim(X_c)[2]
  T = dim(X_c)[3]-1
  N = m + n
  lambda = m/n
  rank = array(0, dim = c(K, T, m+n))
  for (i in 1:K){
    for (j in 1:T) rank[i, j, ] = rank(c(X_c[, i, j], Y_c[, i, j]), na.last = T,
                                       ties.method = 'average')
  }

  ## compute the average rank difference over t
  bar_rx = rep(0, T)
  for (i in 1:T) bar_rx[i] = mean(rank[, i, 1:m])

  bar_ry = rep(0, T)
  for (i in 1:T) bar_ry[i] = mean(rank[, i, (m+1):(m+n)])

  U = (bar_ry - bar_rx)/sqrt(N)

  ## calculate the estimated covariance
  hat_theta = matrix(0, K, T)
  for (i in 1:K){
    for (j in 1:T) hat_theta[i, j] = (mean(rank[i, j, (m+1):(m+n)]) - mean(rank[i, j, 1:m]))*2/N
  }

  P = array(0, dim = c(T, m, K))
  for (t in 1:T){
    for (i in 1:m){
      for (k in 1:K){
        P[t, i, k] = rank(c(X_c[i, k, t], Y_c[, k, t]), na.last = T, ties.method = 'average')[1] -
          1 - n*(1 - hat_theta[k, t])/2
      }
    }
  }

  Q = array(0, dim = c(T, n, K))
  for (t in 1:T){
    for (i in 1:n){
      for (k in 1:K){
        Q[t, i, k] = rank(c(Y_c[i, k, t], X_c[, k, t]), na.last = T, ties.method = 'average')[1] -
          1 - m*(1 + hat_theta[k, t])/2
      }
    }
  }

  hat_C_matrix = matrix(0, T, T)
  for (i in 1:T){
    for (j in 1:T) hat_C_matrix[i, j] = sum((t(P[i, , ]) %*% P[j, , ])/(m*n^2))
  }


  hat_D_matrix = matrix(0, T, T)
  for (i in 1:T){
    for (j in 1:T) hat_D_matrix[i, j] = sum((t(Q[i, , ]) %*% Q[j, , ])/(m^2*n))
  }

  hat_cov_U = (1/K^2)*((1 + 1/lambda)*hat_C_matrix + (1 + lambda)*hat_D_matrix)

  T_stat = sum(U)
  V_T_stat = sum(hat_cov_U)
  T_final = (T_stat/sqrt(V_T_stat))

  p = pnorm(T_final, 0, 1, lower.tail = F, log.p = F)
  return(list("T.stat" = T_stat, "T.sd" = V_T_stat, "p.value" = p))
}





#' Title Linear Mixed Effect Model Efficacy Test
#'
#' @param X_c Longitudinal Data for placebo of dimension (nx, K, T), where nx is the number of
#' patients in the placebo gorup, K is the number of outcomes, and T is the number of time points
#' @param Y_c Longitudinal Data for Treatment of dimension (ny, K, T), where ny is the number of
#' patients in the treatment gorup, K is the number of outcomes, and T is the number of time points
#'
#' @return A vector of p-values, length being same as the number of outcomes
#' @export
#'
#' @examples
lmmFun4 = function(X_c, Y_c){

  nx = dim(X_c)[1]
  ny = dim(Y_c)[1]

  K = dim(X_c)[2]
  T = dim(X_c)[3]-1
  A = numeric(K)
  for(k in 1:K){
    out1 = cbind(1:(nx+ny),rbind(X_c[,k,], Y_c[,k,]), c(rep(0,nx), rep(1,ny)))
    out1 = data.frame(out1)
    names(out1) = c("Subject", 0:T,  "Treatment")
    long1 = tidyr::pivot_longer(out1, -c(Subject, Treatment),
                                values_to = "Value", names_to = "Day")
    long1$Day = as.numeric(long1$Day)
    fit1 = lme4::lmer(Value ~ Treatment + Day + Treatment*Day + (1 + Day| Subject), data = long1)
    fit2 = lme4::lmer(Value ~ Day + (1 + Day|Subject), data = long1)
    A[k] = anova(fit1, fit2)$`Pr(>Chisq)`[2]

  }

  return(A)

}


#' Title nparLD Test of Efficacy
#'
#' @param X_c Longitudinal Data for placebo of dimension (nx, K, T), where nx is the number of
#' patients in the placebo gorup, K is the number of outcomes, and T is the number of time points
#' @param Y_c Longitudinal Data for Treatment of dimension (ny, K, T), where ny is the number of
#' patients in the treatment gorup, K is the number of outcomes, and T is the number of time points
#'
#' @return A vector of p-values, length being same as the number of outcomes
#' @export
#'
#' @examples
npld = function(X_c, Y_c){

  nx = dim(X_c)[1]
  ny = dim(Y_c)[1]
  K = dim(X_c)[2]
  T = dim(X_c)[3]-1
  A = numeric(K)
  for(k in 1:K){
    out1 = cbind(1:(nx+ny),rbind(X_c[,k,], Y_c[,k,]), c(rep(0,nx), rep(1,ny)))
    out1 = data.frame(out1)
    names(out1) = c("Subject", 0:T,  "Treatment")
    long1 = tidyr::pivot_longer(out1, -c(Subject, Treatment),
                                values_to = "Value", names_to = "Day")
    long1$Day = as.numeric(long1$Day)
    # fit1 = nparLD(Value ~ Day, data = data.frame(long1), subject = "Subject",
    # description = FALSE)
    fit = nparLD(Value ~ Day*Treatment, data = data.frame(long1), subject = "Subject",
                 description = FALSE)
    A[k] = fit$ANOVA.test.mod.Box[4]

  }

  return(A)

}

