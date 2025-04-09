#' Title Longitudinal Rank Sum Test for Two-Arm Clinical Trial (Placebo and Treatment)
#'
#' @param X_c Longitudinal Data for placebo of dimension \eqn{(n_x, K, T)}, where \eqn{n_x} is the number of
#' patients in the placebo group, K is the number of outcomes, and T is the number of time points
#' @param Y_c Longitudinal Data for Treatment of dimension \eqn{(n_y, K, T)}, where \eqn{n_y} is the number of
#' patients in the treatment group, K is the number of outcomes, and T is the number of time points
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{T.stat}: Test Statistic
#'   \item \code{T.sd}: Standard Error..
#'   \item \code{p.value}: p-value
#' }
#'
#' @examples
#' library(MASS)
#' xm = -t(cbind(c(0, 0.738, 1.313, 3.109, 4.525, 5.864, 7.338),
#'           c(0, 0.668, 3.975, 5.931, 8.259, 11.989, 13.976)))
#' xs = t(cbind(c(0, 4.79, 5.43, 6.54, 7.37, 8.15, 9.11),
#'           c(0, 10.27, 12.85, 14.95, 15.35, 16.87, 18.19)))
#' N = 100 # Total Number of Patients
#' nx = 2/5*N # Niumber of patients in Placebo
#' ny = 3/5*N # Number of patients in Treatment
#' K = 2 # Number of Outcomes
#' T = 6 # Number of visits other than baseline
#' sigma_u <- c(3.477, 7.232) # Random Effect Standard Deviation for Outcomes
#' corr <- matrix(c(1,0.5, 0.5, 1),2,2)  # Correlation between the outcomes
#' sigma = diag(sigma_u) %*% corr %*% diag(sigma_u)
#' alpha=1
#' beta=1
#' effect_size_1 = 2.21*alpha/6
#' effect_size_2 = 5.38*beta/6
#' ym = xm + t(cbind(c(0,(rep(effect_size_1,6))), c(0,(rep(effect_size_2,6)))))
#' X_placebo = array(NA, dim = c(nx, 2, 7))
#' X_dose_1 = array(NA, dim = c(ny, 2, 7))
#'
#' for(i in 1:nx){
#'   X_placebo[i,,] = gen.fun.dep(xm, xs, corr, sigma)
#' }
#' for(i in 1:ny){
#'   X_dose_1[i,,] = gen.fun.dep(ym, xs, corr, sigma)
#' }
#' X_c = X_placebo
#' Y_c = X_dose_1
#' lrst.2arm(X_c, Y_c)
#'
#' @references
#' Xu, X., Ghosh, D., Luo, S.(2025). A novel longitudinal rank-sum test for multiple primary endpoints in clinical trials: Applications to neurodegenerative disorders. Statistics in Biopharmaceutical Research, 1-17.
#' @export


lrst.2arm = function(X_c,Y_c){
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
#'
#' @examples
#' library(MASS)
#' library(lme4)
#' xm = -t(cbind(c(0, 0.738, 1.313, 3.109, 4.525, 5.864, 7.338),
#'           c(0, 0.668, 3.975, 5.931, 8.259, 11.989, 13.976)))
#' xs = t(cbind(c(0, 4.79, 5.43, 6.54, 7.37, 8.15, 9.11),
#'           c(0, 10.27, 12.85, 14.95, 15.35, 16.87, 18.19)))
#' N = 100 # Total Number of Patients
#' nx = 2/5*N # Niumber of patients in Placebo
#' ny = 3/5*N # Number of patients in Treatment
#' K = 2 # Number of Outcomes
#' T = 6 # Number of visits other than baseline
#' sigma_u <- c(3.477, 7.232) # Random Effect Standard Deviation for Outcomes
#' corr <- matrix(c(1,0.5, 0.5, 1),2,2)  # Correlation between the outcomes
#' sigma = diag(sigma_u) %*% corr %*% diag(sigma_u)
#' alpha=1
#' beta=1
#' effect_size_1 = 2.21*alpha/6
#' effect_size_2 = 5.38*beta/6
#' ym = xm + t(cbind(c(0,(rep(effect_size_1,6))), c(0,(rep(effect_size_2,6)))))
#' X_placebo = array(NA, dim = c(nx, 2, 7))
#' X_dose_1 = array(NA, dim = c(ny, 2, 7))
#'
#' for(i in 1:nx){
#'   X_placebo[i,,] = gen.fun.dep(xm, xs, corr, sigma)
#' }
#' for(i in 1:ny){
#'   X_dose_1[i,,] = gen.fun.dep(ym, xs, corr, sigma)
#' }
#' X_c = X_placebo
#' Y_c = X_dose_1
#' lmmFun4(X_c, Y_c)
#' @export

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
#' @param X_c Longitudinal Data for placebo of dimension \eqn{(n_x, K, T)}, where \eqn{n_x} is the number of
#' patients in the placebo gorup, K is the number of outcomes, and T is the number of time points
#' @param Y_c Longitudinal Data for Treatment of dimension \eqn{(n_y, K, T)}, where \eqn{n_y} is the number of
#' patients in the treatment gorup, K is the number of outcomes, and T is the number of time points
#'
#' @return A vector of p-values, length being same as the number of outcomes
#'
#' @examples
#' library(MASS)
#' library(nparLD)
#' xm = -t(cbind(c(0, 0.738, 1.313, 3.109, 4.525, 5.864, 7.338),
#'           c(0, 0.668, 3.975, 5.931, 8.259, 11.989, 13.976)))
#' xs = t(cbind(c(0, 4.79, 5.43, 6.54, 7.37, 8.15, 9.11),
#'           c(0, 10.27, 12.85, 14.95, 15.35, 16.87, 18.19)))
#' N = 100 # Total Number of Patients
#' nx = 2/5*N # Niumber of patients in Placebo
#' ny = 3/5*N # Number of patients in Treatment
#' K = 2 # Number of Outcomes
#' T = 6 # Number of visits other than baseline
#' sigma_u <- c(3.477, 7.232) # Random Effect Standard Deviation for Outcomes
#' corr <- matrix(c(1,0.5, 0.5, 1),2,2)  # Correlation between the outcomes
#' sigma = diag(sigma_u) %*% corr %*% diag(sigma_u)
#' alpha=1
#' beta=1
#' effect_size_1 = 2.21*alpha/6
#' effect_size_2 = 5.38*beta/6
#' ym = xm + t(cbind(c(0,(rep(effect_size_1,6))), c(0,(rep(effect_size_2,6)))))
#' X_placebo = array(NA, dim = c(nx, 2, 7))
#' X_dose_1 = array(NA, dim = c(ny, 2, 7))
#'
#' for(i in 1:nx){
#'   X_placebo[i,,] = gen.fun.dep(xm, xs, corr, sigma)
#' }
#' for(i in 1:ny){
#'   X_dose_1[i,,] = gen.fun.dep(ym, xs, corr, sigma)
#' }
#' X_c = X_placebo
#' Y_c = X_dose_1
#' npld(X_c, Y_c)
#'
#' @references
#' Noguchi, K., Gel, Y. R., Brunner, E., & Konietschke, F. (2012). nparLD: an R software package for the nonparametric analysis of longitudinal data in factorial experiments. Journal of Statistical software, 50, 1-23.
#'
#' @export
#'

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
    fit = nparLD::nparLD(Value ~ Day*Treatment, data = data.frame(long1), subject = "Subject",
                 description = FALSE)
    A[k] = fit$ANOVA.test.mod.Box[4]

  }

  return(A)

}

