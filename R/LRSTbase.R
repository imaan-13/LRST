utils::globalVariables(c("Subject", "Treatment", "xm", "xs", "prob", "rt2d"))


#'Generate Data for Multiple outcomes.
#'
#'Simulates repeated measures for multiple outcomes across time points per subject, incorporating outcome correlations and random effects.
#'
#' @param xm A mean matrix with element (i,j) giving the mean of outcome i and time j.
#' @param xs A s.d. matrix with element (i,j) giving the s.d. of outcome i and time j.
#' @param corln Correlation matrix between the outcomes.
#' @param Sigma Covariance Matrix for Random Effect
#'
#' @return Longitudinal Data with multiple outcomes
#' @details This function is used inside a for-loop to simulate repeated measures per subject
#'
#' @examples
#' # Example: simulate 1 subject with 7 visits
#' library(MASS)
#' M <- -t(cbind(
#'   c(0, 0.738, 1.313, 3.109, 4.525, 5.864, 7.338),
#'   c(0, 0.668, 3.975, 5.931, 8.259, 11.989, 13.976)
#' ))
#'
#' S <- t(cbind(
#'   c(0, 4.79, 5.43, 6.54, 7.37, 8.15, 9.11),
#'   c(0, 10.27, 12.85, 14.95, 15.35, 16.87, 18.19)
#' ))
#'
#' Sigma <- matrix(c(3.477^2, 0, 0, 3.477^2), nrow = 2)
#' corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#'
#' gen.fun.dep(M, S, corr, Sigma)
#'
#' @export

gen.fun.dep = function(xm, xs, corln, Sigma){

  K = nrow(xm)
  u = mvrnorm(1,rep(0,K),Sigma)

  ar = matrix(0, K, 7)
  for(i in 2:7){
    sigma = matrix(0, K, K)
    for(k1 in 1:K){
      for(k2 in 1:K){
        if(k1 != k2){
          sigma[k1,k2] = corln[k1,k2]*(xs[k1,i]^2)*(xs[k2,1]^2)
        }else{
          sigma[k1,k2] = xs[k1, i]^2
        }
      }
    }

    ar[,i] = MASS::mvrnorm(1, xm[,i], sigma) + u

  }
  ar = t(apply(ar, 1, cumsum))
  ar
}

#' Discretize Continuous Data into discrete bins
#'
#' Converts continuous longitudinal outcomes into discrete bins based on their mean and standard deviation at each time point.
#'
#' @param X Longitudinal Data with continuous outcomes
#'
#' @return Longitudinal Data with discretized outcomes
#'
#' @examples
#' library(MASS)
#' library(dplyr)
#' xm = -t(cbind(
#'   c(0, 0.738, 1.313, 3.109, 4.525, 5.864, 7.338),
#'   c(0, 0.668, 3.975, 5.931, 8.259, 11.989, 13.976)
#' ))
#' xs = t(cbind(
#'   c(0, 4.79, 5.43, 6.54, 7.37, 8.15, 9.11),
#'   c(0, 10.27, 12.85, 14.95, 15.35, 16.87, 18.19)
#' ))
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
#' X_d = discretize(X_c)
#' @export

discretize = function(X){
  X1 = X
  for(k in 1:dim(X)[2]){
    for(i in 2:dim(X)[3]){
      x1 = X[,k,i]
      x1 = x1|>
        cut(breaks = c(-Inf, xm[k,i] - 3*xs[k,i], xm[k,i] - xs[k,i],
                       xm[k,i] + xs[k,i], xm[k,i] + 3*xs[k,i], Inf),
            labels = 1:5) |>
        `dim<-`(dim(x1))
      x1 = as.numeric(x1)-1
      X1[,k,i] = x1
    }
  }
  return(X1)
}


