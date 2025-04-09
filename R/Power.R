library(MASS)

#' Title Theoretical C and D matrix under assumption of Gaussianity
#'
#' @param xm Mean Matrix of Placebo
#' @param ym Mean Matrix of Treatment
#' @param xs Standard Deviation matrix
#' @param corr Correlation between outcomes.
#' @param n.iter Number of iterations
#'
#' @return C: Theoretical C under Gaussianity Assumption, D: Theoretical D under Gaussianity assumption
#' @export
#'
#' @examples
#' xm = -t(cbind(c(0, 0.738, 1.313, 3.109, 4.525, 5.864, 7.338),
#'           c(0, 0.668, 3.975, 5.931, 8.259, 11.989, 13.976)))
#' xs = t(cbind(c(0, 4.79, 5.43, 6.54, 7.37, 8.15, 9.11),
#'           c(0, 10.27, 12.85, 14.95, 15.35, 16.87, 18.19)))
#' alpha=1
#' beta=1
#' effect_size_1 = 2.21*alpha/6
#' effect_size_2 = 5.38*beta/6
#' ym = xm + t(cbind(c(0,(rep(effect_size_1,6))), c(0,(rep(effect_size_2,6)))))
#' C.gen(xm, ym, xs, n.iter=10)

C.gen = function(xm, ym, xs, corr=0.5, n.iter= 1000){

  T = ncol(xm)
  K = nrow(xm)
  C = matrix(0,T-1,T-1)
  D = matrix(0,T-1,T-1)

  for(t1 in 2:T){
    for(t2 in 2:T){
      s1 = 0
      s2 = 0
      for(k1 in 1:K){
        for(k2 in 1:K){
          u1 = numeric(n.iter)
          u2 = numeric(n.iter)
          for(iter in 1:n.iter){

            mu1 = c(xm[k1,t1], xm[k2,t2])
            mu2 = c(ym[k1,t1], ym[k2,t2])
            sigma = matrix(c(xs[k1,t1]^2, corr*(xs[k1,t1])*(xs[k2,t2]),
                             corr*(xs[k1,t1])*(xs[k2,t2]), xs[k2,t2]^2), 2, 2)


            x = MASS::mvrnorm(10000, mu1, sigma)
            y = MASS::mvrnorm(10000, mu2, sigma)

            v11 = pnorm(x[,1], mean = ym[k1,t1], sd = xs[k1,t1])
            v12 = pnorm(x[,2], mean = ym[k2,t2], sd = xs[k2,t2])

            v21 = pnorm(y[,1], mean = xm[k1,t1], sd = xs[k1,t1])
            v22 = pnorm(y[,2], mean = xm[k2,t2], sd = xs[k2,t2])

            u1[iter] = cov(v11, v12)
            u2[iter] = cov(v21, v22)

          }

          s1 = s1 + mean(u1)
          s2 = s2 + mean(u2)


        }
      }

      C[t1-1, t2-1] = 2*s1
      D[t1-1, t2-1] = 2*s2


    }
  }


  return(list("C" = C, "D" = D))
}



#' Title Theoretical Power when C and D are known
#'
#' @param alpha Size of the Test
#' @param effect_size Theoretical Effect Size
#' @param nx Number of patients in the plcabo group
#' @param ny Number of patients in the treatment group
#' @param K Number of outcomes
#' @param T Number of Time Points
#' @param C Theoretical C matrix
#' @param D Theoretica; D matrix
#'
#' @return Theoretical Power
#'
#' @references
#' Ghosh, D., Xu, X., Luo, S. (2025). Power and sample size calculation for multivariate longitudinal trials using the longitudinal rank sum test. arXiv preprint arXiv:2502.07152.
#'
#' @examples
#' N = 500 # Total Number of Patients
#' nx = 2/5*N # Niumber of patients in Placebo
#' ny = 3/5*N # Number of patients in Treatment
#' K = 2 # Number of Outcomes
#' T = 6 # Number of visits other than baseline
#' C = matrix(rnorm(T^2), nrow = T)
#' D = matrix(rnorm(T^2), nrow = T)
#' theoretical_power(alpha = 0.05, effect_size = 1, nx, ny, K, T+1, C, D)
#' @export

theoretical_power = function(alpha = 0.05, effect_size, nx, ny, K, T, C, D){
  N = nx + ny
  lambda = nx/ny
  J <- rep(1, T-1)
  quad_form <- (t(J) %*% (C + lambda * D) %*% J)
  z_alpha <- qnorm(1 - alpha)
  se_term <- sqrt((4 * (1 + lambda) * quad_form) / ( N * lambda * ((T-1)^2)))
  power <- pnorm((effect_size / se_term) - z_alpha)
  return(power)
}


#' Title Minimum Sample Size required to obtain a given power
#'
#' @param pi Required Power
#' @param xm Placebo mean matrix of size K * T.
#' @param ym Treatment mean matrix of size K * T
#' @param alpha Size of Test
#' @param C Theoretical C matrix
#' @param D Theoretical D matrix
#' @param lambda nx/ny where nx and ny are number of patients in the placebo and treatment groups respectively
#'
#' @return Minium Sample Size
#'
#' @export

sampSize = function(pi, xm, ym, alpha, C, D, lambda){
  mu1 = xm
  mu2 = ym
  v = ((qnorm(pi) + qnorm(1-alpha))/theta.bar.theo(mu1, mu2, xs))
  sum(C + lambda * D) * (v^2) * 4 * ((1 + lambda)/lambda) / (T^2)
}



#' Title Theoretical Relative Treatment Effect Size Under Gaussianity Assumption
#'
#' @param xm Placebo Mean
#' @param ym Treatment Mean
#' @param xs Placebo Standard Deviation
#' @param n.iter Number of iteartions
#'
#' @examples
#' xm = -t(cbind(c(0, 0.738, 1.313, 3.109, 4.525, 5.864, 7.338),
#'               c(0, 0.668, 3.975, 5.931, 8.259, 11.989, 13.976)))
#' xs = t(cbind(c(0, 4.79, 5.43, 6.54, 7.37, 8.15, 9.11),
#'               c(0, 10.27, 12.85, 14.95, 15.35, 16.87, 18.19)))
#' alpha=1
#' beta=1
#' effect_size_1 = 2.21*alpha/6
#' effect_size_2 = 5.38*beta/6
#' ym = xm + t(cbind(c(0,(rep(effect_size_1,6))), c(0,(rep(effect_size_2,6)))))
#' theta.bar.theo(xm, ym, xs, n.iter=10)
#' @export

theta.bar.theo = function(xm, ym, xs, n.iter = 1000){
  T = ncol(xm)
  K = nrow(xm)
  s = matrix(0,K,T)
  mux = t(apply(xm,1,cumsum))
  muy = t(apply(ym,1,cumsum))


  for(k in 1:K){
    for(t in 2:T){
      u = numeric(n.iter)
      for(iter in 1:n.iter){
        x = rnorm(1000, mux[k,t], xs[k,t])
        y = rnorm(1000, muy[k,t], xs[k,t])
        u[iter] = 1 - 2 * length(which(x>y))/length(x)
      }

      s[k,t] = mean(u)
    }
  }
  sum(s[,-1])/(K * (T-1))
  # sum(s)
}


#' Title Estimate C and D matrix from the data
#'
#' @param X_c Longitudinal Data for placebo of dimension \eqn{(n_x, K, T)}, where \eqn{n_x} is the number of
#' patients in the placebo group, K is the number of outcomes, and T is the number of time points
#' @param Y_c Longitudinal Data for Treatment of dimension \eqn{(n_y, K, T)}, where \eqn{n_y} is the number of
#' patients in the treatment group, K is the number of outcomes, and T is the number of time points
#'
#' @return C.hat gives estimated C matrix and D.hat gives estimated D matrix
#'
#' @examples
#' library(MASS)
#' xm = -t(cbind(c(0, 0.738, 1.313, 3.109, 4.525, 5.864, 7.338),
#'           c(0, 0.668, 3.975, 5.931, 8.259, 11.989, 13.976)))
#' xs = t(cbind(c(0, 4.79, 5.43, 6.54, 7.37, 8.15, 9.11),
#'           c(0, 10.27, 12.85, 14.95, 15.35, 16.87, 18.19)))
#' N = 500 # Total Number of Patients
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
#' C.gen.hat(X_c, Y_c)
#' @export


C.gen.hat = function(X_c, Y_c){
  m = dim(X_c)[1]
  n = dim(Y_c)[1]
  K = dim(X_c)[2]
  T = dim(X_c)[3]
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

  return(list("C.hat" = hat_C_matrix, "D.hat" = hat_D_matrix))
}


#' Title Estimate of Power from the data
#'
#' @param alpha Size of Test
#' @param X_c Longitudinal Data for placebo of dimension (nx, K, T), where nx is the number of
#' patients in the placebo group, K is the number of outcomes, and T is the number of time points
#' @param Y_c Longitudinal Data for Treatment of dimension (ny, K, T), where ny is the number of
#' patients in the treatment group, K is the number of outcomes, and T is the number of time points
#'
#'
#' @return Estimated Power from the data
#'
#' @examples
#' library(MASS)
#' xm = -t(cbind(c(0, 0.738, 1.313, 3.109, 4.525, 5.864, 7.338),
#'           c(0, 0.668, 3.975, 5.931, 8.259, 11.989, 13.976)))
#' xs = t(cbind(c(0, 4.79, 5.43, 6.54, 7.37, 8.15, 9.11),
#'           c(0, 10.27, 12.85, 14.95, 15.35, 16.87, 18.19)))
#' N = 500 # Total Number of Patients
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
#' C.gen.hat(X_c, Y_c)
#' estimated_Power(0.05, X_c, Y_c)
#' @export


estimated_Power = function(alpha = 0.05, X_c, Y_c){
  m = dim(X_c)[1]
  n = dim(Y_c)[1]
  K = dim(X_c)[2]
  T = dim(X_c)[3]
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

  V_T_stat = sum(hat_cov_U)

  uv = mean(hat_theta)/(sqrt((4*V_T_stat)/(N*(T^2)))) - qnorm(1 - alpha)
  return(pnorm(uv))
}
