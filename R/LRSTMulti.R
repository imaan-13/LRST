#' Title Longitudinal Rank Sum test for multi-arm Clinical Trial
#'
#' @param X_c Longitudinal Data for placebo of dimension \eqn{(n_x, K, T)}, where \eqn{n_x} is the number of
#' patients in the placebo group, K is the number of outcomes, and T is the number of time points
#' @param Y A list of length A, where A is the number of arms. Each element of Y is a longitudinal
#' data for the corresponding treatment arm of dimension \eqn{(n_{y_i}, K, T)}, where \eqn{n_{y_i}} is the number of
#' patients in the \eqn{i^{th}} treatment arm, K is the number of outcomes, and T is the number of time points.
#'
#' @return A list containing:
#' \item{p}{The computed p-value}
#' \item{v}{Proportion of times  higher dose was selected}
#' \item{T}{Final test Statistic}
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
#' X_dose_2 = array(NA, dim = c(ny, 2, 7))
#'
#' for(i in 1:nx){
#'   X_placebo[i,,] = gen.fun.dep(xm, xs, corr, sigma)
#' }
#' for(i in 1:ny){
#'   X_dose_1[i,,] = gen.fun.dep(ym, xs, corr, sigma)
#' }
#' for(i in 1:ny){
#'   X_dose_2[i,,] = gen.fun.dep(ym, xs, corr, sigma)
#' }
#' X_c = X_placebo
#' Y = vector(length = 2, mode = "list")
#' Y[[1]] = X_dose_1
#' Y[[2]] = X_dose_2
#' lrst.MultiArm(X_c, Y)
#'
#' @references
#' Ghosh, D., Luo, S. (2024). A non-parametric U-statistic testing approach for multi-arm clinical trials with multivariate longitudinal data. arXiv preprint arXiv:2408.10149.
#' @export


lrst.MultiArm = function(X_c, Y){

  nx = dim(X_c)[1]
  n = unlist(lapply(Y,function(x){dim(x)[1]}))
  T = dim(X_c)[3]
  K = dim(X_c)[2]

  N = nx + sum(n)
  N_0 = nx + n

  lambda_pair = nx/N_0

  lambda_x = nx/N
  lambda_y = n/N

  # lambda_yx = ny/N_1
  # lambda_zx = nz/N_2

  # alpha = 1 + n/nx

  v=0

  rank_xy = vector(mode = "list", length = length(n))

  for(ni in 1:length(n)){
    rank_xy[[ni]] = array(0, dim = c(K, T, nx+n[ni]))
    for (i in 1:K){
      for (j in 1:T) rank_xy[[ni]][i, j, ] = rank(c(X_c[, i, j], Y[[ni]][, i, j]), na.last = TRUE, ties.method = 'average')
    }
  }



  ## compute the average rank difference over t

  bar_rx = vector(mode = "list", length = length(n))

  for(ni in 1:length(n)){
    bar_rx[[ni]] = rep(0, T)
    for (i in 1:T) bar_rx[[ni]][i] = mean(rank_xy[[ni]][, i, 1:nx])
  }


  bar_ry = vector(mode = "list", length = length(n))

  for(ni in 1:length(n)){
    bar_ry[[ni]] = rep(0, T)
    for (i in 1:T) bar_ry[[ni]][i] = mean(rank_xy[[ni]][, i, (nx+1):(nx+n[ni])])
  }


  U_xy = vector(mode = "list", length = length(n))
  for(ni in 1:length(n)){
    U_xy[[ni]] = (bar_ry[[ni]] - bar_rx[[ni]])/sqrt(N)
  }

  ## calculate the estimated covariance

  hat_theta_xy = vector(mode = "list", length = length(n))
  for(ni in 1:length(n)){
    hat_theta_xy[[ni]] = matrix(0, K, T)
    for (i in 1:K){
      for (j in 1:T) {
        hat_theta_xy[[ni]][i, j] = (mean(rank_xy[[ni]][i, j, (nx+1):(nx+n[ni])]) -
                                      mean(rank_xy[[ni]][i, j, 1:nx]))*2/N_0[ni]
      }
    }
    hat_theta_xy[[ni]] = hat_theta_xy[[ni]][,-1]
  }


  P_xy = vector(mode = "list", length = length(n))
  for(ni in 1:length(n)){
    P_xy[[ni]] = array(0, dim = c((T-1), nx, K))
    for (t in 1:(T-1)){
      for (i in 1:nx){
        for (k in 1:K){
          P_xy[[ni]][t, i, k] = rank(c(X_c[i, k, t+1], Y[[ni]][, k, t+1]), na.last = TRUE,
                                     ties.method = 'average')[1] - 1 -
            n[ni]*(1 - hat_theta_xy[[ni]][k, t])/2
        }
      }
    }
  }

  Q_xy = vector(mode = "list", length = length(n))
  for(ni in 1:length(n)){
    Q_xy[[ni]] = array(0, dim = c(T-1, n[ni], K))
    for (t in 1:(T-1)){
      for (i in 1:n[ni]){
        for (k in 1:K){
          Q_xy[[ni]][t, i, k] = rank(c(Y[[ni]][i, k, t+1], X_c[, k, t+1]), na.last = TRUE,
                                     ties.method = 'average')[1] -
            1 - nx*(1 + hat_theta_xy[[ni]][k, t])/2
        }
      }
    }
  }


  hat_C_matrix_xy = matrix(vector("list", length = length(n)*length(n)), nrow = length(n), ncol = length(n))

  for(ni in 1:length(n)){
    for(nj in 1:length(n)){
      hat_C_matrix_xy[[ni, nj]] = matrix(0, T-1, T-1)
      for (i in 1:(T-1)){
        for (j in 1:(T-1)) {
          if(ni == nj){
            hat_C_matrix_xy[[ni, nj]][i,j] = sum((t(P_xy[[ni]][i, , ]) %*% P_xy[[nj]][j, , ])/(nx*n[ni]^2))
          }else{
            hat_C_matrix_xy[[ni, nj]][i,j] = sum((t(P_xy[[ni]][i, , ]) %*% P_xy[[nj]][j, , ])/(nx*n[ni]*n[nj]))
          }
        }
      }
    }
  }


  hat_D_matrix_xy = vector("list", length = length(n))
  for(ni in 1:length(n)){
    hat_D_matrix_xy[[ni]] = matrix(0, T-1, T-1)
    for (i in 1:(T-1)){
      for (j in 1:(T-1)) {
        hat_D_matrix_xy[[ni]][i, j] = sum((t(Q_xy[[ni]][i, , ]) %*% Q_xy[[ni]][j, , ])/(nx^2*n[ni]))
      }
    }
  }


  # hat_cov_U = (1/K^2)*((1 + 1/lambda)*hat_C_matrix + (1 + lambda)*hat_D_matrix)


  hat_cov_U = matrix(vector("list", length = length(n)*length(n)), nrow = length(n), ncol = length(n))

  for(ni in 1:length(n)){
    for(nj in 1:length(n)){
      if(ni == nj){
        hat_cov_U[[ni, nj]] = (lambda_x + lambda_y[ni])*(1/K^2)*((1 + 1/lambda_pair[ni])*hat_C_matrix_xy[[ni, ni]]
                                                                 + (1 + lambda_pair[ni])*hat_D_matrix_xy[[ni]])
      }else{
        alpha = (1 + 1/n[ni])*(lambda_x + lambda_y[nj])
        hat_cov_U[[ni, nj]] = (1/K^2)*alpha*hat_C_matrix_xy[[ni, nj]]
      }
    }
  }


  # T_stat = sum(U)
  V_T_stat = matrix(0, length(n), length(n))

  for(i in 1:length(n)){
    for(j in 1:length(n)){
      V_T_stat[i,j] = mean(hat_cov_U[[i,j]])
    }
  }

  T_stat = unlist(lapply(U_xy, sum))/sqrt(diag(V_T_stat))


  T_final = max(T_stat)/(T-1)

  v = which(T_stat == max(T_stat))

  sigma = cov_to_corr(V_T_stat)
  max_values = simulate_max_gaussian(rep(0, length(n)), sigma, 100000)

  p_value <- mean(max_values >= T_final)


  return(list("p" = p_value, "v"=v, "T" = T_final))
}


simulate_max_gaussian <- function(mean, sigma, n_simulations) {

  k = length(mean)
  simulations <- mvrnorm(n_simulations, mean, sigma)
  max_values <- apply(simulations, 1, max)

  return(max_values)
}


cov_to_corr <- function(cov_matrix) {

  std_devs <- sqrt(diag(cov_matrix))
  outer_std_devs <- outer(std_devs, std_devs)
  corr_matrix <- cov_matrix / outer_std_devs
  diag(corr_matrix) <- 1
  return(corr_matrix)
}

