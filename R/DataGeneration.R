#############################   Data Generation ##############################################

######################  Linear Set-up   #######################################################


data.gen.lmm = function(N, alpha, beta){

  effect_size_1 = 2.21*alpha/6
  effect_size_2 = 5.38*beta/6

  nx = 2/5*N
  ny = 3/5*N

  sigma_e_1 <- 3.477  # Error standard deviation
  sigma_u_11 <- 3.665  # Random effect Intercept standard deviation
  sigma_u_12 <- 1.396 # Random effect Slope standard deviation
  correlation_1 <- 0.2  # Correlation between random intercepts and slopes
  rep.col<-function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
  }

  ##### DAD

  sigma_e_2 <- 7.232  # Error standard deviation
  sigma_u_21 <- 9.036  # Random effect Intercept standard deviation
  sigma_u_22 <- 2.871 # Random effect Slope standard deviation
  correlation_2 <- 0.08  # Correlation between random intercepts and slopes

  sigma = matrix(c(sigma_e_1^2, 0.5*sigma_e_1*sigma_e_2, 0.5*sigma_e_1*sigma_e_2, sigma_e_2^2), 2, 2)


  u1 <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(sigma_u_11^2,
                                                  correlation_1 * sigma_u_11 * sigma_u_12,
                                                  correlation_1 * sigma_u_11 * sigma_u_12,
                                                  sigma_u_12^2), nrow = 2))
  u2 <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(sigma_u_21^2,
                                                  correlation_2 * sigma_u_21 * sigma_u_22,
                                                  correlation_2 * sigma_u_21 * sigma_u_22,
                                                  sigma_u_22^2), nrow = 2))




  X_c = array(NA, dim = c(nx, 2, T))
  Y_c = array(NA, dim = c(ny, 2, T))

  for(i in 1:nx){

    X_c[i,,] = rbind(rep(-0.36211, T), rep(-0.83016, T)) + t(t((c(-1.38507, -2.65461)))) %*% t((0:(T-1))) +
      rbind(rep(u1[i,1], T), rep(u2[i,1], T)) +  t(t((c(u1[i,2], u2[i,2])))) %*% t(((0:(T-1)))) +
      (t(mvrnorm(T, c(0,0), sigma)))
    X_c[i,,] = X_c[i,,] - rep.col(X_c[i,,1], T)
  }


  for(i in 1:ny){
    Y_c[i,,] = rbind(rep(-0.36211, T), rep(-0.83016, T)) + t(t((c(-1.38507, -2.65461)))) %*% t((0:(T-1))) +
      rbind(rep(u1[nx + i,1], T), rep(u2[nx + i,1], T)) +  t(t((c(u1[nx + i,2], u2[nx + i,2])))) %*% t(((0:(T-1)))) +
      t(t((c(effect_size_1, effect_size_2)))) %*% t((0:(T-1))) +
      (t(mvrnorm(T, c(0,0), sigma)))
    Y_c[i,,] = Y_c[i,,] - rep.col(Y_c[i,,1], T)
  }

  return(list("Placebo" = X_c, "Treatment" = Y_c))
}


######################   Non-linear Set-up   #######################################################

data.gen.nl = function(N, f.nl, f.dev, alpha, beta){

  effect_size_1 = 2.21*alpha/6
  effect_size_2 = 5.38*beta/6

  nx = 2/5*N
  ny = 3/5*N

  sigma_e_1 <- 3.477  # Error standard deviation
  sigma_u_11 <- 3.665  # Random effect Intercept standard deviation
  sigma_u_12 <- 1.396 # Random effect Slope standard deviation
  correlation_1 <- 0.2  # Correlation between random intercepts and slopes
  rep.col<-function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
  }

  ##### DAD

  sigma_e_2 <- 7.232  # Error standard deviation
  sigma_u_21 <- 9.036  # Random effect Intercept standard deviation
  sigma_u_22 <- 2.871 # Random effect Slope standard deviation
  correlation_2 <- 0.08  # Correlation between random intercepts and slopes

  sigma = matrix(c(sigma_e_1^2, 0.5*sigma_e_1*sigma_e_2, 0.5*sigma_e_1*sigma_e_2, sigma_e_2^2), 2, 2)


  u1 <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(sigma_u_11^2,
                                                  correlation_1 * sigma_u_11 * sigma_u_12,
                                                  correlation_1 * sigma_u_11 * sigma_u_12,
                                                  sigma_u_12^2), nrow = 2))
  u2 <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(sigma_u_21^2,
                                                  correlation_2 * sigma_u_21 * sigma_u_22,
                                                  correlation_2 * sigma_u_21 * sigma_u_22,
                                                  sigma_u_22^2), nrow = 2))




  X_c = array(NA, dim = c(nx, 2, T))
  Y_c = array(NA, dim = c(ny, 2, T))

  for(i in 1:nx){

    X_c[i,,] = rbind(rep(-0.36211, T), rep(-0.83016, T)) + t(t((c(-1.38507, -2.65461)))) %*% t(f.nl((0:(T-1)))) +
      rbind(rep(u1[i,1], T), rep(u2[i,1], T)) +  t(t((c(u1[i,2], u2[i,2])))) %*% t(((0:(T-1)))) +
      t(as.matrix(rt2d(T, rho = 0.5, nu = 4)))/sqrt(3)
    # (t(mvrnorm(T, c(0,0), sigma)))
    X_c[i,,] = X_c[i,,] - rep.col(X_c[i,,1], T)
  }


  for(i in 1:ny){
    Y_c[i,,] = rbind(rep(-0.36211, T), rep(-0.83016, T)) + t(t((c(-1.38507, -2.65461)))) %*% t(f.nl((0:(T-1)))) +
      rbind(rep(u1[nx + i,1], T), rep(u2[nx + i,1], T)) +  t(t((c(u1[nx + i,2], u2[nx + i,2])))) %*% t(((0:(T-1)))) +
      t(t((c(effect_size_1, effect_size_2)))) %*% t(f.dev((0:(T-1)))) +
      t(as.matrix(rt2d(T, rho = 0.5, nu = 4)))/sqrt(3)
    # (t(mvrnorm(T, c(0,0), sigma)))
    Y_c[i,,] = Y_c[i,,] - rep.col(Y_c[i,,1], T)
  }

  return(list("Placebo" = X_c, "Treatment" = Y_c))
}



######################   General Set-up with covariates   #######################################################

data.gen.nl.cov = function(N, f.nl, f.dev, alpha, beta){

  effect_size_1 = 2.21*alpha/6
  effect_size_2 = 5.38*beta/6

  nx = 2/5*N
  ny = 3/5*N

  sigma_e_1 <- 3.477  # Error standard deviation
  sigma_u_11 <- 3.665  # Random effect Intercept standard deviation
  sigma_u_12 <- 1.396 # Random effect Slope standard deviation
  correlation_1 <- 0.2  # Correlation between random intercepts and slopes
  rep.col<-function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
  }

  ##### DAD

  sigma_e_2 <- 7.232  # Error standard deviation
  sigma_u_21 <- 9.036  # Random effect Intercept standard deviation
  sigma_u_22 <- 2.871 # Random effect Slope standard deviation
  correlation_2 <- 0.08  # Correlation between random intercepts and slopes

  sigma = matrix(c(sigma_e_1^2, 0.5*sigma_e_1*sigma_e_2, 0.5*sigma_e_1*sigma_e_2, sigma_e_2^2), 2, 2)


  cov_time = matrix(0, nx + ny, T)
  treatment = c(rep(0,nx), rep(1, ny))
  allele_Prob = c(prob, 2*sqrt(prob)*(1- sqrt(prob)) , (1 - sqrt(prob))^2)
  allele = sample(c(0,1,2), N, prob = allele_Prob, replace = T)
  for(i in 1:N){
    cov_time[i,] = 0.2 * (0:(T-1)) + rnorm(T)
  }

  u1 <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(sigma_u_11^2,
                                                  correlation_1 * sigma_u_11 * sigma_u_12,
                                                  correlation_1 * sigma_u_11 * sigma_u_12,
                                                  sigma_u_12^2), nrow = 2))
  u2 <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(sigma_u_21^2,
                                                  correlation_2 * sigma_u_21 * sigma_u_22,
                                                  correlation_2 * sigma_u_21 * sigma_u_22,
                                                  sigma_u_22^2), nrow = 2))




  X_c = array(NA, dim = c(nx, 2, T))
  Y_c = array(NA, dim = c(ny, 2, T))

  for(i in 1:nx){

    X_c[i,,] = rbind(rep(-0.36211, T), rep(-0.83016, T)) + t(t((c(-1.38507, -2.65461)))) %*% (t((0:(T-1)))) +
      rbind(rep(u1[i,1], T), rep(u2[i,1], T)) +  t(t((c(u1[i,2], u2[i,2])))) %*% (t((0:(T-1)))) +
      rbind(rep(0.3 * allele[i], T), rep(0.2 * allele[i], T)) +
      t(t(c(0.1 * allele[i] , 0.5 * allele[i]))) %*% (t((0 : (T-1)))) +
      0.3 * cov_time[i,] + (t(mvrnorm(T, c(0,0), sigma)))
    X_c[i,,] = X_c[i,,] - rep.col(X_c[i,,1], T)

    X_c_1[i,,] = rbind(rep(-0.36211, T), rep(-0.83016, T)) + t(t((c(-1.38507, -2.65461)))) %*% (t((0:(T-1)))) +
      rbind(rep(u1[i,1], T), rep(u2[i,1], T)) +  t(t((c(u1[i,2], u2[i,2])))) %*% (t((0:(T-1)))) +
      (t(mvrnorm(T, c(0,0), sigma)))
    X_c_1[i,,] = X_c_1[i,,] - rep.col(X_c_1[i,,1], T)


  }


  for(i in 1:ny){

    Y_c[i,,] = rbind(rep(-0.36211, T), rep(-0.83016, T)) + t(t((c(-1.38507, -2.65461)))) %*% (t((0:(T-1)))) +
      rbind(rep(u1[nx + i,1], T), rep(u2[nx + i,1], T)) + t(t((c(u1[nx + i,2], u2[nx + i,2])))) %*% (t((0:(T-1)))) +
      rbind(rep(effect_size_1, T), rep(effect_size_2, T)) +
      t(t((c(effect_size_1, effect_size_2)))) %*% (t((0:(T-1)))) +
      rbind(rep(0.3 * allele[nx + i], T), rep(0.2 * allele[nx + i], T)) +
      t(t(c(0.1 * allele[nx + i] , 0.5 * allele[nx + i]))) %*% (t((0 : (T-1)))) +
      0.3 * cov_time[nx + i,] +
      (t(mvrnorm(T, c(0,0), sigma)))
    Y_c[i,,] = Y_c[i,,] - rep.col(Y_c[i,,1], T)


    Y_c_1[i,,] = rbind(rep(-0.36211, T), rep(-0.83016, T)) + t(t((c(-1.38507, -2.65461)))) %*% (t((0:(T-1)))) +
      rbind(rep(u1[nx + i,1], T), rep(u2[nx + i,1], T)) + t(t((c(u1[nx + i,2], u2[nx + i,2])))) %*% (t((0:(T-1)))) +
      rbind(rep(effect_size_1, T), rep(effect_size_2, T)) +
      t(t((c(effect_size_1, effect_size_2)))) %*% (t((0:(T-1)))) +
      (t(mvrnorm(T, c(0,0), sigma)))
    Y_c_1[i,,] = Y_c_1[i,,] - rep.col(Y_c_1[i,,1], T)

  }


  return(list("Placebo" = X_c, "Treatment" = Y_c))
}




data.gen.nl.cov.combo = function(N, f.nl, f.dev, alpha, beta){

  effect_size_1 = 2.21*alpha/6
  effect_size_2 = 5.38*beta/6

  nx = 2/5*N
  ny = 3/5*N

  sigma_e_1 <- 3.477  # Error standard deviation
  sigma_u_11 <- 3.665  # Random effect Intercept standard deviation
  sigma_u_12 <- 1.396 # Random effect Slope standard deviation
  correlation_1 <- 0.2  # Correlation between random intercepts and slopes
  rep.col<-function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
  }

  ##### DAD

  sigma_e_2 <- 7.232  # Error standard deviation
  sigma_u_21 <- 9.036  # Random effect Intercept standard deviation
  sigma_u_22 <- 2.871 # Random effect Slope standard deviation
  correlation_2 <- 0.08  # Correlation between random intercepts and slopes

  sigma = matrix(c(sigma_e_1^2, 0.5*sigma_e_1*sigma_e_2, 0.5*sigma_e_1*sigma_e_2, sigma_e_2^2), 2, 2)


  cov_time = matrix(0, nx + ny, T)
  treatment = c(rep(0,nx), rep(1, ny))
  prob = 0.7
  allele_Prob = c(prob, 2*sqrt(prob)*(1- sqrt(prob)) , (1 - sqrt(prob))^2)
  allele = sample(c(0,1,2), N, prob = allele_Prob, replace = T)
  for(i in 1:N){
    cov_time[i,] = 0.2 * (0:(T-1)) + rnorm(T)
  }

  u1 <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(sigma_u_11^2,
                                                  correlation_1 * sigma_u_11 * sigma_u_12,
                                                  correlation_1 * sigma_u_11 * sigma_u_12,
                                                  sigma_u_12^2), nrow = 2))
  u2 <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(sigma_u_21^2,
                                                  correlation_2 * sigma_u_21 * sigma_u_22,
                                                  correlation_2 * sigma_u_21 * sigma_u_22,
                                                  sigma_u_22^2), nrow = 2))




  X_c = array(NA, dim = c(nx, 2, T))
  Y_c = array(NA, dim = c(ny, 2, T))

  X_c_1 = array(NA, dim = c(nx, 2, T))
  Y_c_1 = array(NA, dim = c(ny, 2, T))

  for(i in 1:nx){


    X_c_1[i,,] = rbind(rep(-0.36211, T), rep(-0.83016, T)) + t(t((c(-1.38507, -2.65461)))) %*% t(f.nl((0:(T-1)))) +
      rbind(rep(u1[i,1], T), rep(u2[i,1], T)) +  t(t((c(u1[i,2], u2[i,2])))) %*% t(((0:(T-1)))) +
      t(as.matrix(rt2d(T, rho = 0.5, nu = 4)))/sqrt(3)
    # (t(mvrnorm(T, c(0,0), sigma)))

    X_c_1[i,,] = X_c_1[i,,] - rep.col(X_c_1[i,,1], T)

    X_c[i,,] = X_c_1[i,,] + rbind(rep(0.3 * allele[i], T), rep(0.2 * allele[i], T)) +
      t(t(c(0.1 * allele[i] , 0.5 * allele[i]))) %*% (t((0 : (T-1)))) +
      0.3 * cov_time[i,]
    X_c[i,,] = X_c[i,,] - rep.col(X_c[i,,1], T)

  }


  for(i in 1:ny){

    Y_c_1[i,,] = rbind(rep(-0.36211, T), rep(-0.83016, T)) + t(t((c(-1.38507, -2.65461)))) %*% t(f.nl((0:(T-1)))) +
      rbind(rep(u1[nx + i,1], T), rep(u2[nx + i,1], T)) +  t(t((c(u1[nx + i,2], u2[nx + i,2])))) %*% t(((0:(T-1)))) +
      t(t((c(effect_size_1, effect_size_2)))) %*% t(f.dev((0:(T-1)))) +
      t(as.matrix(rt2d(T, rho = 0.5, nu = 4)))/sqrt(3)
    # (t(mvrnorm(T, c(0,0), sigma)))
    Y_c_1[i,,] = Y_c_1[i,,] - rep.col(Y_c_1[i,,1], T)

    Y_c[i,,] = Y_c_1[i,,] + rbind(rep(0.3 * allele[nx + i], T), rep(0.2 * allele[nx + i], T)) +
      t(t(c(0.1 * allele[nx + i] , 0.5 * allele[nx + i]))) %*% (t((0 : (T-1)))) +
      0.3 * cov_time[nx + i,]
    Y_c[i,,] = Y_c[i,,] - rep.col(Y_c[i,,1], T)


  }


  return(list("Placebo" = X_c_1, "Treatment" = Y_c_1,
              "Placebo.Cov" = X_c, "Treatment.Cov" = Y_c,
              "allele" = allele, "cov_time" = cov_time))
}




