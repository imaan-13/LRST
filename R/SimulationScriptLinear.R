#### ADAS
# Define parameters
N <- 300  # Number of patients
T <- 7   # Number of time points
sigma_e_1 <- 3.477  # Error standard deviation
sigma_u_11 <- 3.665  # Random effect Intercept standard deviation
sigma_u_12 <- 1.396 # Random effect Slope standard deviation
correlation_1 <- 0.2  # Correlation between random intercepts and slopes
effect_size_1 <-  2.21/6  # Effect size
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

##### DAD

sigma_e_2 <- 7.232  # Error standard deviation
sigma_u_21 <- 9.036  # Random effect Intercept standard deviation
sigma_u_22 <- 2.871 # Random effect Slope standard deviation
correlation_2 <- 0.08  # Correlation between random intercepts and slopes
effect_size_2 <- 5.38/6  # Effect size



#### Adding covariates


beta11 = 0.5
beta12 = 0.9

beta21 = 0.9
beta22 = 0.1

prob = 0.7
allele_Prob = c(prob, 2*sqrt(prob)*(1- sqrt(prob)) , (1 - sqrt(prob))^2)
times = seq(0, 1, length = 10)
n = 1000
allele = sample(c(0,1,2), N, prob = allele_Prob, replace = T)
alleleFull <- rep(allele, each = T)  # Treatment condition for each patient


#### Generate

sigma = matrix(c(sigma_e_1^2, 0.5*sigma_e_1*sigma_e_2, 0.5*sigma_e_1*sigma_e_2, sigma_e_2^2), 2, 2)
# Assign patients to placebo and treatment groups in a 2:3 ratio
num_placebo <- round(N * 2 / 5)
num_treatment <- N - num_placebo
treatment_assign <- c(rep("Placebo", num_placebo), rep("Treatment", num_treatment))

# Generate random data for patients
group <- rep(1:N, each = T)  # Each patient has T observations
time <- rep(1:T, times = N)   # Time variable
treatment <- rep(treatment_assign, each = T)  # Treatment condition for each patient

library(MASS)
# Generate random intercepts and slopes for each patient


lmm = numeric(1000)
lmm2 =numeric(1000)
lrst = numeric(1000)
nfars = numeric(1000)
nfars2 = numeric(1000)
lrst.discrete = numeric(1000)
nfars.discrete = numeric(1000)
lmm.cov = numeric(1000)
lmm.cov.2 = numeric(1000)
lrst.cov = numeric(1000)
nfars.cov = numeric(1000)

for(iter in 1:1000){

  u1 <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(sigma_u_11^2,
                                                  correlation_1 * sigma_u_11 * sigma_u_12,
                                                  correlation_1 * sigma_u_11 * sigma_u_12,
                                                  sigma_u_12^2), nrow = 2))
  u2 <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(sigma_u_21^2,
                                                  correlation_2 * sigma_u_21 * sigma_u_22,
                                                  correlation_2 * sigma_u_21 * sigma_u_22,
                                                  sigma_u_22^2), nrow = 2))


  # u_intercept_1 <- rep(u1[,1], each = T)
  # u_intercept_2 <- rep(u2[,1], each = T)
  #
  # u_slope_1 <- rep(u1[,2], each = T)
  # u_slope_2 <- rep(u2[,2], each = T)

  # # y <- cbind(rep(-0.36211, N*T), rep(-0.83016, N*T)) +  cbind(-1.38507 * time, -2.65461 * time) +
  # #   cbind(effect_size_1 * as.numeric(treatment == "Treatment") ,
  # #         effect_size_2 * as.numeric(treatment == "Treatment")) +
  # #   cbind((effect_size_1 * as.numeric(treatment == "Treatment") * time),
  # #         (effect_size_2 * as.numeric(treatment == "Treatment") * time))  +
  # #   c(u_intercept_1, u_intercept_2) + cbind(u_slope_1 * time, u_slope_2 * time) + mvrnorm(N * T, c(0,0), sigma)
  #
  # data <- data.frame(
  #   y = y,
  #   time = time,
  #   treatment = factor(treatment),
  #   patient = factor(group)
  # )
  #
  # data1 = data[,c(1,3,4,5)]
  # data2 = data[,c(2,3,4,5)]

  nx = num_placebo
  ny = num_treatment
  X_c = array(NA, dim = c(nx, 2, 7))
  Y_c = array(NA, dim = c(ny, 2, 7))
  #
  # d1 = data1 %>%
  #   pivot_wider(names_from = time, values_from = y.1)
  # d2 = data2 %>%
  #   pivot_wider(names_from = time, values_from = y.2)


  for(i in 1:nx){
    # X_c[i,,] = as.matrix((rbind(d1[i,], d2[i,]))[, -c(1,2)])
    # X_c[i,,] = X_c[i,,] - cbind(t(t(X_c[i,,1])), matrix(0, 2, T-1))

    X_c[i,,] = rbind(rep(-0.36211, T), rep(-0.83016, T)) + t(t((c(-1.38507, -2.65461)))) %*% (t((0:(T-1)))) +
      rbind(rep(u1[i,1], T), rep(u2[i,1], T)) +  t(t((c(u1[i,2], u2[i,2])))) %*% (t((0:(T-1)))) +
      t(mvrnorm(T, c(0,0), sigma))
    X_c[i,,] = X_c[i,,] - rep.col(X_c[i,,1], T)
  }


  for(i in 1:ny){
    Y_c[i,,] = rbind(rep(-0.36211, T), rep(-0.83016, T)) + t(t((c(-1.38507, -2.65461)))) %*% (t((0:(T-1)))) +
      rbind(rep(u1[i,1], T), rep(u2[i,1], T)) +  t(t((c(u1[i,2], u2[i,2])))) %*% (t((0:(T-1)))) +
      rbind(rep(effect_size_1, T), rep(effect_size_2, T)) + t(t((c(effect_size_1, effect_size_2)))) %*% (t((0:(T-1)))) +
      t(mvrnorm(T, c(0,0), sigma))
    Y_c[i,,] = Y_c[i,,] - rep.col(Y_c[i,,1], T)
  }

  lmm[iter] = max(lmmFun4(X_c, Y_c))
  lmm2[iter] = min(lmmFun4(X_c, Y_c))
  lrst[iter] = uniUstat(X_c, Y_c)
  nfars[iter] = min(npld(X_c, Y_c))
  nfars2[iter] = max(npld(X_c, Y_c))

  X_c_d = discretize(X_c)
  Y_c_d = discretize(Y_c)

  lrst.discrete[iter] = uniUstat(X_c_d, Y_c_d)
  nfars.discrete[iter] = npld(X_c_d, Y_c_d)

  #
  # y.cov = cbind(rep(-0.36211, N*T), rep(-0.83016, N*T)) +  cbind(-1.38507 * time, -2.65461 * time) +
  #   beta4 * alleleFull + beta5 * alleleFull * time +
  #   cbind(effect_size_1 * as.numeric(treatment == "Treatment") ,
  #         effect_size_2 * as.numeric(treatment == "Treatment")) +
  #   cbind((effect_size_1 * as.numeric(treatment == "Treatment") * time),
  #         (effect_size_2 * as.numeric(treatment == "Treatment") * time))  +
  #   c(u_intercept_1, u_intercept_2) + cbind(u_slope_1 * time, u_slope_2 * time) + mvrnorm(N * T, c(0,0), sigma)


  X_c_cov = array(NA, dim = c(nx, 2, 7))
  Y_c_cov = array(NA, dim = c(ny, 2, 7))
  #
  # d1 = data1 %>%
  #   pivot_wider(names_from = time, values_from = y.1)
  # d2 = data2 %>%
  #   pivot_wider(names_from = time, values_from = y.2)


  for(i in 1:nx){
    # X_c[i,,] = as.matrix((rbind(d1[i,], d2[i,]))[, -c(1,2)])
    # X_c[i,,] = X_c[i,,] - cbind(t(t(X_c[i,,1])), matrix(0, 2, T-1))

    X_c_cov[i,,] = rbind(rep(-0.36211, T), rep(-0.83016, T)) + t(t((c(-1.38507, -2.65461)))) %*% (t((0:(T-1)))) +
      rbind(rep(u1[i,1], T), rep(u2[i,1], T)) +  t(t((c(u1[i,2], u2[i,2])))) %*% (t((0:(T-1)))) +
      t(mvrnorm(T, c(0,0), sigma)) + rbind(rep(beta11 * allele[i], T) , rep(beta21 * allele[i], T)) +
      t(t((c(beta12 * allele[i], beta22 * allele[i])))) %*% (t((0:(T-1))))
    # X_c[i,,] = X_c[i,,] - cbind(t(t(X_c[i,,1])), matrix(0, 2, T-1))
  }


  for(i in 1:ny){
    Y_c_cov[i,,] = rbind(rep(-0.36211, T), rep(-0.83016, T)) + t(t((c(-1.38507, -2.65461)))) %*% (t((0:(T-1)))) +
      rbind(rep(u1[i,1], T), rep(u2[i,1], T)) +  t(t((c(u1[i,2], u2[i,2])))) %*% (t((0:(T-1)))) +
      rbind(rep(effect_size_1, T), rep(effect_size_2, T)) + t(t((c(effect_size_1, effect_size_2)))) %*% (t((0:(T-1)))) +
      t(mvrnorm(T, c(0,0), sigma))  + rbind(rep(beta11 * allele[i], T) , rep(beta21 * allele[i], T)) +
      t(t((c(beta12 * allele[i], beta22 * allele[i])))) %*% (t((0:(T-1))))
    # Y_c[i,,] = Y_c[i,,] - cbind(t(t(Y_c[i,,1])), matrix(0, 2, T-1))
  }
  A = numeric(2)

  X_c_n = array(NA, dim = c(nx, 2, 7))
  Y_c_n = array(NA, dim = c(ny, 2, 7))


  for(k in 1:2){
    out1 = cbind(1:(nx+ny),rbind(X_c_cov[,k,], Y_c_cov[,k,]), c(rep(0,nx), rep(1,ny)), allele)
    out1 = data.frame(out1)
    names(out1) = c("Subject", 0:(T-1),  "Treatment", "Allele")
    long1 = tidyr::pivot_longer(out1, -c(Subject, Treatment, Allele),
                                values_to = "Value", names_to = "Day")
    long1$Day = as.numeric(long1$Day)
    fit1 = lme4::lmer(Value ~  Allele + Day +  Allele*Day + Treatment + Treatment*Day + (1 + Day | Subject), data = long1)
    fit2 = lme4::lmer(Value ~  Allele + Day +  Allele*Day + (1 + Day | Subject), data = long1)

    long1$fitted = fitted.values(fit2)
    ld = long1[,c(1,4,6)]
    ld = ld %>%
      pivot_wider(names_from = Day, values_from = fitted)
    X_c_n[,k,] = as.matrix(ld[1:nx,-1])
    Y_c_n[,k,] = as.matrix(ld[(nx+1):(nx+ny),-1])
    A[k] = anova(fit1, fit2)$`Pr(>Chisq)`[2]
  }

  lmm.cov[iter] = max(A)
  lmm.cov.2[iter] = min(A)
  lrst.cov[iter] = uniUstat(X_c_n, Y_c_n)
  nfars.cov[iter] = npld(X_c_n, Y_c_n)

  print(iter)
}

length(which(lmm < 0.05))/1000
length(which(lmm2 < 0.05))/1000
length(which(lrst< 0.05))/1000
length(which(nfars< 0.05))/1000
length(which(lrst.discrete < 0.05))/1000
length(which(nfars.discrete< 0.05))/1000
length(which(lrst.cov< 0.05))/1000
length(which(nfars.cov< 0.05))/1000
length(which(lmm.cov< 0.05))/1000
length(which(lmm.cov.2< 0.05))/1000
