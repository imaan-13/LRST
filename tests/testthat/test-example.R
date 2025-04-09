library(LRST)
library(MASS)
library(nparLD)


test_that("lrst.2arm returns correct structure on simulated data", {
  xm = -t(cbind(c(0, 0.738, 1.313, 3.109, 4.525, 5.864, 7.338),
                c(0, 0.668, 3.975, 5.931, 8.259, 11.989, 13.976)))

  xs = t(cbind(c(0, 4.79, 5.43, 6.54, 7.37, 8.15, 9.11),
               c(0, 10.27, 12.85, 14.95, 15.35, 16.87, 18.19)))
  alpha=1
  beta=-1

  N = 300
  nx = 2/5*N
  ny = 3/5*N
  K = 2
  T = 6

  sigma_u <- c(3.477, 7.232) # Random Effect Standard Deviation for Outcomes
  corr <- matrix(c(1,0.5, 0.5, 1),2,2)  # Correlation between the outcomes

  sigma <- diag(sigma_u) %*% corr %*% diag(sigma_u)

  X_placebo <- array(NA, dim = c(nx, 2, 7))
  X_dose_1 <- array(NA, dim = c(ny, 2, 7))

  effect_size_1 <- 2.21*alpha/6
  effect_size_2 <- 5.38*beta/6

  ym <- xm + t(cbind(c(0,(rep(effect_size_1,6))), c(0,(rep(effect_size_2,6)))))



  for(i in 1:nx){
    X_placebo[i,,] <- gen.fun.dep(xm, xs, corr, sigma)
  }

  for(i in 1:ny){
    X_dose_1[i,,] <- gen.fun.dep(ym, xs, corr, sigma)
  }

  X_c <- X_placebo
  Y_c <- X_dose_1
  result <- lrst.2arm(X_c, Y_c)
  # lmmFun4(X_c, Y_c)
  power <- estimated_Power(0.05, X_placebo, X_dose_1)

  expect_type(result, "list")
  expect_named(result, c("T.stat", "T.sd", "p.value"))

  expect_true(is.numeric(result$T.stat))
  expect_true(is.numeric(result$T.sd))
  expect_true(is.numeric(result$p.value))

  expect_length(result$p.value, 1)
  expect_true(result$p.value >= 0 && result$p.value <= 1)

  expect_length(power, 1)
  expect_true(power >= 0 && result$p.value <= 1)



})




