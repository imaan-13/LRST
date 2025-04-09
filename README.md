LRST_tutorial
================

# Introduction

The `LRST` package implements the Longitudinal Rank-Sum Test for
analyzing longitudinal clinical trial data.

# Installation

You can install the package using:

``` r
knitr::opts_chunk$set(echo = TRUE)
library(LRST)
library(MASS)
```

The `uniUstat` function perform LRST on placebo and treatment data,
while the `multiGenUStat` perform LRST on multi-arm clinical trials.
Before perfomring LRST let us see what the data should look like.

# Simulated Data Generation

The below code generates the mean and standard deviation curves for the
plaecbo group. The numbers are taken from the BAPI 302 trial.

``` r
xm = -t(cbind(c(0, 0.738, 1.313, 3.109, 4.525, 5.864, 7.338),
              c(0, 0.668, 3.975, 5.931, 8.259, 11.989, 13.976)))

xs = t(cbind(c(0, 4.79, 5.43, 6.54, 7.37, 8.15, 9.11),
             c(0, 10.27, 12.85, 14.95, 15.35, 16.87, 18.19)))
```

Now let’s provide the dimentions for the simulation.

``` r
N = 300 # Total Number of Patients
nx = 2/5*N # Niumber of patients in Placebo
ny = 3/5*N # Number of patients in Treatment
K = 2 # Number of Outcomes
T = 6 # Number of visits other than baseline
```

Let’s define the temporal covariance structure for the two outcomes.

``` r
sigma_e_1 <- 3.477  # Error standard deviation
sigma_u_11 <- 3.665  # Random effect Intercept standard deviation
sigma_u_12 <- 1.396 # Random effect Slope standard deviation
correlation_1 <- 0.2  # Correlation between random intercepts and slopes
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

sigma_e_2 <- 7.232  # Error standard deviation
sigma_u_21 <- 9.036  # Random effect Intercept standard deviation
sigma_u_22 <- 2.871 # Random effect Slope standard deviation
correlation_2 <- 0.08  # Correlation between random intercepts and slopes

sigma = matrix(c(sigma_e_1^2, 0.5*sigma_e_1*sigma_e_2, 0.5*sigma_e_1*sigma_e_2, sigma_e_2^2), 2, 2)
Sigma1 = sigma
Sigma2 = sigma
```

Suppose the correlation between the outcomes is 0.5.

``` r
corr = 0.5
```

Let’s get the treatment mean curve, the standard deviations are assumed
same as the placebo group.

``` r
alpha=0.5
beta=0.5
effect_size_1 = 2.21*alpha/6
effect_size_2 = 5.38*beta/6

ym = xm + t(cbind(c(0,(rep(effect_size_1,6))), c(0,(rep(effect_size_2,6)))))
```

Finally Let’s simulate the data.

``` r
X_placebo = array(NA, dim = c(nx, 2, 7))
X_dose_1 = array(NA, dim = c(ny, 2, 7))

for(i in 1:nx){
  X_placebo[i,,] = gen.fun.dep(xm, xs, corr, Sigma1, Sigma2)
}

for(i in 1:ny){
  X_dose_1[i,,] = gen.fun.dep(ym, xs, corr, Sigma1, Sigma2)
}

X_c = X_placebo
Y_c = X_dose_1
dim(X_c)
#> [1] 120   2   7
dim(Y_c)
#> [1] 180   2   7
```

Let’s perform LRST on the generated data.

``` r
uniUstat(X_c, Y_c)
#> $T.stat
#> [1] 25.62192
#> 
#> $T.sd
#> [1] 15.67714
#> 
#> $p.value
#> [1] 4.864514e-11
```

Let’s get the estimated power.

``` r
estimated_Power(0.05, X_c, Y_c)
#> [1] 0.9999995
```
