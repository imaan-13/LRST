# Introduction

The `LRST` package implements the Longitudinal Rank-Sum Test for
analyzing longitudinal clinical trial data.

# Installation

You can install the package using:

    knitr::opts_chunk$set(echo = TRUE)
    library(LRST)
    library(MASS)

The `lrst.2arm` function perform LRST on placebo and treatment data,
while the `lrst.MultiArm` perform LRST on multi-arm clinical trials.
Before perfomring LRST let us see what the data should look like.

# Simulated Data Generation

The below code generates the mean and standard deviation curves for the
plaecbo group. The numbers are taken from the BAPI 302 trial.

    xm = -t(cbind(c(0, 0.738, 1.313, 3.109, 4.525, 5.864, 7.338),
                  c(0, 0.668, 3.975, 5.931, 8.259, 11.989, 13.976)))

    xs = t(cbind(c(0, 4.79, 5.43, 6.54, 7.37, 8.15, 9.11),
                 c(0, 10.27, 12.85, 14.95, 15.35, 16.87, 18.19)))

Now let’s provide the dimentions for the simulation.

    N = 500 # Total Number of Patients
    nx = 2/5*N # Niumber of patients in Placebo
    ny = 3/5*N # Number of patients in Treatment
    K = 2 # Number of Outcomes
    T = 6 # Number of visits other than baseline

Let’s define the temporal covariance structure for the two
outcomes.Suppose the correlation between the outcomes is 0.5.

    sigma_u <- c(3.477, 7.232) # Random Effect Standard Deviation for Outcomes 
    corr <- matrix(c(1,0.5, 0.5, 1),2,2)  # Correlation between the outcomes

    sigma = diag(sigma_u) %*% corr %*% diag(sigma_u)

<!-- Suppose the correlation between the outcomes is 0.5. -->
<!-- ```{r Correlation} -->
<!-- corr = matrix(c(1,0.5,0.5,1),2,2) -->
<!-- ``` -->

Let’s get the treatment mean curve, the standard deviations are assumed
same as the placebo group.

    alpha=1
    beta=1
    effect_size_1 = 2.21*alpha/6
    effect_size_2 = 5.38*beta/6

    ym = xm + t(cbind(c(0,(rep(effect_size_1,6))), c(0,(rep(effect_size_2,6)))))

Finally Let’s simulate the data.

    X_placebo = array(NA, dim = c(nx, 2, 7))
    X_dose_1 = array(NA, dim = c(ny, 2, 7))

    for(i in 1:nx){
      X_placebo[i,,] = gen.fun.dep(xm, xs, corr, sigma)
    }

    for(i in 1:ny){
      X_dose_1[i,,] = gen.fun.dep(ym, xs, corr, sigma)
    }

    X_c = X_placebo
    Y_c = X_dose_1
    dim(X_c)
    #> [1] 200   2   7
    dim(Y_c)
    #> [1] 300   2   7

Let’s perform LRST on the generated data.

    lrst.2arm(X_c, Y_c)
    #> $T.stat
    #> [1] 3.150061
    #> 
    #> $T.sd
    #> [1] 4.598027
    #> 
    #> $p.value
    #> [1] 0.07091118

Let’s get the estimated power.

    estimated_Power(0.05, X_c, Y_c)
    #> [1] 0.5138214

Let’s conduct the Multi-Arm LRST.

    X_dose_2 = array(NA, dim = c(ny, 2, 7))
    for(i in 1:ny){
      X_dose_2[i,,] = gen.fun.dep(ym, xs, corr, sigma)
    }
    X_c = X_placebo
    Y = vector(length = 2, mode = "list")
    Y[[1]] = X_dose_1
    Y[[2]] = X_dose_2
    lrst.MultiArm(X_c, Y)
    #> $p
    #> [1] 0.11167
    #> 
    #> $v
    #> [1] 1
    #> 
    #> $T
    #> [1] 1.564431
