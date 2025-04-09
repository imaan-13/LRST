# LRST: Longitudinal Rank-Sum Test for Multivariate Longitudinal Data

The **LRST** R package implements the Longitudinal Rank-Sum Test, a
nonparametric, rank-based omnibus test statistic designed for
comprehensive assessment of treatment efficacy across multiple
longitudinal primary endpoints in clinical trials. It effectively
controls Type I error while enhancing statistical power, offering
flexibility against various data distributions encountered in clinical
research.

## Key Features

-   **Comprehensive Assessment**: Evaluates treatment effects across
    multiple endpoints and time points without requiring multiplicity
    adjustments.
-   **Enhanced Statistical Power**: Effectively controls Type I error,
    reducing the need for larger sample sizes.
-   **Flexibility**: Robust against various data distributions, making
    it suitable for complex clinical trial data.
-   **Efficient Data Utilization**: Maximizes the use of longitudinal
    data, capturing dynamic changes across multiple endpoints.

# Installation

You can install the package using:

    knitr::opts_chunk$set(echo = TRUE)
    # install.packages("devtools")
    devtools::install_github("djghosh1123/LRST")
    #> Using GitHub PAT from the git credential store.
    #> Skipping install of 'LRST' from a github remote, the SHA1 (406a0e07) has not changed since last install.
    #>   Use `force = TRUE` to force installation
    library(LRST)

# Example

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
    #> [1] 1.080766
    #> 
    #> $T.sd
    #> [1] 4.800963
    #> 
    #> $p.value
    #> [1] 0.3109178

Let’s get the estimated power.

    estimated_Power(0.05, X_c, Y_c)
    #> [1] 0.1274444

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
    #> [1] 0.07858
    #> 
    #> $v
    #> [1] 2
    #> 
    #> $T
    #> [1] 1.730265

# Citation

If you use this paper, please cite:

Xu, X., Ghosh, D., Luo, S., & CPP Integrated Parkinson’s Database.
(2025). A novel longitudinal rank-sum test for multiple primary
endpoints in clinical trials: Applications to neurodegenerative
disorders. Statistics in Biopharmaceutical Research, (just-accepted),
1-17.
