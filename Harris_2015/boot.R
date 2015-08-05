# Appendix 4: bootstrap estimates of model stability

library(rosalia)

set.seed(1)

n_boot = 50

# Load a simulated data set with 20 species at 200 sites.
# Multiply by one to avoid problems with FALSE/TRUE vs 0/1
d = readRDS("fakedata/20-200-1.rds")$observed * 1

# Fit a Markov network to the full data set
full_fit = rosalia(d, prior = make_logistic_prior(scale = 2), hessian = TRUE,
                   maxit = 200)


# Repeat the analysis on each bootstrap replicate
boots = replicate(
  n_boot,
  {
    boot_sample = d[sample.int(nrow(d), replace = TRUE), ]
    
    # Initializing the optimization at full_fit$opt$par should save 
    # time and shouldn't change the final outcome because the prior
    # and likelihood are both convex
    boot_fit = rosalia(
      boot_sample, 
      prior = make_logistic_prior(scale = 2), 
      parlist = relist(full_fit$opt$par),
      trace = 0
    )
    
    # Output the upper triangle of the fitted model
    boot_fit$beta[upper.tri(boot_fit$beta)]
  }
)

# Correlations among the bootstrapped estimates
# (upper triangle avoids the diagonal, which is always 1)
boot_cors = cor(boots)[upper.tri(cor(boots))]

# Mean squared correlation among replicates = 0.494
mean(boot_cors^2)
