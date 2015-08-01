library(rosalia)

n_boot = 50

d = readRDS("fakedata/20-200-1.rds")$observed * 1

full_fit = rosalia(d, prior = make_logistic_prior(scale = 2), hessian = TRUE)



boots = replicate(
  n_boot,
  {
    boot_sample = d[sample.int(nrow(d), replace = TRUE), ]
    boot_fit = rosalia(
      boot_sample, 
      prior = make_logistic_prior(scale = 2), 
      parlist = full_fit$opt$par
    )
    boot_fit$beta[upper.tri(boot_fit$beta)]
  }
)
