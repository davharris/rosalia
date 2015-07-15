set.seed(1)

n_spp = 15
K = 15

maxit = 1E7
thin = 2E4

seed_rain = rep(.1, n_spp)
growth = rep(1, n_spp)

competition_vector = ifelse(
  rbinom(choose(n_spp, 2), size = 1, prob = .5),
  rexp(choose(n_spp, 2), 2),
  rexp(choose(n_spp, 2), 1/2)
)

competition = matrix(0, nrow = n_spp, ncol = n_spp)
competition[upper.tri(competition)] = competition_vector
competition = competition + t(competition)
diag(competition) = 1

pops = matrix(nrow = n_spp, ncol = maxit / thin)

population = rep(ceiling(K / n_spp / mean(competition)), n_spp)

for(i in 1:maxit){
  local_growth = growth * population
  birthrates = seed_rain + local_growth
  deathrates = local_growth * (population %*% competition) / K
  
  rates = c(birthrates, deathrates)
  
  event = sample.int(length(rates), 1, prob = rates)
  if(event <= length(birthrates)){
    # Event from first batch, i.e. birth.  Increase the population by one
    index = event
    population[index] = population[index] + 1
  }else{
    # Event from the second batch, i.e. death. If population is nonzero,
    # decrease it by one for the corresponding species
    index = event - n_spp
    if(population[index] > 0){
      population[index] = population[index] - 1
    }
  }
  if(i%%thin == 0){
    pops[ , i %/% thin] = population
  }
}

matplot(t(pops), type = "l", lty = 1)


library(rosalia)
fit = rosalia(
  t(pops > 0) * 1, 
  hessian = TRUE, 
  maxit = 1E3,
  prior = make_logistic_prior(location = 0, scale = 2)
)



library(BayesComm)
bc = BC(Y = t(pops > 0), model = "community", its = 1000, thin = 10)


# Markov network performance
cor(fit$beta[upper.tri(fit$beta)], -competition_vector, method = "spearman")

# Correlation performance
cor(
  cor(t(pops > 0))[upper.tri(cor(t(pops > 0)))], 
  -competition_vector, method = "spearman"
)

# Partial correlation performance
cor( 
  -solve(cor(t(pops > 0)))[upper.tri(cor(t(pops > 0)))], 
  -competition_vector, 
  method = "spearman"
)

# BayesComm performance
cor( 
  colMeans(bc$trace$R), 
  -competition_vector, 
  method = "spearman"
)

mles = fit$beta[upper.tri(fit$beta)]
ses = sqrt(diag(solve(fit$opt$hessian)))[-(1:n_spp)]
is_significant = ifelse(mles > 0, mles - 2 * ses > 0, mles + 2 * ses < 0)

plot(
  fit$beta[upper.tri(fit$beta)], 
  -competition_vector,
  col = 1 + is_significant
)
abline(v = 0, h = 0)

#xy = matrix(numeric(0), nrow = 0, ncol = 2)
xy = rbind(xy, cbind(fit$beta[upper.tri(fit$beta)], -competition_vector))
plot(xy, col = 1 + (0:(nrow(xy)-1) %/% choose(n_spp, 2)))
