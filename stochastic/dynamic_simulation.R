n_spp = 15
K = 20

maxit = 1E7
thin = 4E3

seed_rain = rep(.1, n_spp)
growth = rep(1, n_spp)

competition_vector = ifelse(
  rbinom(choose(n_spp, 2), size = 1, prob = .5),
  rexp(choose(n_spp, 2), 5),
  rexp(choose(n_spp, 2), 1/5)
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


library(corpcor)
pcor = pcor.shrink(t(pops))
class(pcor) = "matrix"
lattice::levelplot(pcor^2)

plot(pcor[upper.tri(pcor)], -competition_vector)

library(rosalia)
fit = rosalia(
  t(pops > 0) * 1, 
  hessian = TRUE, 
  maxit = 1E3
)

mles = fit$beta[upper.tri(fit$beta)]
ses = sqrt(diag(solve(fit$opt$hessian)))[-(1:n_spp)]
is_significant = ifelse(mles > 0, mles - 2 * ses > 0, mles + 2 * ses < 0)


cor(fit$beta[upper.tri(fit$beta)], -competition_vector, method = "spearman")
plot(
  fit$beta[upper.tri(fit$beta)], 
  -competition_vector,
  col = 1 + is_significant
)
abline(v = 0, h = 0)
cor(
  cor(t(pops > 0))[upper.tri(cor(t(pops > 0)))], 
  -competition_vector, method = "spearman"
)


