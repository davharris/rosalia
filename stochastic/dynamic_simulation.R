set.seed(1)

library(progress)

n_spp = 20
K = n_spp / 3

maxit = 1E6
thin = 2E3

seed_rain = rep(.1, n_spp)
growth = rep(1, n_spp)

interaction_vector = rexp(choose(n_spp, 2), 1/2)

interaction_signs = ifelse(
  rbinom(choose(n_spp, 2), size = 1, prob = .75),
  -1, 
  1
)

signed_interaction_vector = interaction_vector * interaction_signs
competition_vector = interaction_vector * (interaction_signs == -1)
mutualism_vector = interaction_vector * (interaction_signs == +1)

mutualism_scaler = function(x){
  (1 + exp(-4 * x / K / n_spp)) / 2
}

competition = matrix(0, nrow = n_spp, ncol = n_spp)
competition[upper.tri(competition)] = competition_vector
competition = competition + t(competition)

mutualism = matrix(0, nrow = n_spp, ncol = n_spp)
mutualism[upper.tri(mutualism)] = mutualism_vector
mutualism = mutualism + t(mutualism)



pops = matrix(nrow = n_spp, ncol = maxit / thin)

population = rep(ceiling(K / n_spp / mean(competition)), n_spp)

pb = progress_bar$new(total = maxit / thin)

for(i in 1:maxit){
  local_growth = growth * population
  interactions = population * mutualism_scaler(population %*% mutualism) + 
    population %*% competition
  
  birthrates = seed_rain + local_growth
  deathrates = local_growth * interactions / K
  
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
    pb$tick()
  }
}

matplot(t(pops), type = "l", lty = 1)
sort(rowMeans(pops >0))

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
cor(fit$beta[upper.tri(fit$beta)], signed_interaction_vector, method = "spearman")

# Correlation performance
cor(
  cor(t(pops > 0))[upper.tri(cor(t(pops > 0)))], 
  signed_interaction_vector, 
  method = "spearman"
)

# Partial correlation performance
cor( 
  -solve(cov(t(pops > 0)))[upper.tri(cor(t(pops > 0)))], 
  signed_interaction_vector, 
  method = "spearman"
)

# BayesComm performance
cor( 
  colMeans(bc$trace$R), 
  signed_interaction_vector, 
  method = "spearman"
)

mles = fit$beta[upper.tri(fit$beta)]
ses = sqrt(diag(solve(fit$opt$hessian)))[-(1:n_spp)]
is_significant = ifelse(mles > 0, mles - 2 * ses > 0, mles + 2 * ses < 0)

plot(
  fit$beta[upper.tri(fit$beta)], 
  signed_interaction_vector,
  col = 1 + is_significant
)
abline(v = 0, h = 0)

#xy = matrix(numeric(0), nrow = 0, ncol = 2)
xy = rbind(xy, cbind(fit$beta[upper.tri(fit$beta)], signed_interaction_vector))
plot(xy, col = 1 + (0:(nrow(xy)-1) %/% choose(n_spp, 2)))
