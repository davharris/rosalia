library(dplyr)
library(corpcor)
library(rosalia)
library(arm)
library(BayesComm)

set.seed(1)

f = function(filename){
  # Multiplying by one is necessary to prevent silly errors
  # regarding TRUE/FALSE versus 1/0
  raw_obs = readRDS(filename)[["observed"]] * 1
  
  # Identify species that are never present (or never absent) so they 
  # can be dropped
  species_is_variable = diag(var(raw_obs)) > 0
  pair_is_variable = tcrossprod(species_is_variable) > 0
  
  x = raw_obs[ , species_is_variable]
  truth = readRDS(filename)[["truth"]][pair_is_variable[upper.tri(pair_is_variable)]]
  
  p_corr = pcor.shrink(x)
  
  splitname = strsplit(filename, "/|-|\\.")[[1]]
  n_sites = as.integer(splitname[[3]])
  rep = as.integer(splitname[[4]])
  
  corr = cor(x)
  
  coef_matrix = matrix(0, ncol(x), ncol(x))
  for(i in 1:ncol(x)){
    if(var(x[,i]) > 0){
      coefs = coef(bayesglm(x[,i] ~ x[ , -i], family = binomial))[-1]
      coef_matrix[i, -i] = coefs
    }
  }
  coef_matrix = (coef_matrix + t(coef_matrix)) / 2
  
  rosie = rosalia(x, maxit = 200, trace = 0)
  
  bc = BC(Y = x, model = "community", its = 1000)
  
  data.frame(
    truth = truth, 
    n_sites = n_sites,
    rep = rep,
    `partial correlation` = p_corr[upper.tri(p_corr)],
    correlation = corr[upper.tri(corr)],
    `Markov network` = rosie$beta[upper.tri(rosie$beta)],
    GLM = coef_matrix[upper.tri(coef_matrix)],
    `BayesComm` = colMeans(bc$trace$R)
  )
}



files = dir("fakedata", pattern = "\\.rds$", full.names = TRUE)

z = lapply(files, f) %>% bind_rows %>% as.data.frame
colnames(z) = gsub("\\.", " ", colnames(z))

library(ggplot2)
ggplot(z[z$n_sites >0, ], aes(x = `Markov network`, y = truth)) + 
  stat_binhex(bins = 50) + 
  xlab("Estimated coefficient value") + 
  ylab("\"True\" coefficient value") + 
  stat_hline(yintercept = 0, size = 1/8) + 
  stat_vline(xintercept = 0, size = 1/8) + 
  scale_fill_gradient(low = "#F0F0F0", high = "darkblue", trans = "identity") +
  stat_smooth(method = "lm", se = FALSE, size = 1) +
  theme_bw()

resids = sapply(
  colnames(z)[-(1:3)],
  function(i){resid(lm(z$truth ~ z[,i] + 0))}
)

sizes = sort(unique(z$n_sites))
results = as.data.frame(
  t(
    sapply(
      sizes,
      function(n){
        total_ss = mean(resid(lm(truth ~ 0, data = z))[z$n_sites == n]^2)
        1 - colMeans(resids[z$n_sites == n , ]^2) / total_ss
      }
    )
  )
)


library(dichromat)
colors = colorschemes$Categorical.12[c(2, 6, 8, 10, 12)]

matplot(
  sizes, 
  results, 
  type = "o", 
  ylab = "", 
  xlim = c(25 * .9, 1600 * 5),
  ylim = c(0, .52),
  pch = 19, 
  lty = 1,
  log = "x",
  xlab = "Number of sites (log scale)",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i",
  col = colors,
  lwd = 3
)
axis(1, c(1, 25, 200, 1600))
axis(2, seq(0, 1, .1), las = 1)
mtext(expression(R^2), side = 2, line = 3)
text(1625, results[nrow(results), ], colnames(results), pos = 4)
