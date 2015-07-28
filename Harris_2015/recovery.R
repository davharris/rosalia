library(dplyr)
library(corpcor)
library(rosalia)
library(arm)
library(BayesComm)

set.seed(1)


# Load in the results from the `pairs` program
pairs_txt = readLines("fakedata/Pairs.txt")

# Find areas of the data file that correspond
# to species pairs' results
beginnings = grep("Sp1", pairs_txt) + 1
ends = c(
  grep("^[^ ]", pairs_txt)[-1],
  length(pairs_txt)
) - 1
# The above code breaks on the very last line of the file
ends[length(ends)] = ends[length(ends)] + 1


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
  
  splitname = strsplit(filename, "/|-|\\.")[[1]]
  n_sites = as.integer(splitname[[3]])
  rep = as.integer(splitname[[4]])
  
  
  ####### Partial correlations
  p_corr = pcor.shrink(x)
  
  ####### Correlations
  corr = cor(x)
  
  ####### GLM
  coef_matrix = matrix(0, ncol(x), ncol(x))
  for(i in 1:ncol(x)){
    if(var(x[,i]) > 0){
      coefs = coef(bayesglm(x[,i] ~ x[ , -i], family = binomial))[-1]
      coef_matrix[i, -i] = coefs
    }
  }
  coef_matrix = (coef_matrix + t(coef_matrix)) / 2
  
  ####### Markov network
  rosie = rosalia(x, maxit = 200, trace = 0, prior = make_logistic_prior(scale = 2))
  
  ####### BayesComm and partial BayesComm
  bc = BC(Y = x, model = "community", its = 1000)
  
  `partial BayesComm` = 0
  
  for(i in 1:nrow(bc$trace$R)){
    Sigma = matrix(0, nrow = ncol(x), ncol = ncol(x))
    Sigma[upper.tri(Sigma)] <- bc$trace$R[i, ]  # Fill in upper triangle
    Sigma <- Sigma + t(Sigma)                   # Fill in lower triangle
    diag(Sigma) <- 1  # Diagonal equals 1 in multivariate probit model
    
    `partial BayesComm` = `partial BayesComm` + cor2pcor(Sigma) / nrow(bc$trace$R)
  }
  
  ####### Pairs
  # Find the line where the current data set is mentioned in
  # pairs.txt
  filename_line = grep(
    paste0(
      gsub("fakedata/(.*)\\.rds", "\\1", filename),
      "-"
    ),
    pairs_txt
  )
  
  # Which chunk of the data file corresponds to this file?
  chunk = min(which(beginnings > filename_line))
  
  # Split the chunk on whitespace.
  splitted = strsplit(pairs_txt[beginnings[chunk]:ends[chunk]], " +")
  
  # Pull out the species numbers and their Z-scores
  pairs_results = lapply(
    splitted,
    function(x){
      spp = sort(as.integer(x[3:4]))
      data.frame(
        sp1 = spp[1], 
        sp2 = spp[2], 
        z = as.numeric(x[14])
      )
    }
  ) %>% 
    bind_rows %>%
    mutate(spp = paste(sp1, sp2, sep = "-"))
  
  # Re-order the pairs_results to match the other methods
  m = matrix(NA, ncol(x), ncol(x))
  new_order = match(
    paste(row(m)[upper.tri(m)], col(m)[upper.tri(m)], sep = "-"),
    pairs_results$spp
  )
  ordered_pairs_results = pairs_results[new_order, ]
  
  ####### Output
  data.frame(
    truth = truth, 
    n_sites = n_sites,
    rep = rep,
    sp1 = ordered_pairs_results$sp1,
    sp2 = ordered_pairs_results$sp2,
    `partial correlation` = p_corr[upper.tri(p_corr)],
    correlation = corr[upper.tri(corr)],
    `Markov network` = rosie$beta[upper.tri(rosie$beta)],
    GLM = coef_matrix[upper.tri(coef_matrix)],
    `BayesComm` = colMeans(bc$trace$R),
    `partial BayesComm` = `partial BayesComm`[upper.tri(`partial BayesComm`)],
    Pairs = ordered_pairs_results$z
  )
}



files = dir("fakedata", pattern = "\\.rds$", full.names = TRUE)

z = lapply(files, f) %>% bind_rows %>% as.data.frame
colnames(z) = gsub("\\.", " ", colnames(z))

library(ggplot2)
ggplot(z[z$n_sites >0, ], aes(x = `Markov network`, y = truth)) + 
  stat_binhex(bins = 100) + 
  xlab("Estimated coefficient value") + 
  ylab("\"True\" coefficient value") + 
  stat_hline(yintercept = 0, size = 1/8) + 
  stat_vline(xintercept = 0, size = 1/8) + 
  scale_fill_gradient(low = "#F0F0F0", high = "darkblue", trans = "identity") +
  stat_smooth(method = "lm", se = FALSE, size = 1) +
  theme_bw()

ggplot(z[z$n_sites >0, ], aes(x = -`Pairs`, y = truth)) + 
  stat_binhex(bins = 100) + 
  xlab("Estimated coefficient value") + 
  ylab("\"True\" coefficient value") + 
  stat_hline(yintercept = 0, size = 1/8) + 
  stat_vline(xintercept = 0, size = 1/8) + 
  scale_fill_gradient(low = "#F0F0F0", high = "darkblue", trans = "identity") +
  stat_smooth(method = "lm", se = FALSE, size = 1) +
  theme_bw()




resids = sapply(
  colnames(z)[-(1:5)],
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
results = results[order(colMeans(results), decreasing = TRUE)]

library(dichromat)
library(RColorBrewer)
colors = brewer.pal(8, "Dark2")[c(1, 2, 3, 4, 8, 5, 7)]

pdf("manuscript-materials/figures/performance.pdf", height = 4, width = 5)
par(mar = c(5, 4, 1, 0) + .1)
matplot(
  sizes, 
  results, 
  type = "o", 
  ylab = "", 
  xlim = c(25 * .9, 1600 * 8),
  ylim = c(0, .5 + 1E-9),
  pch = c(15, 1, 17, 0, 16, 2, 3), 
  lty = 1,
  log = "x",
  xlab = "Number of sites (log scale)",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i",
  col = colors,
  lwd = 1.25,
  cex = 0.9
)
axis(1, c(1, 25, 200, 1600))
axis(2, seq(0, 1, .1), las = 1)
mtext(expression(R^2), side = 2, line = 3, las = 1)
heights = results[nrow(results), ]

heights$Pairs = heights$Pairs - .01
heights$correlation = heights$correlation + .01

heights$`partial correlation` = heights$`partial correlation` + .01
heights$`partial BayesComm` = heights$`partial BayesComm` - .01

heights$GLM = heights$GLM - .01


text(1625, heights, colnames(results), pos = 4, cex = 0.75, col = colors)
dev.off()


total_ss = mean(resid(lm(truth ~ 0, data = z))^2)
round(100 * sort(1 - colMeans(resids^2) / total_ss), 1)


