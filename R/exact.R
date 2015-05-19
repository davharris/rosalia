#' Make a logistic prior with specified scale
#' 
#' @param location Numeric vector; see \code{\link{Logistic}}
#' @param scale Numeric vector; see \code{\link{Logistic}}
#' @param ...   Currently not used
#' @return A `prior` object, which is a list of two functions:
#' \item{log_d}{A function for calculating the log of the prior density for 
#' each element of a vector. Calculated with \code{\link{dlogis}} with \code{log = TRUE}}
#' \item{log_grad}{A function for calculating the gradient of the log-density for each element 
#' of a vector.}
#' @examples 
#' p = make_logistic_prior(location = 0, scale = 2)
#' curve(p$log_d(x), from = -5, to = 5)
#' curve(p$log_grad(x), from = -5, to = 5)
#' @export
make_logistic_prior = function(location = 0, scale, ...){
  structure(
    list(
      log_d = function(x){dlogis(x, location = location, scale = scale, log = TRUE)},
      log_grad = function(x){
        -tanh((x - location)/scale/2) / scale
      }
    ),
    class = "prior"
  )
}


#' Make a flat prior for maximum likelihood estimation
#' 
#' @param ...   Currently not used
#' @return A `prior` object, which is a list of two functions:
#' \item{log_d}{A function for calculating the log of the prior density for each element of a vector. 
#' With a uniform prior, the log-density is always zero.}
#' \item{log_grad}{A function for calculating the gradient of the log-density for each element of 
#' a vector. For the uniform prior, this is always 0. With a uniform prior, this gradient is always zero.}
#' @examples 
#' p = make_flat_prior()
#' curve(p$log_d(x), from = -5, to = 5)
#' curve(p$log_grad(x), from = -5, to = 5)
#' @export
make_flat_prior = function(...){
  structure(
    list(
      log_d = function(x){rep(0, length(x))},
      log_grad = function(x){rep(0, length(x))}
    ),
    class = "prior"
  )
}

# Add up the likelihoods without losing precision during exponentiation of log-likelihoods
logSumExp = function(x){
  biggest = max(x)
  log(sum(exp(x - biggest))) + biggest
}


# Generate all possible binary vectors of length n_spp
generate_possibilities = function(n_spp){
  possibilities = expand.grid(replicate(n_spp, c(0L, 1L), simplify = FALSE))
  as.matrix(possibilities[ , n_spp:1])
}

# Where do the observed vectors occur in the possibilities matrix?
find_rows = function(x){
  # string-to-integer
  # +1 because R indexes from 1 not 0
  strtoi(apply(x, 1, paste0, collapse = ""), 2) + 1
}


# Negative log-likelihood
#' @importFrom assertthat assert_that
nll = function(par, rows, possible_cooc, ...){    
  assert_that(!missing(rows))
  
  # Find energy for all possible co-occurrence patterns
  E = -c(par %*% possible_cooc)
  
  # Log of the partition function
  logZ = logSumExp(-E)
  
  # positive log-likelihood is sum(-energy - logZ)
  # This just distributes the minus sign
  sum(E[rows] + logZ)
}

# Negative log posterior
nlp = function(par, rows, possible_cooc, prior, ...){    
  nll(par, rows, possible_cooc, ...) - sum(prior$log_d(par))
}

# Gradient of the negative log-likelihood
nll_grad = function(par, possible_cooc, observed_cooc, n_sites, ...){
  
  # Find energy for all possible co-occurrence patterns
  E = -c(par %*% possible_cooc)
  
  # Log of the partition function
  logZ = logSumExp(-E)
  
  p = exp(-(E + logZ))
  
  expected_cooc = c(possible_cooc %*% p)
  
  n_sites * (expected_cooc - observed_cooc)
}

# Negative log posterior gradient
nlp_grad = function(
  par, 
  possible_cooc, 
  observed_cooc, 
  n_sites, 
  prior, 
  ...
){
  nll_grad(par, possible_cooc, observed_cooc, n_sites, ...) - prior$log_grad(par)
}

# Find the number of observed co-occurrences (and occurrences)
find_observed_cooc = function(x){
  cp = crossprod(x) / nrow(x)
  
  c(cp[upper.tri(cp)], diag(cp))
}

#' Fit a Markov network to binary data
#' 
#' @param x        A binary matrix.
#' @param prior    An object of class \code{prior}. By default, the prior is flat for maximum likelihood estimation
#' @param maxit,trace,hessian,...    Arguments passed to \code{\link{optim}}
#' @param parlist  user-specified starting values (optional).
#' @export
rosalia = function(
  x, 
  prior = make_flat_prior(), 
  maxit = 100, 
  trace = 1, 
  hessian = FALSE,
  parlist,
  ...
){
  n_spp = ncol(x)
  n_sites = nrow(x)
  
  if (n_spp > 25) {
    message(
      "Note: exact computation with more than 25 species is time-consuming and memory-intensive"
    )
  }
  
  possibilities = generate_possibilities(n_spp = n_spp)
  
  possible_cooc = sapply(
    1:2^n_spp,
    function(i){
      tcp = tcrossprod(possibilities[i, ])
      c(tcp[upper.tri(tcp)], diag(tcp))
    }
  )
  
  rm(possibilities) # no longer needed and possibly large drag on memory
  
  if (missing(parlist)) {
    parlist = as.relistable(list(
      upper = rep(0, choose(n_spp, 2)),
      diagonal = qlogis((colSums(x) + 1) / (nrow(x) + 2))
    ))
  }
  
  optim(
    unlist(parlist), 
    fn = nlp, 
    gr = nlp_grad, 
    method = "BFGS", 
    hessian = hessian,
    possible_cooc = possible_cooc,
    rows = find_rows(x),
    observed_cooc = find_observed_cooc(x),
    control = list(trace = trace, maxit = maxit, REPORT = 1, ...),
    n_sites = n_sites,
    prior = prior
  )
}

reform = function(par){
  parlist = relist(par)
  
  dims = length(parlist$diagonal)
  
  out = matrix(0, nrow = dims, ncol = dims)
  
  out[upper.tri(out)] = parlist$upper
  out = out + t(out)
  diag(out) = parlist$diagonal
  
  out
}
