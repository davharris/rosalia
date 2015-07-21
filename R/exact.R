#' Make a logistic prior with specified scale
#' 
#' @param location Numeric vector; see \code{\link[stats]{Logistic}}
#' @param scale Numeric vector; see \code{\link[stats]{Logistic}}
#' @param ...   Currently not used
#' @return A \code{prior} object, which is a list of two functions:
#' \item{log_d}{A function for calculating the log of the prior density for 
#' each element of a vector. Calculated with \code{\link[stats]{dlogis}} using \code{log = TRUE}}
#' \item{log_grad}{A function for calculating the gradient of the log-density for each element 
#' of a vector.}
#' @seealso \code{\link{make_flat_prior}}
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

#' @export
make_cauchy_prior = function(location = 0, scale, ...){
  structure(
    list(
      log_d = function(x){dcauchy(x, location = location, scale = scale, log = TRUE)},
      log_grad = function(x){
        (-2 * x + 2 * location)/(scale^2 + (x - location)^2)
      }
    ),
    class = "prior"
  )
}

#' @export
make_t_prior = function(df){
  structure(
    list(
      log_d = function(x){dt(x, df, log = TRUE)},
      log_grad = function(x){
        -(x * (1 + df)) / (x^2 + df)
      }
    ),
    class = "prior"
  )
}



#' Make a flat prior for maximum likelihood estimation
#' 
#' @param ...   Currently not used
#' @return A \code{prior} object, which is a list of two functions:
#' \item{log_d}{A function for calculating the log of the prior density for each element of a vector. 
#' With a uniform prior, the log-density is always zero.}
#' \item{log_grad}{A function for calculating the gradient of the log-density for each element of 
#' a vector. For the uniform prior, this is always 0. With a uniform prior, this gradient is always zero.}
#' @seealso \code{\link{make_logistic_prior}}
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


# Generate all possible binary vectors of length n_nodes
generate_possibilities = function(n_nodes){
  possibilities = expand.grid(replicate(n_nodes, c(0L, 1L), simplify = FALSE))
  as.matrix(possibilities[ , n_nodes:1])
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
nll_grad = function(par, possible_cooc, observed_cooc, n_samples, ...){

  # Find energy for all possible co-occurrence patterns
  E = -c(par %*% possible_cooc)
  
  # Log of the partition function
  logZ = logSumExp(-E)
  
  p = exp(-(E + logZ))
  
  expected_cooc = c(possible_cooc %*% p)
  
  n_samples * (expected_cooc - observed_cooc)
}

# Negative log posterior gradient
nlp_grad = function(
  par, 
  possible_cooc, 
  observed_cooc, 
  n_samples, 
  prior, 
  ...
){
  nll_grad(par, possible_cooc, observed_cooc, n_samples, ...) - prior$log_grad(par)
}

# Find the number of observed co-occurrences (and occurrences)
find_observed_cooc = function(x){
  cp = crossprod(x) / nrow(x)
  
  c(diag(cp), cp[upper.tri(cp)])
}

#' Fit a Markov network to binary data
#' 
#' @description 
#' The optimization is performed using the \code{BFGS} method in the \code{\link[stats]{optim}}
#' function. Note that the likelihood function includes \code{2^N} terms for \code{N} nodes,
#' so parameter estimation can be very slow with more than about 20-25 nodes.
#' 
#' @param x A binary matrix.  Each row corresponds to an observed sample and each
#' column corresponds to one node in the observed network.
#' @param prior    An object of class \code{prior}. The default prior, produced by
#' \code{\link{make_logistic_prior}}, regularizes the parameter estimates. This is
#' important in Markov networks because the number of parameters will often be large
#' compared to the number of independent observations. The flat prior produced by 
#' \code{\link{make_logistic_prior}} is unbiased, but can result in unbounded estimates
#' if two nodes always share the same state in the data set.
#' @param maxit,trace,hessian,...    Arguments passed to \code{\link{optim}}
#' @param parlist  user-specified starting values (optional).
#' @return a \code{list} with the following elements:
#' \item{alpha}{A vector of estimated univariate potentials ("intercepts"). Larger values
#' support higher probabilities for the corresponding node.}
#' \item{beta}{A symmetric matrix of bivariate potentials ("interaction strengths"). Larger values
#' support higher probabilities for the corresponding pair of nodes.}
#' \item{prior}{The prior object that was used during model fitting.}
#' \item{opt}{The list returned by \code{\link[stats]{optim}}.}
#' @examples 
#' # Simulate a random binary matrix with 1000 observations of five binary variables
#' m = matrix(rbinom(5000, size = 1, prob = .5), nrow = 1000, ncol = 5)
#' 
#' fit = rosalia(m)
#' 
#' fit$alpha
#' fit$beta
#' 
#' @export
#' @importFrom assertthat assert_that
rosalia = function(
  x, 
  prior = make_logistic_prior(location = 0, scale = 1), 
  maxit = 100, 
  trace = 1, 
  hessian = FALSE,
  parlist,
  ...
){
  n_nodes = ncol(x)
  n_samples = nrow(x)
  
  assert_that(all(x == 1 | x == 0))
  if(!is.matrix(x)){
    message("coercing x to matrix")
    x = as.matrix(x)
  }
  
  if (n_nodes > 25) {
    message(
      "Note: exact computation with more than 25 nodes is time-consuming and memory-intensive"
    )
  }
  
  possibilities = generate_possibilities(n_nodes = n_nodes)
  
  possible_cooc = sapply(
    1:2^n_nodes,
    function(i){
      tcp = tcrossprod(possibilities[i, ])
      c(diag(tcp), tcp[upper.tri(tcp)])
    }
  )
  
  rm(possibilities) # no longer needed and possibly large drag on memory
  
  if (missing(parlist)) {
    parlist = as.relistable(list(
      alpha = qlogis((colSums(x) + 1) / (nrow(x) + 2)),
      beta = rep(0, choose(n_nodes, 2))
    ))
  }
  
  opt = optim(
    unlist(parlist), 
    fn = nlp, 
    gr = nlp_grad, 
    method = "BFGS", 
    hessian = hessian,
    possible_cooc = possible_cooc,
    rows = find_rows(x),
    observed_cooc = find_observed_cooc(x),
    control = list(trace = trace, maxit = maxit, REPORT = 1, ...),
    n_samples = n_samples,
    prior = prior
  )
  
  betas = relist(opt$par)$beta
  
  beta_matrix = matrix(0, n_nodes, n_nodes)
  beta_matrix[upper.tri(beta_matrix)] = betas
  beta_matrix = beta_matrix + t(beta_matrix)
  
  list(
    alpha = relist(opt$par)$alpha,
    beta = beta_matrix,
    prior = prior,
    opt = opt
  )
}

