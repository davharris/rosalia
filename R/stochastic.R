logistic = binomial()$linkinv

rbern = function(p){rbinom(length(p), size = 1, prob = p)}

#' @importFrom progress progress_bar
f = function(
  x, 
  y, 
  maxit, 
  intercept_prior, 
  alpha_x_prior, 
  beta_prior, 
  initial_learning_rate, 
  rate_decay
){
  n_nodes = ncol(y)
  n_samples = nrow(y)
  
  # Calculate sufficient statistics for the observed data
  y_stats = crossprod(y)
  xy_stats = t(x) %*% y
  
  alpha_x = delta_alpha_x = matrix(0, nrow = n_env, ncol = n_spp)
  intercept = delta_intercept = qlogis(colMeans(y))
  beta = delta_beta = matrix(0, nrow = n_nodes, ncol = n_nodes)
  alpha = matrix(0, nrow = n_nodes, ncol = n_nodes) # no delta alpha to initialize
  # b/c alpha not optimized directly
  
  pb = progress_bar$new(
    format = "  Fitting [:bar] :percent eta: :eta",
    total = maxit,
    clear = FALSE
  )
  
  # Initialize simulated matrix
  y_sim = matrix(0.5, nrow = nrow(y), ncol = ncol(y))
  
  for(i in 1:maxit){
    
    #### Update learning rates ####
    learning_rate = initial_learning_rate * rate_decay / (rate_decay + i)
    
    #### Gibbs sampling ####
    # Update alpha
    alpha = x %*% alpha_x %plus% intercept
    
    # Sample entries in y_sim from their conditional distribution (Gibbs sampling)
    for(j in sample.int(n_nodes)){
      y_sim[,j] = rbern(logistic(alpha[ , j] + y_sim %*% beta[ , j]))
    }
  }
  
  #### Stochastic gradients ####
  # Calculate sufficient statistics
  y_sim_stats = crossprod(y_sim)
  xy_sim_stats = t(x) %*% y_sim
  
  # Calculate the gradient with respect to alpha and beta
  stats_difference = y_stats - y_sim_stats
  intercept_grad = (diag(stats_difference) + alpha_prior$log_grad(alpha_x)) / n_samples
  beta_grad = (stats_difference + beta_prior$log_grad(beta)) / n_samples
  diag(beta_grad) = 0
  
  xy_difference = xy_stats - xy_sim_stats
  alpha_x_grad = (xy_difference + alpha_x_prior$species(alpha_x))  / n_samples
  
  
  #### Calculate parameter updates ####
  delta_beta = beta_grad * learning_rate / 10 + 
    momentum * delta_beta
  delta_alpha_species = alpha_species_grad * learning_rate + 
    momentum  * delta_alpha_species
  delta_alpha_env = alpha_env_grad * learning_rate + 
    momentum  * delta_alpha_env
  
  
  beta = beta + delta_beta
  alpha_species = alpha_species + delta_alpha_species
  alpha_env = alpha_env + delta_alpha_env
}


