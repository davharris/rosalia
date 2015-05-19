test_that("possibilities are correct", {
  n_nodes = 12
  expect_equal(
  strtoi(apply(generate_possibilities(n_nodes), 1, paste0, collapse = ""), 2),
    seq(0, 2^n_nodes - 1)
  )
})

test_that("extracting rows from possibilities works",{
  n_nodes = 12
  n_samples = 55
  x = matrix(rbinom(n_nodes * n_samples, size = 1, prob = 0.5), ncol = n_nodes)
  
  possibilities = generate_possibilities(n_nodes)
  
  rows = find_rows(x)
  
  expect_equal(x,  possibilities[rows, ], check.attributes = FALSE)
})

test_that("loss works with empty theta",{
  n_nodes = 12
  n_samples = 55
  
  x = matrix(rbinom(n_nodes * n_samples, size = 1, prob = 0.5), ncol = n_nodes)
  
  parlist = as.relistable(list(
    upper = rep(0, choose(n_nodes, 2)),
    diagonal = rep(0, n_nodes)
  ))
  
  possibilities = generate_possibilities(n_nodes)
  
  possible_cooc = sapply(
    1:2^n_nodes,
    function(i){
      tcp = tcrossprod(possibilities[i, ])
      c(tcp[upper.tri(tcp)], diag(tcp))
    }
  )
  
  rows = find_rows(x)
  
  expect_equal(
    nll(rows, possible_cooc = possible_cooc, par = unlist(parlist), skeleton = parlist),
    -log(1/2^n_nodes) * length(rows)
  )
})

test_that("loss works on a toy example", {
  n_nodes = 2
  n_samples = 55
  
  theta = matrix(c(1, -1, -1, 2), nrow = 2)
  par = c(upper = -1, diagonal = c(1, 2))
  
  
  # Known values calculated by hand in Figure X of associated paper
  correct_E = c(0, -2, -1, -2)
  correct_logZ = logSumExp(-correct_E)
  
  
  # Simulate 4 possibilities
  possibilities = generate_possibilities(n_nodes)
  possible_cooc = sapply(
    1:2^n_nodes,
    function(i){
      tcp = tcrossprod(possibilities[i, ])
      c(tcp[upper.tri(tcp)], diag(tcp))
    }
  )
  
  # Simulate small landscape
  x = possibilities[sample.int(nrow(possibilities), n_samples, replace = TRUE, prob = exp(-correct_E)), ]
  
  correct_p = exp(-correct_E[find_rows(x)]) / exp(correct_logZ)
  correct_nll = -sum(log(correct_p))
  expect_equal(
    nll(par = par, rows = find_rows(x), possible_cooc = possible_cooc),
    correct_nll
  )
})

test_that("pair likelihood gradient is correct", {
  n_nodes = 3
  n_samples = 1E3
  
  pre_theta = matrix(rnorm(9), ncol = 3)
  theta = pre_theta + t(pre_theta)
  
  x = matrix(rbinom(n_samples * n_nodes, size = 1, prob = .5), ncol = n_nodes)
  
  eps = 1E-5
  
  possibilities = generate_possibilities(n_nodes = n_nodes)
  
  possible_cooc = sapply(
    1:2^n_nodes,
    function(i){
      tcp = tcrossprod(possibilities[i, ])
      c(diag(tcp), tcp[upper.tri(tcp)])
    }
  )
  
  
  diff = nll(par = c(0, eps, 0, 0, 0, 0), 
             rows = find_rows(x),
             possible_cooc = possible_cooc
  ) - nll(par = c(0, -eps, 0, 0, 0, 0), 
          rows = find_rows(x),
          possible_cooc = possible_cooc
  )
  
  observed_cooc = find_observed_cooc(x)
  
  grad = nll_grad(par = c(0, 0, 0, 0, 0, 0), 
                  rows = find_rows(x),
                  possible_cooc = possible_cooc,
                  observed_cooc = observed_cooc,
                  n_samples = n_samples
  )
  
  
  expect_equal(
    diff / 2 / eps,
    grad[2],
    check.attributes = FALSE
  )
})

test_that("diagonal likelihood gradient is correct", {
  n_nodes = 3
  n_samples = 1E3
  
  pre_theta = matrix(rnorm(9), ncol = 3)
  theta = pre_theta + t(pre_theta)
  
  x = matrix(rbinom(n_samples * n_nodes, size = 1, prob = .5), ncol = n_nodes)
  
  eps = 1E-5
  
  possibilities = generate_possibilities(n_nodes = n_nodes)
  
  possible_cooc = sapply(
    1:2^n_nodes,
    function(i){
      tcp = tcrossprod(possibilities[i, ])
      c(diag(tcp), tcp[upper.tri(tcp)])
    }
  )
  
  
  diff = nll(par = c(1, 1, 1, 1, 1, 1 + eps), 
             rows = find_rows(x),
             possible_cooc = possible_cooc
  ) - nll(par = c(1, 1, 1, 1, 1, 1 - eps), 
          rows = find_rows(x),
          possible_cooc = possible_cooc
  )
  
  observed_cooc = find_observed_cooc(x)
  
  grad = nll_grad(par = c(1, 1, 1, 1, 1, 1), 
                  rows = find_rows(x),
                  possible_cooc = possible_cooc,
                  observed_cooc = observed_cooc,
                  n_samples = n_samples
  )
  
  expect_equal(
    diff / 2 / eps,
    grad[6],
    check.attributes = FALSE
  )
})

test_that("optimization works on small data", {
  n_nodes = 2
  n_samples = 1E6
  
  theta = matrix(c(1, -1, -1, 2), nrow = 2)
  
  # Known values calculated by hand in Figure X of associated paper
  correct_E = c(0, -2, -1, -2)
  correct_logZ = logSumExp(-correct_E)
  
  p = exp(-correct_E) / exp(correct_logZ )
  
  
  # Simulate 4 possibilities
  possibilities = generate_possibilities(n_nodes)
  
  rows = sample.int(nrow(possibilities), prob = p, size = n_samples, replace = TRUE)
  x = possibilities[rows, ]
  
  out = rosalia(x, maxit = 1E4)

  estimated_theta = out$beta
  diag(estimated_theta) = out$alpha
  
  expect_equal(
    theta,
    estimated_theta,
    tolerance = 0.01
  )
})

test_that("logistic prior works", {
  x = c(pi, log(7))
  eps = 1E-6
  s = c(exp(2), sqrt(12))
  
  diff = dlogis(x + eps, scale = s, log = TRUE) - 
    dlogis(x - eps, scale = s, log = TRUE)
  
  expect_equal(
    make_logistic_prior(scale = s)$log_grad(x),
    (diff) / 2 / eps
  )
})

