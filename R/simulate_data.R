simulate_data = function(par, possibilities, possible_cooc, n_sites){
  # Assemblage probabilities under specified model --------------------------
  
  # Energies of all communities
  E = colSums(-par * possible_cooc)
  
  # Partition function
  Z = logSumExp(-E)
  
  # Probabilities of each community
  p = exp(-E) / exp(Z)
  
  # generate x --------------------------------------------------------------
  
  # Randomly select row indices from the distribution described above
  rows = sample.int(nrow(possibilities), size = n_sites, prob = p, replace = TRUE)
  
  # Return the selected rows
  possibilities[rows, ]
}
