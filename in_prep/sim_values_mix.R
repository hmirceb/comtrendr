aaa <- function(n_sp = 10,
years = 25,
tot_abu = 200 * n_sp,
power = 1.6,# values closer to 1 indicate more stable dominants and thus more dominance effect
bound_pos = TRUE,
corr = 0.5, # correlation between time series
p = 0.8, # Higher values mean less dominance
trend = FALSE,
trend_sd = 0.01,
switch_trend = "on",
mean_trend_resp = 1,
sd_trend_resp = 1,
bimodal_trend = FALSE){
  
  # Vector of mean abundances. First get a vector of relative abundances 
  # that add up to 1 using the dirchlet distribtion. Parameter alpha (p)
  # controls the spread of the values, with higher values leading to 
  # more even relative abundances and thus less dominance
  mean_abu <- sort(tot_abu * gtools::rdirichlet(1, alpha = rep(p, n_sp))[1,], decreasing = TRUE)
  
  # Create a matrix of abundances
  abu_matrix <- matrix(rep(mean_abu, times = years),
                       nrow = years, ncol = n_sp, byrow = TRUE)
  
  trend = seq(-1, 1, length.out = years)
  trend_resp <- response(switch_trend,
                         n_sp,
                         mean = mean_trend_resp,
                         sd = sd_trend_resp,
                         bimodal = bimodal_trend)
  
  # Get Sd of abundances from TPL
  sd_abu <- sqrt(mean_abu ** power)
  
  # Variance-covariance matrix for MVN distribution
  k <- 3
  # This makes positive definitive matrix (necessary for MVN) 
  # with correlation between species
  Lambda <- matrix(rnorm(n_sp * k, sd = sqrt(corr)), n_sp, k)
  Psi <- diag(sd_abu^2 * (1 - corr)) 
  ss <- Lambda %*% t(Lambda) + Psi
  
  # Simulate random variation around mean cover for each species
  # drawn from multivariate normal so species can correlate
  simcom <- matrix(0, years, n_sp)
  for (j in 1:years) {
    abi <- MASS::mvrnorm(n = 1, 
                         mu = abu_matrix[j,] * (1 + trend[j] * trend_resp), 
                         Sigma = ss)
    # Force positive values
    if (bound_pos) {
      abi[abi < 0] <- 0
    }
    simcom[j, ] <- abi
  }
  
  res <- list(sim_data = as.data.frame(simcom),
              av_trend = colMeans(apply(log(simcom),1,diff)))
}

par(mfrow = c(1,2))
plot_com(aaa(years = 25, corr = 0, bimodal_trend = T)$sim_data)

plot_com(syngenr(n_sp = 10, years = 25,
        tot_abu = 200 * n_sp,
        power = 1.6,# values closer to 1 indicate more stable dominants and thus more dominance effect
        bound_pos = TRUE,
        switch_trend = "on",
        mean_trend_resp = 1,
        sd_trend_resp = 1,
        bimodal_trend = T)$time_species_matrix)
