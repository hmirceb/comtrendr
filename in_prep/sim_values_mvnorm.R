#' Title
#'
#' @param n_sp 
#' @param years 
#' @param tot_abu 
#' @param power 
#' @param bound_pos 
#' @param corr 
#' @param p 
#' @param trend 
#' @param trend_sd 
#'
#' @returns
#' @export
#'
#' @examples
simmvnrom <- function(n_sp = 10,
         years = 25,
         tot_abu = 200 * n_sp,
         power = 1.6, # values closer to 1 indicate more stable dominants and thus more dominance effect
         bound_pos = TRUE,
         corr = 0.5, # correlation between time series
         p = 0.8, # Higher values mean less dominance
         trend = FALSE,
         trend_sd = 0.01)  {
  
  # Vector of mean abundances. First get a vector of relative abundances 
  # that add up to 1 using the dirchlet distribtion. Parameter alpha (p)
  # controls the spread of the values, with higher values leading to 
  # more even relative abundances and thus less dominance
  mean_abu <- sort(tot_abu * gtools::rdirichlet(1, alpha = rep(p, n_sp))[1,], decreasing = TRUE)
  
  # Create a matrix of abundances
  abu_matrix <- matrix(rep(mean_abu, times = years),
                       nrow = years, ncol = n_sp, byrow = TRUE)
  
  # # Draw growth rates for each species
  # if ( isTRUE(trend) ) {
  #   av_grates <- rnorm(n_sp, mean = 0, sd = 1)
  #   g_rates <- sapply(av_grates, function(r) {
  #     rnorm(years-1, mean = r, sd = trend_sd)
  #   })
  # } else {
  #   g_rates <- matrix(0, ncol = n_sp, nrow = years)
  # }
  # 
  # # Simulate trends of each species using Ricker model to avoid overshooting populations
  # for (t in 2:years) {
  #   abu_matrix[t,] <- exp( log( abu_matrix[t-1,] ) + g_rates[t-1,])
  #   
  #   idx <- abu_matrix[t, ] > mean_abu
  #   abu_matrix[t, idx] <- mean_abu[idx]
  # }
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
    abi <- MASS::mvrnorm(n = 1, mu = abu_matrix[j,], Sigma = ss)
    # Force positive values
    if (bound_pos) {
      abi[abi < 0] <- 0
    }
    simcom[j, ] <- abi
  }
  
  res <- list(sim_data = as.data.frame(simcom),
              av_trend = colMeans(apply(log(simcom),1,diff)))
  return(res)
} 

d <- simmvnrom(trend_sd = 0.1, trend = T)
plot_com(d$sim_data)
d$av_trend


apply(d$sim_data, 2, trend_doak, 1:25)

params <- expand.grid(n_sp = seq(from = 10, to = 30, by = 2),
            power = seq(from = 1.2, to = 1.8, by = 0.2),
            corr = seq(from = 0.1, to = 0.8, by = 0.2),
            p = seq(from = 0.1, to = 0.8, by = 0.2))

niter <- 1
sims <- apply(params, 1, function(y, nsamp = nsamp) {
  n_sp <- y[[1]]
  power <- y[[2]]
  corr <- y[[3]]
  p <- y[[4]]
  
  jj <- list()
  for (t in seq_len(niter)) {
    d <- simmvnrom(n_sp = n_sp,
                   power = power,
                   corr = corr,
                   p = p,
                   bound_pos = TRUE,
                   years = 20,
                   tot_abu = 100 * n_sp,
                   trend = FALSE)
    
    b <- comstab_term(as.data.frame(d$sim_data))
    c <- data.frame(delta_rel = b$Relative["Delta_cont"],
                    psi_rel = b$Relative["Psi_cont"],
                    omega_rel = b$Relative["omega_cont"],
                    n_sp = n_sp, power = power, corr = corr, p = p)
    jj[[t]] <- c
    
  }
  jj <- do.call("rbind", jj)
  return(jj)
}
)
sims2<- do.call("rbind", sims)

sims2 %>% 
  mutate(id = row_number()) %>% 
  pivot_longer(cols = delta_rel:omega_rel) %>% 
  ggplot(aes(x = n_sp, y = value))+
  geom_point()+
  facet_wrap(~name)

sims2 %>% 
  pivot_longer(cols = delta_rel:omega_rel) %>% 
  ggplot(aes(x = power, y = value))+
  geom_point()+
  facet_wrap(~name)

sims2 %>% 
  pivot_longer(cols = delta_rel:omega_rel) %>% 
  ggplot(aes(x = p, y = value))+
  geom_point()+
  facet_wrap(~name)

sims2 %>% 
  mutate(id = row_number()) %>% 
  pivot_longer(cols = delta_rel:omega_rel) %>% 
  ggplot(aes(x = id, y = value, group = id, fill = name))+
  geom_bar(stat = "identity")+
  theme_void()

