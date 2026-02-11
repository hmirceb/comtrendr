trend_doak  <- function(x, time = NULL){
  # Make time variable if missing
  if ( is.null(time) ) {
    time <- seq_len(length(x))
  }
  # Get good years
  g_y <- ( x >= 0 | !is.na(x) )
  # Keep only good years
  x <- x[g_y]
  # Add 1 percent of median abundance
  x <- x + (0.001 * median(x))
  # Keep only good years
  time <- time[g_y]
  # Square root of difference between time steps
  d_t <- sqrt(diff(time))
  # Log ratios between years
  d_n <- diff(log(x))
  # Average log ratios
  y <- d_n / d_t
  # Lm for log ratios and corrected time
  t_m <- lm(y~0+d_t)
  # Stochastic variance of growth rates
  stoc_v <- sigma(t_m)^2
  # Return mean growth rate
  return(coef(t_m))
}

trend_logl  <- function(x, time = NULL){
  # Make time variable if missing
  if ( is.null(time) ) {
    time <- seq_len(length(x))
  }
  # Get good years
  g_y <- ( x >= 0 | !is.na(x) )
  # Keep only good years
  x <- x[g_y]
  # Add 1 percent of median abundance
  x <- x + (0.001 * median(x))
  # Keep only good years
  time <- time[g_y]
  
  d_t <- time[g_y]-min(time[g_y])
  d_n <- log(x)
  t_m <- lm(d_n~d_t)
  return(coef(t_m)[2])
}

mean(diff(log(d)))
