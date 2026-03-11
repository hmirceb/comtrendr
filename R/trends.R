#' Estimate mean annual growth rate following Dennis et al. (2001). 
#'
#' This function estimates the mean annual growth rate of a time series of abundance data following Dennis et al. (2001) (see Details).
#'
#' @param x Numeric. A vector of abundances.
#' @param time Numeric. A vector with the time steps corresponding to each value in x.
#'
#' @details For a given a time series of abundances \eqn{n_{t}} the function estimates its mean growth rate \eqn{\mu} and variance \eqn{\sigma^2} using a linear regression model without intercept such as:
#' \deqn{y_{i} \sim \mu t_{i} + \epsilon_{i} }
#' \deqn{\epsilon_{i} \sim Normal(0, \sigma^2)}
#'  where \eqn{ y_{i} = \dfrac{ \ln{ ( n_{t}/n_{t-1}) } }{ t_{i} } } and \eqn{ t_{i} = \sqrt{ t_{t}-t_{t-1} } }
#'  
#'  Note that the confidence intervals for the estimate of \eqn{\mu} are not reliable.
#'  
#' @returns A named list with the mean annual growth rate in the natural logarithm scale and its variance.
#' 
#' @references
#' - Dennis, B., Munholland, P. L., & Scott, J. M. (1991). Estimation of Growth and Extinction Parameters for Endangered Species. Ecological Monographs, 61(2), 115–143.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @examples
#' require(detrending)
#' 
#' # Simulate time series
#' ts <- 5^seq(1, 2, by = 0.01)
#' mean(diff(log(ts))) # True trend = 0.016 (~1.62%)
#' # Simulate some random noise
#' noise <- rnorm(length(ts))
#' 
#' # Estimate trend
#' trend_dennis(ts+noise)
#' @export
trend_dennis  <- function(x, time = NULL){
  # Make time variable if missing
  if ( is.null(time) ) {
    time <- seq_along(x)
  }
  # Get good years
  g_y <- ( x >= 0 | !is.na(x) )
  # Keep only good data
  x <- x[g_y]
  # Keep only good times
  time <- time[g_y]
  # Square root of difference between time steps
  d_t <- sqrt(diff(time))
  # Log ratios between time steps
  d_n <- diff(log(x))
  # Average log ratios
  y <- d_n / d_t
  # Lm for log ratios and corrected time
  t_m <- stats::lm(y~0+d_t)
  # Stochastic variance of growth rates
  stoc_v <- stats::sigma(t_m)^2
  # Return mean growth rate, variance and confidence interval
  res <- c(stats::coef(t_m), stoc_v, stats::confint(t_m))
  names(res) <- c("trend", "var", "l95", "u95")
  return(res)
}

#' Estimate mean annual growth rate using linear regression 
#'
#' @param x Numeric. A vector of abundances.
#' @param time Numeric. A vector with the time steps corresponding to each value in x.
#'
#' @details For a given a time series of abundances \eqn{n_{t}} the function estimates a linear regression model with log-transformed abundances as the response variable and time steps (with \eqn{t_{0} = 0, t_{1} = 1...}) as the explanatory variable.
#' 
#'
#' @returns A numeric value with the slope of a linear regression of the log transformed abundances.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @examples
#' require(detrending)
#' 
#' # Simulate time series
#' ts <- 5^seq(1, 2, by = 0.01)
#' mean(diff(log(ts))) # True trend = 0.016 (~1.62%)
#' # Simulate some random noise
#' noise <- rnorm(length(ts))
#' 
#' # Estimate trend
#' trend_loglinear(ts+noise)
#' @export
trend_loglinear  <- function(x, time = NULL){
  # Make time variable if missing
  if ( is.null(time) ) {
    time <- seq_len(length(x))
  }
  # Get good years
  g_y <- ( x >= 0 | !is.na(x) )
  # Keep only good years
  x <- x[g_y]
  # Keep only good years
  time <- time[g_y]
  # Set time as increasing from 0
  d_t <- time[g_y]-min(time[g_y])
  # log transform
  d_n <- log(x)
  # linear regression
  t_m <- stats::lm(d_n~d_t)
  # Return slope and confidence interval
  res <- c(stats::coef(t_m)[2], stats::confint(t_m)[2,])
  names(res) <- c("trend", "l95", "u95")
  
  return(res)
}


#' Estimate population trends in a community
#'
#' @param x A data.frame. A community matrix of species abundances with time in rows and taxa in columns. Optionally it can include community and time columns. 
#' @param time_col Character. Name of the column with time variable. Optional with default "time".
#' @param method Character. Method to estimate the trends, one of "dennis" or "loglinear". Default "dennis".
#' @param plot Boolean. Plot species abundances and their estimated trends. Default FALSE. 
#'
#' @returns A data.frame with the trend (in the natural logarithm scale) for each species in the community along with its variance and 95% confidence interval.
#' 
#' @examples
#' require(detrending)
#' 
#' # Simulate community data with trends
#' comm_df <- sim_mvcomm(switch_trend = T, bimodal_trend = T)
#' 
#' # Estimate trend for each species and plot estimated trends
#' comm_trend(comm_df$sim_data, method = "loglinear", plot = T)
#' @export
comm_trend <- function(x, time_col = "time", method = "loglinear", plot = FALSE){
  # Match variance function
  trend_func <- switch(
    method,
    dennis = trend_dennis,
    loglinear = trend_loglinear
  )
  
  # Check if a time column was specified for detrending methods and order rows
  x <- check_time(x, time_col = time_col, term = "var", rm = TRUE)
  
  # Replace NAs with 0 and remove columns (species) with 0 abundance across all years 
  x <- remove_empty_sps(x = x, time_col = time_col)
  
  trends <- as.data.frame(
    cbind(
      taxa = colnames(x), 
      as.data.frame(
        do.call("rbind", 
                apply(x, MARGIN = 2, trend_func, simplify = F)
        )
      )
    )
  )
  rownames(trends) <- NULL
  
  # Plot abundances
  if (plot) {
    graphics::par(mfrow = c(1,2))
    plot_com(x)
    
    plot(x = trends[1,]$trend, y = seq_along(trends$taxa)[1],
         xlim = c(min(trends$l95), max(trends$u95)),
         ylim = c(min(seq_along(trends$taxa)), max(seq_along(trends$taxa))),
         pch = 19,
         col = 1,
         xlab = "trend (log)",
         ylab = "taxa", 
         yaxt = "n")
    graphics::arrows(x0 = trends$l95, x1 = trends$u95, y0 = seq_along(trends$taxa),
           code = 3, length = 0.05, angle = 90)
    for (i in 2:nrow(trends)) {
      graphics::points(x = trends[i,]$trend, y = seq_along(trends$taxa)[i], 
                       pch = 19,
                       col = i)
    }
    
    graphics::axis(2, at = seq_along(trends$taxa), labels = trends$taxa, las = 2)
    graphics::abline(v = 0, lty = "dashed")
    
    graphics::par(mfrow = c(1,1), 
        xpd=FALSE, 
        mar=c(5.1, 4.1, 4.1, 2.1))
  }
  
  return(trends)
}