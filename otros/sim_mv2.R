#' Simulate community from a multivariate normal
#'
#' @param n_sp Numeric. Number of species in the community.
#' @param years Numeric. NUmber of years to simulate.
#' @param tot_abu Numeric. Total abundance of the community.
#' @param power Numeric. Exponent of the Taylor's Power Law to estimate variance from mean abundance.
#' @param corr Numeric. Average correlation between populations.
#' @param even Numeric. The relative abundance (between 0 and 1) of the most abundant species. Controls the evenness of the community with lower values indicating more even communities. Alternatively, a vector of relative abundaces of length = n_sp.
#' @param trend_mean Numeric. Mean of the trend. Can be a single value (for a shared trend across species) or a vector of length = `n_sp` (for individual trends). Positive values indicate growth and negative ones, decline. Default 0 (no trend).
#' @param trend_sd Numeric. Standard deviation of the trend. Can be a single value or a vector of the same length as `trend_mean`.
#' 
#' @returns A named list with three elements:
#' - `sim_data`: A data.frame with the simulated data including any possible trends. Species in columns and time steps in rows.
#' 
#' - `baseline`: A data.frame with the simulated data without any trends. Species in columns and time steps in rows.
#' 
#' - `true_trend`: A named vector with the true mean trends of each simulated species. 
#' 
#' - `params`: A named vector with the parameters used to simulate the data.
#' 
#' @author Lars Götzenberger, \email{jsegrestin@@gmail.com}
#' @author Jan Lepš, \email{suspa@@prf.jcu.cz}
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @examples
#' require(comtrendr)
#' 
#' sim_mvcomm(n_sp = 15, years = 30)
#' @export
sim_mvcomm <- function(n_sp = 10,
                       years = 25,
                       tot_abu = 200 * n_sp,
                       power = 1.8,
                       corr = 0,
                       even = 0.5,
                       trend_mean = 0,
                       trend_sd = 0.01) {
  
  # Vector of mean abundances.
  # check evenness values
  if (length(even) == 1) {
    even <- ifelse(even == 0, even + 0.01, even) # values == 0 give error
    mean_abu <- sort(
      tot_abu * geom_seq(max_rel_abu = even, n_sp = n_sp),
      decreasing = TRUE
    )
  } else {
    if (length(even) != n_sp) {
      stop("The length of the vector of relative abundances and the number of species differ.")
    }
    if (sum(even) != 1) {
      stop("The vector of relative abundances do not add up to 1.")
    }
    mean_abu <- sort(
      tot_abu * even,
      decreasing = TRUE
    )
  }
  
  # Simulate trends 
  # Check vector of trends
  if (length(trend_mean) == 1) {
    trend_resp <- stats::rnorm(n = n_sp, mean = trend_mean, sd = trend_sd)
  } else {
    if (length(trend_mean) != length(trend_sd)) {
      stop("Lengths of vectors of means and SD differ.")
    }
    trend_resp <- sapply(seq_along(trend_mean), FUN = function(z) {
      stats::rnorm(n = 1, mean = trend_mean[z], sd = trend_sd[z])
    })
  }
  
  # Check that correlation is feasible
  eta_min  <- -1 / (n_sp - 1)
  n_sp_max <- ceiling((-1 / corr) + 1)
  if (corr < eta_min | corr > 1) {
    stop(paste0("correlation value must be between ", round(eta_min, 3),
                " and 1 or n_sp lower than ", n_sp_max))
  }
  
  # Simulate random variation around mean abundance for each species
  # drawn from multivariate normal so species correlate
  simcom <- matrix(0, years, n_sp)
  for (j in 1:years) {
    mu <- mean_abu * exp(trend_resp * j)
    mu <- pmax(mu, 0.01)
    # Get SD of abundances from TPL
    sd_abu <- sqrt(mu ** power)
    abi <- unlist(faux::rnorm_multi(n = 1,
                                    mu = mu,
                                    sd = sd_abu,
                                    r = corr))
    # Force positive values
    abi <- pmax(abi, 0)
    simcom[j, ] <- abi
  }
  
  # Add a small offset (1% of the mean abundance of each species) to avoid having 0s
  off    <- colMeans(simcom) * 0.01
  simcom <- as.data.frame(sweep(x = simcom, MARGIN = 2, STATS = off, FUN = "+"))
  
  p <- matrix(trend_resp, ncol = n_sp, nrow = years, byrow = T)
  p <- sweep(p, 1, seq_len(years), "*")
  baseline <- simcom / exp(p)
  
  # Set species names
  colnames(simcom) <- paste(sep = "_", "sp", seq_len(n_sp))
  colnames(baseline) <- paste(sep = "_", "sp", seq_len(n_sp))
  
  # Results into list
  res <- list(
    sim_data = simcom,
    baseline = baseline,
    true_trend = colMeans(apply(log(simcom), 2, diff)),
    params = c(n_sp = n_sp,
               years = years,
               tot_abu = tot_abu,
               power = power,
               corr = corr,
               even = even,
               trend_mean = unique(trend_mean),
               trend_sd = unique(trend_sd))
  )
}
