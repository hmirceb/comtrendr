# Functions to decompose community stability

#' Decompose community stability
#'
#' `comstab_term()` partitions the temporal coefficient of variation of a community into the variability of the average species and three stabilizing effects: the dominance, asynchrony and averaging effects. It allows standard estimates of variance and CV as well as their detrended versions using Hill's two and three term local quadratic variance estimates (see Details).
#' 
#' @usage comstab_term(x, term = "var", ...)
#' 
#' @param x A data.frame. A community matrix of abundances with time in rows and taxa in columns. Optionally it can include community and time columns. 
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#' @param community_col Character. Name of the column with the community identifier. Optional with default "comm".
#' @param time_col Character. Name of the column with time variable. Optional with default "time".
#' 
#' @details The analytic framework is described in details in Segrestin *et al.* (2024).
#' In short, the partitioning relies on the following equation: \deqn{CV_{com} = CV_e \Delta \Psi \omega} 
#' where \eqn{CV_{com}} is the community coefficient of variation (reciprocal of community stability), 
#' \eqn{CV_e} is the expected community CV when controlling for the dominance structure and species temporal synchrony,
#' \eqn{ \Delta} is the dominance effect, \eqn{ \Psi} is the asynchrony effect, and \eqn{ \omega} is the averaging effect.
#' 
#' @returns An object of class `'comstab'`, a list of named vectors containing the following components:
#'  * `'CVs'` a named vector of calculated coefficient of variations. `CVe` is the CV of an average species,
#'  `CVtilde` is the mean of species CVs weighted by their relative abundances, `CVa` is the expected community CV if 
#'   the community was stabilized by species asynchrony only, and `CVc` is the observed community CV.
#'  * `'Stabilization'` a named vector of the stabilizing effects. `tau` is the total stabilization, `Delta` is
#'  the dominance effect, `Psi` is the asynchrony effect, and `omega` is the averaging effect.
#'  * `'Relative'` a named vector of the relative contributions of each stabilizing effect to the total stabilization.
#'  `Delta_cont`, `Psi_cont`, and `omega_cont` are the relative contribution of respectively, the dominance, asynchrony, and averaging effects to the total stabilization.
#'  Returns a vector of NAs if any Stabilizing effect is higher than 1.
#'  
#' @references
#'   * Segrestin, J., Götzenberger, L., Valencia, E., de Bello, F., & Lepš, J. (2024). A unified framework for partitioning the drivers of stability of ecological communities. Global Ecology and Biogeography, 33(5), e13828.
#' 
#' @examples
#' 
#' @author Jules Segrestin, \email{jsegrestin@@gmail.com}
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @export
comstab_term = function(x, 
                        term = "var",
                        community_col = "comm",
                        time_col = "time") {
  # Match argument for variance function to use
  # Table of options
  options <- data.frame(term = c("var", "two", "three", "linear"),
                        var = c("var", "var_t2", "var_t3", "var_linear"))
  
  # Get choice
  opt <- options[options$term == term,]
  
  # Match functions
  var_func <- match.fun(opt$var) # variance
  
  # Errors if data is not properly formated
  # if ( !is.matrix(x) | !is.data.frame(x)) {
  #   stop("Error: x is not a matrix")
  # }
  # if ( !is.numeric(x) ) {
  #   stop("Error: non-numerical values in x")
  # } 
  # if ( any(x < 0) ) {
  #   stop("Error: negative values in x")
  # } 
  if ( nrow(x) == 1 ) {
    stop("Only one year provided.")
  }
  
  # Check if a time column was specified for detrending methods
  if ( !time_col %in% colnames(x) ) {
    x <- cbind(time = seq_len(nrow(x)), x)
    colnames(x)[1] <- time_col
  } 
  if ( !time_col %in% colnames(x) & 
       term %in% c("two", "three") ) {
    warning("Missing time column. Rows are assumed to be in order for detrending.")
  }
  # Order by year if time column provided
  x <- x[with(x, order(x[, time_col])),]
  
  # Remove time column once df is ordered
  id_cols <- colnames(x) %in% c(community_col, time_col)
  x <- x[,!id_cols]
  
  # Replace NAs with 0
  x[is.na(x)] <- 0
  
  # Remove columns (species) with 0 abundance across all years 
  x <- x[, colSums(x) > 0, drop = FALSE]
  
  # Remove columns (species) with constant abundance (min abundance == max abundance)
  x <- x[, apply(X = x, MARGIN = 2, FUN = min) != 
           apply(X = x, MARGIN = 2, FUN = max), 
         drop = FALSE]
  
  # Turn into matrix
  x <- as.matrix(x)
  
  # number of species
  n <- ncol(x) 
  
  ## Community metrics ##
  varsum <- var_func(rowSums(x)) # variance of sum of abundances
  meansum <- mean(rowSums(x)) # mean of sum of abundances
  CV <- sqrt(varsum) / meansum # CV of sum of abundances
  
  # Check if community fluctuates
  if (CV == 0) {
    stop("The community CV is zero. This analysis does not apply to \n perfectly stable communities.")
  }
  
  # Stop if there is only one species
  if (ncol(x) == 1) {
    stop("This analysis is not relevant for single-species communities.")
  }
  
  #--------------#
  # Partitioning #
  #--------------#
  
  # Estimate CVe from TPL of all species
  vari <- apply(X = x, MARGIN = 2, FUN = var_func) # Variance of each species
  meani <- colMeans(x) # Mean of each species
  CVi <- sqrt(vari) / meani  # CV of each species
    
  if (any(CVi == 0)) { # Warn if constant species present
    warning("Non-fluctuating species found in the data.")
    }
    
  CV0 <- which(CVi > 0) # Use only species with CV != 0
  TPL <- stats::coef(stats::lm(log10(CVi[CV0]) ~ log10(meani[CV0]))) # LM of CVs and means on log scale
  CVe <- 10^TPL[1] * (mean(x)^TPL[2]) # Predict CVe from mean abundance and LM (backtransformed from log scale)
    
  # Test correlation between individual CVs and mean abundances (if there are more than 5 species with variation)
  if (sum(CV0) > 5) {
    testcor <- stats::cor.test(log10(CVi[CV0]), log10(meani[CV0]))$p.value > 0.05
    if (testcor) {
      warning("No significant power law between species CVs and abundances.")
      }
  } else {
    warning("Low number of species. The power law between species CVs and abundances cannot be tested.")
  }
    
  ## Dominance effect #
  sumsd <- sum(sqrt(vari)) # sum of individual SDs
  CVtilde <- sumsd / meansum # CV tilde. Weighted mean of individual CVs. sum(pi * sdi/mui) = sum(mui/meansum * sdi/mui) = sum(sdi/meansum) = sum(sdi)/meansum
  Delta <- CVtilde / CVe # Ratio CVtilde / CVe "average species"
  if (Delta > 1) {
    warning("Destabilizing effect of dominants. Relative effects cannot be computed.")
    }
  
  ## Compensatory dynamics ##
  sdsum <- sqrt(varsum) # Square root of sum of variances (SD). Equivalent to SD of the sum of yearly abundances (whole community)
  rootPhi <- sdsum / sumsd # Ratio between SD of whole community vs sum of individual SDs
  sumvar <- sum(vari) # Sum of individual variances
  alpha <- log10(1/2) / (log10(sumvar/(sumsd^2))) # Scaling coefficient (eq. 7)
  Psi <- rootPhi^alpha # Asynchrony effect
  omega <- rootPhi / Psi # Diversity effect
  if (omega > 1) {
    warning("Community diversity is lower than the null diversity. Relative effects cannot be computed.")
  }
  
  ## Partitioning ##
  tau <- Delta * Psi * omega
  CVs <- stats::setNames(object = c(CVe, CVtilde, CVtilde * 
                                      Psi, CV), nm = c("CVe", "CVtilde", "CVa", "CVc"))
  Stabilization <- stats::setNames(object = c(tau, Delta, 
                                              Psi, omega), nm = c("tau", "Delta", "Psi", "omega"))
  if (any(Stabilization > 1)) {
    Relative <- stats::setNames(object = rep(NA, 3), 
                                nm = c("Delta_cont", "Psi_cont", "omega_cont"))
  } else { # Return relative importance of each component
    Relative <- stats::setNames(object = c(log10(Delta) / log10(tau), # dominance
                                           log10(Psi) / log10(tau), # asynchrony
                                           log10(omega) / log10(tau)), # averaging 
                                nm = c("Delta_cont", "Psi_cont", "omega_cont"))
  }
  # Results into a list
  res <- list(CVs = CVs, Stabilization = Stabilization, 
              Relative = Relative)
  class(res) <- "comstab"
  return(res)
}

#' Print method for comstab objects
#'
#' @export
print.comstab <- function(x, ...){
  cat("\nPartitionning of the community temporal variability (CV)")
  cat("\n")
  cat(paste0("Community CV = ", round(x$CVs["CVc"], 2),
             "\nTotal stabilization = ", round(x$Stabilization["tau"], 2),
             "\nDominance effect = ", round(x$Stabilization["Delta"], 2),
             "\nAsynchrony effect = ", round(x$Stabilization["Psi"], 2),
             "\nAveraging effect = ", round(x$Stabilization["omega"], 2)))
  cat("\n")
  cat("\nRelatives contributions:")
  cat(paste0("\n% Dominance = ", round(x$Relative["Delta_cont"], 2)),
      paste0("\n% Asynchrony = ", round(x$Relative["Psi_cont"], 2)),
      paste0("\n% Averaging = ", round(x$Relative["omega_cont"], 2)))
}

as.data.frame.comstab <- function(x, ...) {
  data.frame(CVc = x$CVs, 
             tau = x$Stabilization["tau"],
             delta = x$Stabilization["Delta"],
             psi = x$Stabilization["Psi"],
             omega = x$Stabilization["omega"],
             delta_rel = x$Relative["Delta_cont"],
             psi_rel = x$Relative["Psi_cont"],
             omega_rel = x$Relative["omega_cont"])
}