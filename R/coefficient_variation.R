#' Calculate the coefficient of variation (CV) of a variable
#'
#' @param x A numeric vector
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#'
#' @returns
#' @export
#'
cv_term <- function(x, term = "var") {
  # Match argument for variance function to use
  options <- data.frame(term = c("var", "two", "three"),
                        var = c("var", "var_t2", "var_t3"))
  # Get choice
  var_func <- match.fun(options[options$term == term,]$var)
  
  # Remove NAs
  x <- x[!is.na(x)]
  
  vari <- var_func(x) # Compute variance
  sdi <- sqrt(vari) # Transform to SD
  cvi <- sdi/mean(x) # Calculate CV
  
  names(cvi) <- paste0("CV_", term)
  return(cvi)
}

#' CV of a community
#'
#' @param x A data.frame. Community matrix with time in rows and taxa in columns.
#' @param weighted Boolean. Weight the CV of each population by its average relative abundance per species across years.
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#'
#' @returns
#' @export
#'
cv_com_term <- function(x, total = TRUE, weighted = FALSE, term = "var") {
  # Check valid data for weighted option
  if( isTRUE(weighted) & ncol(as.matrix(x)) == 1 ) {
    stop("Weights can only be applied to communities")
  }
  
  if( isTRUE(total) ) {
    com_t <- rowSums(x, na.rm = T) # Sum of species abundance per year
    cvc <- cv_term(com_t, term = term)
    names(cvc) <- paste0("CVt_", term) 
  } else {
    # Calculate weighted CV
    if( isTRUE(weighted) ) {
      ps <- colMeans(x/rowSums(x, na.rm = T), na.rm = T) # Calculate weigths (average relative abundance per species across years)
      cvc <- sum( # Weighted mean
        apply(x, 2, cv_term, term = term) * # CVs if each species
          ps # Multiply by weight
      ) # Sum
      names(cvc) <- paste0("CVw_", term) 
    } else {
      # Calculate unweighted mean CV
      cvc <- mean(
        apply(x, 2, cv_term, term = term)
      )
      names(cvc) <- paste0("CVm_", term)
    }
  }
  return(cvc)
}