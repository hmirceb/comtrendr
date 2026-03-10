#' Calculate the coefficient of variation (CV) of a variable
#'
#' @param x A numeric vector
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#'
#' @returns A numeric value.
#' 
#' @examples
#' require(detrending)
#' d <- rnorm(1000, mean = 10, sd = 30)
#' cv_term(d) # ~3
#' 
#' @export
cv_term <- function(x, term = "var") {
  # Match variance function
  var_func <- switch(
    term,
    var = stats::var,
    two = var_t2,
    three = var_t3
  )
  
  # Remove NAs
  x <- x[!is.na(x)]
  
  vari <- var_func(x) # Compute variance with desired function
  sdi <- sqrt(vari) # Transform to SD
  cvi <- sdi/mean(x) # Calculate CV
  
  names(cvi) <- paste0("CV_", term)
  return(cvi)
}

#' CV of a community
#'
#' @param x A data.frame. Community matrix with time in rows and taxa in columns.
#' @param total Boolean. If TRUE, compute CV of the sum of annual abundances. If FALSE computes the average of the CV of each species. Default TRUE. 
#' @param weighted Boolean. Weight the CV of each population by its average relative abundance per species across years. Default FALSE.
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#' @param time_col Character. Name of the column with time variable. Optional with default "time".
#'
#' @returns If total = TRUE, a numeric value with the CV of the sum of annual abundances. If total = FALSE, a named list with the average CV across the populations in the community and the CVs of each population.
#' 
#' @examples
#' require(detrending)
#' 
#' # Load and clean data
#' data(example_data_wide)
#' metacomm_df <- clean_community_wide(x = example_data_wide)
#' comm_df <- metacomm_df[metacomm_df$comm == 1,][,-c(1:2)]
#' 
#' # Calculate CV of total abundance
#' cv_com_term(x = comm_df, 
#'            total = TRUE,
#'            term = "var")
#'            
#' # Calculate average CV of populations
#' cv_com_term(x = comm_df, 
#'            total = FALSE,
#'            term = "var")
#'            
#' @export
cv_com_term <- function(x, total = TRUE, weighted = FALSE, term = "var", time_col = "time") {
  
  # Check if a time column was specified for detrending methods and order rows
  x <- check_time(x, time_col = time_col, term = term, rm = TRUE)
  
  # Replace NAs with 0 and remove columns (species) with 0 abundance across all years 
  x <- remove_empty_sps(x = x, time_col = time_col)
  
  # Check valid data for weighted option
  if( isTRUE(weighted) & ncol(as.matrix(x)) == 1 ) {
    stop("Weights cannot be applied to a single species")
  }
  
  if( isTRUE(total) ) {
    # Sum of species abundance per year
    com_t <- rowSums(x, na.rm = T) 
    cvc <- cv_term(com_t, term = term)
    names(cvc) <- paste0("CVt_", term) 
    
  } else {
    # CVs of each species
    cv_sps <- apply(x, 2, cv_term, term = term)
    
    # Calculate weighted CV
    if( isTRUE(weighted) ) {
      # Calculate weigths (average relative abundance per species across years)
      ps <- colMeans(x/rowSums(x, na.rm = T), na.rm = T)
      cv_p <- sum( # Weighted mean
        cv_sps * ps # Multiply by weight
      )
      names(cv_p) <- paste0("CV_popw_", term) 
      
    } else {
      # Calculate unweighted mean CV
      cv_p <- mean(
        cv_sps
      )
      names(cv_p) <- paste0("CV_pop_", term)
    }
    cvc <- list(cv_p, cv_sps)
    names(cvc) <- c("CV", "CV_species")
  }
  
  return(cvc)
}