# Functions to calculate correlation based on detrended variance

#' Pearson's correlation index between two variables based on variance or Hill's 2 and 3 term variance
#'
#' @param x,y A pair of numeric vectors.
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariane), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#'
#' @returns
#' 
#' @export
cor_term <- function(x, y, term = "var") {
  # Match argument for variance function to use
  options <- data.frame(term = c("var", "two", "three", "linear"),
                        var = c("var", "var_t2", "var_t3", "var_linear"))
  # Get choice
  var_func <- match.fun(options[options$term == term,]$var)
  
  cov <- cov_term(x, y, term = term) # Covariance between both variables
  sd1 = sqrt(var_func(x)) # SD of x
  sd2 = sqrt(var_func(y)) # SD of y
  
  cor <- cov/(sd1*sd2) # Correlation
  return(cor)
}

#' Pearson's correlation matrix
#'
#' @param x A dataframe. A community matrix of species abundance with years as rows and species as columns.
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariane), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#'
#' @returns A symmetric correlation matrix.
#' 
#' @export
cormat_term = function(x, term = "var") {
  # Get number of species
  nc <- ncol(x)
  # Make empty correlation matrix
  S <- matrix(NA, nc, nc)
  
  # If there are only two species return only one correlation value
  if (nc == 2) {
    S <- cor_term(x = x[,1], y = x[,2], term = term)
    
  } else {
    # Populate matrix
    for (i in 1:nc) {
      # Diagonal
      S[i, i] <- cor_term(x = x[,i], y = x[,i], term = term)
      
      # Off diagonals
      for (j in (i + 1):nc) {
        if (i < nc) {
          S[i, j] <- cor_term(x = x[,i], y = x[,j], term = term)
          S[j, i] <- S[i, j] # make it symmetric
        }
      }
    }
  }
  return(S)
}