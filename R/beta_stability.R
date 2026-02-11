#' Chord transformation for a matrix of abundances
#'
#' Applies chord transformation to a matrix (see Details).
#'
#' @param x A community abundance matrix.
#'
#' @details
#' Chord transformation standardizes values by dividing values in each row by the sample norm (*SN*) of the corresponding row following: 
#' \deqn{SN = \sqrt{\sum_{i=1}^{S} x_i^2}}
#'
#' @returns A matrix of community abundance after chord transformation.
#'
#' @references 
#' - Orlóci, L. (1967). An agglomerative method for classification of plant communities. The Journal of Ecology, 193-206.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' @export
chord_transform <- function(x) {
  c_t <- x / sqrt(rowSums(x^2))
  return(c_t)
}

#' Multivariate variance of community composition
#'
#' @param x A community abundance matrix.
#' @param d Character. Distance metric to use. One of 'euclidean' or 'chord'. Default 'euclidean'.
#'
#' @details
#' The multivariate variance of community composition is defined as the average square Euclidean distance between annual observations and the average community composition following: 
#' \deqn{var_{mv} = \dfrac{ \sum_{i=1}^{t}{ED(X_{i}, X_{mean})^2} }{t-1}}
#' Where \eqn{ED(X_{i}, X_{mean})} is the Euclidean distance between the composition of the community (\eqn{X_{i}}) at time \eqn{i}.
#'
#' @returns
#' @author Jan Lepš, \email{}
#' @author Aleš Lisner, \email{lisnea00@@jcu.cz}
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' @export
#'
var_mv <- function(x, d = c("euclidean", "chord")){
  if( !d %in% c("euclidean", "chord") ){
    stop("Unsuitable distance metric. Please choose one of 'euclidean' or 'chord'")
  }
  
  # Create DF with average abundance values per species in first row
  mean_vec <- colMeans(x)
  dd <- rbind(mean_vec, x)
  
  # Apply chord transformation if necessary
  if( d == "chord" ) {
    dd <- chord_transform(dd)
  }
  # Compute distances
  dis <- as.matrix(dist(dd))
  # First column is distances between average community and each year (remove first obs)
  sum_sq <- sum(dis[-1,1]^2)
  mv_var <- sum_sq / (nrow(x)-1)
  
  return(mv_var)
}

#' Multivariate two term local quadratic variance of community composition (TTLQV_2)
#'
#' @param x A community abundance matrix.
#' @param d Character. Distance metric to use. One of 'euclidean' or 'chord'. Default 'euclidean'.
#'
#' @details
#' The multivariate two term local quadratic variance (\eqn{TTLQV_{mv2}}) of community composition is the detrended version of multivariate variance (\eqn{var_{mv}}): 
#' \deqn{TTLQV_{mv2} = \dfrac{ \sum_{i=1}^{t-1}{ED(X_{i} - X_{i+1})^2} }{2(t-1)}}
#' Where \eqn{X_{i}} is the composition of the community at time \eqn{i}.
#'
#' @returns
#' @export
#'
#' @author Aleš Lisner, \email{lisnea00@@jcu.cz}
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
var_t2mv <- function(x, d = c("euclidean", "chord")){
  if( !d %in% c("euclidean", "chord") ){
    stop("Unsuitable distance metric. Please choose one of 'euclidean' or 'chord'")
  }
  
  # Apply chord transformation if necessary
  if( d == "chord" ) {
    x <- chord_transform(x)
  }
  # Compute distances
  dis <- as.matrix(dist(x))
  # Get superdiagonal (distance between consecutive years)
  dis <- dis[row(dis) == col(dis) + 1]
  mv_var <- var_t2(dis)
  
  return(mv_var)
}

#' Multivariate three term local quadratic variance of community composition
#'
#' @param x A community abundance matrix.
#' @param d Character. Distance metric to use. One of 'euclidean' or 'chord'. Default 'euclidean'.
#'
#' @returns
#'
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#'
#' @export
#'
var_t3mv <- function(x, method = c("euclidean", "chord")){
  if( !d %in% c("euclidean", "chord") ){
    stop("Unsuitable distance metric. Please choose one of 'euclidean' or 'chord'")
  }
  
  # Apply chord transformation if necessary
  if( method == "chord" ) {
    x <- chord_transform(x)
  }
  # Compute distances
  dis <- as.matrix(dist(x))
  # Get superdiagonal (distance between consecutive years)
  dis <- dis[row(dis) == col(dis) + 1]
  mv_var <- var_t3(dis)
  
  return(mv_var)
}

#' Multivariate coefficient of variation
#'
#' @param x A community abundance matrix.
#' @param d Character. Distance metric to use. One of 'euclidean' or 'chord'. Default 'euclidean'.
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance) "two" for Hills' two term local quadrat variance and covariance. Default "var".
#'
#' @returns
#' @export
#'
#' @examples
cv_mv <- function(x, d = "euclidean", term = "var") {
  if( !d %in% c("euclidean", "chord") ){
    stop("Unsuitable distance metric. Please choose one of 'euclidean' or 'chord'")
  }
  if( !term %in% c("var", "two") ){
    stop("Unsuitable variance term. Please choose one of 'var' or 'two'")
  }
  
  # Match argument for variance function to use
  options <- data.frame(term = c("var", "two"),
                        var = c("var_mv", "var_t2mv"))
  # Get choice
  var_func <- match.fun(options[options$term == term,]$var)
  
  vv <- var_func(x, d = d)
  mean_vec <- colMeans(x)
  
  # Sample norm of average composition
  SN <- sqrt(sum(mean_vec^2))
  # CVmv
  if (d == "chord") {
    CVmv <- sqrt(vv)
  } else {
    CVmv <- sqrt(vv) / SN
  }
  return(CVmv)
}