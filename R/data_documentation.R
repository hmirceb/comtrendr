#' Simulated community data in long format
#'
#' @name example_data_long
#' @format A data frame with 300 rows and 4 columns:
#' \describe{
#'   \item{comm}{Community ID}
#'   \item{time}{Time step}
#'   \item{species}{Species names}
#'   \item{abundance}{Abundance value}
#' }
#' @source Simulated
#' @export
load("data/example_data_long.rda")
NULL

#' Simulated community data in wide format
#'
#' @name example_data_wide
#' @format A data frame with 30 rows and 12 columns:
#' \describe{
#'   \item{comm}{Community ID}
#'   \item{time}{Time step}
#'   \item{abundance.sp_X}{Abundance of species X}
#' }
#' @source Simulated
#' @export
load("data/example_data_wide.rda")
NULL