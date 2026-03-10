#' Simulated community data in long format
#'
#' @format A data frame with 300 rows and 4 columns:
#' \describe{
#'   \item{comm}{Community ID}
#'   \item{time}{Time step}
#'   \item{species}{Species names}
#'   \item{abundance}{Abundance value}
#' }
#' @source Simulated
"example_data_long"

#' Simulated community data in wide format
#'
#' @format A data frame with 30 rows and 12 columns:
#' \describe{
#'   \item{comm}{Community ID}
#'   \item{time}{Time step}
#'   \item{abundance.sp_1, abundance.sp_2, abundance.sp_3,
#'   abundance.sp_4, abundance.sp_5, abundance.sp_6, 
#'   abundance.sp_7, abundance.sp_8, abundance.sp_9, 
#'   abundance.sp_10}{Abundance of species X}
#' }
#' @source Simulated
"example_data_wide"
