#' Example Regression Dataset
#'
#' A simulated dataset for demonstrating safe Bayesian regression methods.
#' Contains a continuous response variable and three predictor variables.
#'
#' @format A data frame with 100 rows and 4 variables:
#' \describe{
#'   \item{y}{Response variable, continuous}
#'   \item{x1}{First predictor variable, standard normal}
#'   \item{x2}{Second predictor variable, standard normal}
#'   \item{x3}{Third predictor variable, standard normal}
#' }
#'
#' @details
#' The data was generated as:
#' \itemize{
#'   \item y = 0.5 * x1 + 0.3 * x2 + epsilon
#'   \item x1, x2, x3 ~ N(0, 1)
#'   \item epsilon ~ N(0, 1)
#' }
#' Thus x1 and x2 have true effects while x3 has no effect.
#'
#' @examples
#' data(example_data)
#' head(example_data)
#' 
#' # Fit a model
#' \dontrun{
#' fit <- fit_safe_brms(y ~ x1 + x2 + x3, data = example_data)
#' }
"example_data"
