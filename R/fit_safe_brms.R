#' Fit a Safe BRMS Model
#'
#' Fits a Bayesian regression model using brms with safe testing in mind.
#' This is a wrapper around brms::brm that ensures proper setup for subsequent
#' e-value calculations.
#'
#' @param formula A formula describing the model to be fitted in brms syntax.
#' @param data A data frame containing the variables in the model.
#' @param family A description of the response distribution and link function.
#'   Can be a family function, a call to a family function or a character string
#'   naming the family. Default is gaussian().
#' @param prior Prior specifications for model parameters. Default uses brms default priors.
#' @param chains Number of Markov chains. Default is 4.
#' @param iter Number of total iterations per chain. Default is 2000.
#' @param warmup Number of warmup iterations per chain. Default is 1000.
#' @param ... Additional arguments passed to brms::brm.
#'
#' @return A brmsfit object containing the fitted model.
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' n <- 100
#' x <- rnorm(n)
#' y <- 0.5 * x + rnorm(n)
#' data <- data.frame(y = y, x = x)
#' 
#' # Fit model
#' fit <- fit_safe_brms(y ~ x, data = data)
#' }
#'
#' @export
#' @importFrom methods is
fit_safe_brms <- function(formula, data, family = gaussian(),
                          prior = NULL, chains = 4, iter = 2000,
                          warmup = 1000, ...) {
  
  # Check if brms is available
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required but not installed. Please install it with install.packages('brms')")
  }
  
  # Fit the brms model
  fit <- brms::brm(
    formula = formula,
    data = data,
    family = family,
    prior = prior,
    chains = chains,
    iter = iter,
    warmup = warmup,
    ...
  )
  
  # Store original data for later use
  attr(fit, "safe_data") <- data
  
  return(fit)
}
