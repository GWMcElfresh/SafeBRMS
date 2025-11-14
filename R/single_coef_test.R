#' Single Coefficient Test
#'
#' Performs a safe test for a single coefficient being equal to a null value.
#' Returns an e-value that can be interpreted as evidence against the null.
#'
#' @param fit A brmsfit object from fit_safe_brms or brms::brm.
#' @param coef_name Name of the coefficient to test.
#' @param null_value Value under the null hypothesis. Default is 0.
#' @param alpha Significance level for reporting. Default is 0.05.
#'
#' @return A list containing:
#' \itemize{
#'   \item evalue: The e-value for the test
#'   \item reject: Logical indicating if null is rejected at level alpha
#'   \item coef_estimate: Posterior mean of the coefficient
#'   \item coef_sd: Posterior standard deviation
#' }
#'
#' @details
#' The test rejects the null hypothesis when evalue > 1/alpha. For example,
#' at alpha = 0.05, we reject when evalue > 20.
#'
#' @examples
#' \dontrun{
#' # Fit model
#' data <- data.frame(y = rnorm(100), x = rnorm(100))
#' fit <- fit_safe_brms(y ~ x, data = data)
#' 
#' # Test single coefficient
#' result <- single_coef_test(fit, "x")
#' print(result)
#' }
#'
#' @export
single_coef_test <- function(fit, coef_name, null_value = 0, alpha = 0.05) {
  
  # Compute e-values
  evals <- compute_evalues(fit, coef_name = coef_name, null_value = null_value)
  
  # Extract the e-value for this coefficient
  evalue <- evals$evalues[1]
  
  # Critical value for rejection
  critical_value <- 1 / alpha
  
  # Check if we reject
  reject <- evalue > critical_value
  
  # Get posterior statistics
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required")
  }
  
  post_samples <- brms::as_draws_df(fit)
  target_col <- paste0("b_", coef_name)
  posterior <- post_samples[[target_col]]
  
  result <- list(
    evalue = evalue,
    reject = reject,
    coef_estimate = mean(posterior),
    coef_sd = stats::sd(posterior),
    null_value = null_value,
    alpha = alpha
  )
  
  class(result) <- "safe_test_result"
  
  return(result)
}


#' Print Method for Safe Test Results
#'
#' @param x A safe_test_result object
#' @param ... Additional arguments (not used)
#'
#' @export
print.safe_test_result <- function(x, ...) {
  cat("Safe Hypothesis Test\n")
  cat("--------------------\n")
  cat(sprintf("Null value: %.3f\n", x$null_value))
  cat(sprintf("Coefficient estimate: %.3f (SD = %.3f)\n", 
              x$coef_estimate, x$coef_sd))
  cat(sprintf("E-value: %.3f\n", x$evalue))
  cat(sprintf("Significance level: %.3f\n", x$alpha))
  cat(sprintf("Decision: %s\n", 
              ifelse(x$reject, "REJECT null", "FAIL TO REJECT null")))
}
