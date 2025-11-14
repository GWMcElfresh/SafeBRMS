#' Multiple Testing with Benjamini-Hochberg Correction
#'
#' Performs multiple hypothesis tests with Benjamini-Hochberg FDR control
#' using e-values instead of p-values.
#'
#' @param fit A brmsfit object from fit_safe_brms or brms::brm.
#' @param coef_names Vector of coefficient names to test. If NULL, tests all coefficients.
#' @param null_value Value under the null hypothesis for all coefficients. Default is 0.
#' @param fdr_level Target false discovery rate. Default is 0.05.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item coefficient: Name of each coefficient
#'   \item evalue: E-value for each test
#'   \item adjusted_threshold: BH-adjusted threshold for each test
#'   \item reject: Whether to reject the null hypothesis
#'   \item estimate: Posterior mean of coefficient
#' }
#'
#' @details
#' Applies the Benjamini-Hochberg procedure to control the false discovery rate
#' when testing multiple coefficients. E-values are converted to pseudo-p-values
#' for the BH procedure, then decisions are made based on adjusted thresholds.
#'
#' @examples
#' \dontrun{
#' # Fit model with multiple predictors
#' data <- data.frame(
#'   y = rnorm(100),
#'   x1 = rnorm(100),
#'   x2 = rnorm(100),
#'   x3 = rnorm(100)
#' )
#' fit <- fit_safe_brms(y ~ x1 + x2 + x3, data = data)
#' 
#' # Multiple testing with BH correction
#' results <- multiple_testing_bh(fit)
#' print(results)
#' }
#'
#' @export
#' @importFrom stats p.adjust
multiple_testing_bh <- function(fit, coef_names = NULL, null_value = 0, fdr_level = 0.05) {
  
  # Compute e-values for all coefficients
  evals <- compute_evalues(fit, coef_name = coef_names, null_value = null_value)
  
  # Get coefficient estimates
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required")
  }
  
  post_samples <- brms::as_draws_df(fit)
  
  n_tests <- length(evals$evalues)
  
  # Convert e-values to pseudo p-values
  # Using the relationship: p = 1 / evalue (conservative)
  pseudo_pvalues <- pmin(1, 1 / evals$evalues)
  
  # Apply Benjamini-Hochberg procedure
  adjusted_pvals <- stats::p.adjust(pseudo_pvalues, method = "BH")
  
  # Make decisions
  reject <- adjusted_pvals < fdr_level
  
  # Get coefficient estimates
  estimates <- numeric(n_tests)
  for (i in seq_along(evals$coef_names)) {
    coef_col <- paste0("b_", evals$coef_names[i])
    if (coef_col %in% names(post_samples)) {
      estimates[i] <- mean(post_samples[[coef_col]])
    }
  }
  
  # Create results data frame
  results <- data.frame(
    coefficient = evals$coef_names,
    evalue = evals$evalues,
    pseudo_pvalue = pseudo_pvalues,
    adjusted_pvalue = adjusted_pvals,
    reject = reject,
    estimate = estimates,
    stringsAsFactors = FALSE
  )
  
  # Sort by adjusted p-value
  results <- results[order(results$adjusted_pvalue), ]
  
  class(results) <- c("safe_multiple_test", "data.frame")
  
  return(results)
}


#' Print Method for Multiple Testing Results
#'
#' @param x A safe_multiple_test object
#' @param ... Additional arguments (not used)
#'
#' @export
print.safe_multiple_test <- function(x, ...) {
  cat("Safe Multiple Testing with Benjamini-Hochberg Correction\n")
  cat("=========================================================\n\n")
  
  # Print summary
  n_total <- nrow(x)
  n_reject <- sum(x$reject)
  cat(sprintf("Total tests: %d\n", n_total))
  cat(sprintf("Rejected: %d\n", n_reject))
  cat(sprintf("FDR level: %.3f\n\n", max(x$adjusted_pvalue[x$reject], 0)))
  
  # Print results table
  print(as.data.frame(x), row.names = FALSE)
}
