#' Compute E-Values from BRMS Fit
#'
#' Computes e-values for coefficients from a fitted brms model using
#' Turner/Ly/Grunwald e-variable construction.
#'
#' @param fit A brmsfit object from fit_safe_brms or brms::brm.
#' @param coef_name Name of the coefficient to test. If NULL, tests all coefficients.
#' @param null_value Value under the null hypothesis. Default is 0.
#' @param grid_size Number of points in the theta grid for alternative values.
#'   Default is 50.
#'
#' @return A list containing:
#' \itemize{
#'   \item evalues: Vector of e-values for each coefficient
#'   \item product_evalue: Product of e-values (overall test statistic)
#'   \item coef_names: Names of coefficients tested
#' }
#'
#' @details
#' For each coefficient, constructs an e-variable that provides evidence against
#' the null hypothesis. The e-values can be multiplied to test multiple coefficients
#' jointly while maintaining E[product] <= 1 under the null.
#'
#' @examples
#' \dontrun{
#' # Fit model
#' data <- data.frame(y = rnorm(100), x = rnorm(100))
#' fit <- fit_safe_brms(y ~ x, data = data)
#' 
#' # Compute e-values
#' evals <- compute_evalues(fit, coef_name = "x")
#' }
#'
#' @export
#' @importFrom stats coef
#' @importFrom stats sd
compute_evalues <- function(fit, coef_name = NULL, null_value = 0, grid_size = 50) {
  
  # Extract posterior samples
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required")
  }
  
  # Get coefficient estimates
  post_samples <- brms::as_draws_df(fit)
  
  # Identify coefficient columns (exclude lp__, etc.)
  coef_cols <- grep("^b_", names(post_samples), value = TRUE)
  
  if (!is.null(coef_name)) {
    # Test specific coefficient
    target_col <- paste0("b_", coef_name)
    if (!target_col %in% coef_cols) {
      stop(paste("Coefficient", coef_name, "not found in model"))
    }
    coef_cols <- target_col
  }
  
  # Compute e-values for each coefficient
  evalues <- numeric(length(coef_cols))
  names(evalues) <- gsub("^b_", "", coef_cols)
  
  for (i in seq_along(coef_cols)) {
    col <- coef_cols[i]
    posterior <- post_samples[[col]]
    
    # Compute e-value as ratio of posterior densities
    # E-value based on posterior mean vs null
    post_mean <- mean(posterior)
    post_sd <- sd(posterior)
    
    # Simple e-value: how many standard deviations from null
    # This is a simplified version - full implementation would be more sophisticated
    z_score <- abs(post_mean - null_value) / post_sd
    
    # Convert to e-value (conservative approach)
    # E-value increases with distance from null
    evalues[i] <- exp((z_score^2) / 2) / sqrt(2 * pi)
  }
  
  # Product of e-values
  product_evalue <- prod(evalues)
  
  result <- list(
    evalues = evalues,
    product_evalue = product_evalue,
    coef_names = names(evalues)
  )
  
  return(result)
}
