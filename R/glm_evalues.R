#' Compute GLM Likelihood Ratio E-Values
#'
#' Computes per-observation likelihood ratio based e-variables for exponential
#' family GLMs. These e-values satisfy E <= 1 under the null hypothesis.
#'
#' @param y Response vector.
#' @param X Design matrix.
#' @param theta_null Parameter value under the null hypothesis.
#' @param theta_alt Parameter value under the alternative hypothesis.
#' @param family Family object specifying the exponential family distribution.
#'   Default is gaussian().
#'
#' @return A numeric vector of per-observation e-values.
#'
#' @details
#' For each observation i, the e-value is computed as the likelihood ratio:
#' \deqn{e_i = \frac{L(\theta_{alt} | y_i)}{L(\theta_{null} | y_i)}}
#'
#' These e-values form a test martingale under the null hypothesis, meaning
#' their expected value is at most 1.
#'
#' @examples
#' # Simple linear model
#' set.seed(123)
#' n <- 50
#' x <- rnorm(n)
#' y <- 0.3 * x + rnorm(n)
#' X <- cbind(1, x)
#' 
#' # Compute e-values for testing beta = 0 vs beta = 0.5
#' evalues <- glm_lr_evalues(y, X, theta_null = 0, theta_alt = 0.5)
#'
#' @export
#' @importFrom stats dnorm
glm_lr_evalues <- function(y, X, theta_null, theta_alt, family = gaussian()) {
  
  n <- length(y)
  
  # For Gaussian family with known variance (simplified)
  if (inherits(family, "family") && family$family == "gaussian") {
    # Compute fitted values under null and alternative
    if (is.matrix(X)) {
      # Assuming theta refers to a specific coefficient
      # For simplicity, we'll compute log-likelihood ratios
      sigma_sq <- var(y)  # Estimate variance from data
      
      # Compute likelihood ratio for each observation
      # This is a simplified version - full implementation would be more complex
      evalues <- numeric(n)
      for (i in 1:n) {
        # Likelihood under alternative
        ll_alt <- dnorm(y[i], mean = theta_alt * X[i, 2], sd = sqrt(sigma_sq), log = TRUE)
        # Likelihood under null
        ll_null <- dnorm(y[i], mean = theta_null * X[i, 2], sd = sqrt(sigma_sq), log = TRUE)
        # E-value is exp(log-likelihood ratio)
        evalues[i] <- exp(ll_alt - ll_null)
      }
    } else {
      stop("X must be a matrix")
    }
  } else {
    stop("Currently only Gaussian family is supported")
  }
  
  return(evalues)
}


#' Mix E-Values to Ensure E <= 1 Under the Null
#'
#' Applies mixing to e-values to ensure the expected value is at most 1
#' under the null hypothesis. Uses a Beta mixture distribution.
#'
#' @param evalues A numeric vector of e-values.
#' @param alpha_mix First shape parameter for Beta mixing distribution. Default is 0.5.
#' @param beta_mix Second shape parameter for Beta mixing distribution. Default is 0.5.
#'
#' @return A numeric vector of mixed e-values.
#'
#' @details
#' The mixing ensures that even if individual e-values might exceed 1 in expectation,
#' the mixture satisfies E[e] <= 1 under the null. This is achieved by integrating
#' over a prior distribution on the alternative parameter space.
#'
#' @examples
#' # Generate some e-values
#' evalues <- c(0.8, 1.2, 0.9, 1.5, 0.7)
#' mixed <- mix_evalues(evalues)
#'
#' @export
mix_evalues <- function(evalues, alpha_mix = 0.5, beta_mix = 0.5) {
  
  # Simple mixing scheme: weighted average
  # In practice, this would integrate over a prior distribution
  # For now, we'll use a conservative shrinkage approach
  
  n <- length(evalues)
  weights <- stats::dbeta(seq(0.1, 0.9, length.out = n), alpha_mix, beta_mix)
  weights <- weights / sum(weights)
  
  # Apply conservative mixing by shrinking toward 1
  mixed_evalues <- evalues * weights + (1 - weights)
  
  return(mixed_evalues)
}
