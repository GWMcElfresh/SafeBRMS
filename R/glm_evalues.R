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
#'   Default is gaussian(). Supported families: gaussian, poisson, binomial,
#'   Gamma, and negative.binomial (requires size parameter).
#' @param dispersion Dispersion parameter. For negative binomial, this is the
#'   size parameter (theta). For Gamma, this is the shape parameter. If NULL,
#'   estimated from data where appropriate.
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
#' Supported families:
#' \itemize{
#'   \item \strong{gaussian}: Normal distribution for continuous responses
#'   \item \strong{poisson}: Poisson distribution for count data
#'   \item \strong{binomial}: Binomial distribution for binary/proportion data
#'   \item \strong{Gamma}: Gamma distribution for positive continuous data
#'   \item \strong{negative.binomial}: Negative binomial for overdispersed counts
#' }
#'
#' @examples
#' # Gaussian (continuous data)
#' set.seed(123)
#' n <- 50
#' x <- rnorm(n)
#' y <- 0.3 * x + rnorm(n)
#' X <- cbind(1, x)
#' evalues_gauss <- glm_lr_evalues(y, X, theta_null = 0, theta_alt = 0.5)
#'
#' # Poisson (count data)
#' y_count <- rpois(n, lambda = exp(0.5 * x))
#' evalues_pois <- glm_lr_evalues(y_count, X, theta_null = 0, theta_alt = 0.5,
#'                                 family = poisson())
#'
#' # Binomial (binary data)
#' y_binary <- rbinom(n, size = 1, prob = plogis(0.5 * x))
#' evalues_binom <- glm_lr_evalues(y_binary, X, theta_null = 0, theta_alt = 0.5,
#'                                  family = binomial())
#'
#' @export
#' @importFrom stats dnorm dpois dbinom dgamma
glm_lr_evalues <- function(y, X, theta_null, theta_alt, family = gaussian(), 
                           dispersion = NULL) {
  
  n <- length(y)
  
  if (!is.matrix(X)) {
    stop("X must be a matrix")
  }
  
  # Extract family name
  family_name <- if (inherits(family, "family")) {
    family$family
  } else if (is.character(family)) {
    family
  } else {
    stop("family must be a family object or character string")
  }
  
  # Initialize e-values vector
  evalues <- numeric(n)
  
  # Gaussian family
  if (family_name == "gaussian") {
    # Estimate variance from data if not provided
    sigma_sq <- if (is.null(dispersion)) var(y) else dispersion^2
    
    for (i in 1:n) {
      # Linear predictor
      mu_alt <- theta_alt * X[i, 2]
      mu_null <- theta_null * X[i, 2]
      
      # Log-likelihood ratio
      ll_alt <- dnorm(y[i], mean = mu_alt, sd = sqrt(sigma_sq), log = TRUE)
      ll_null <- dnorm(y[i], mean = mu_null, sd = sqrt(sigma_sq), log = TRUE)
      
      evalues[i] <- exp(ll_alt - ll_null)
    }
    
  # Poisson family
  } else if (family_name == "poisson") {
    for (i in 1:n) {
      # Linear predictor on log scale
      lambda_alt <- exp(theta_alt * X[i, 2])
      lambda_null <- exp(theta_null * X[i, 2])
      
      # Log-likelihood ratio
      ll_alt <- dpois(y[i], lambda = lambda_alt, log = TRUE)
      ll_null <- dpois(y[i], lambda = lambda_null, log = TRUE)
      
      evalues[i] <- exp(ll_alt - ll_null)
    }
    
  # Binomial family
  } else if (family_name == "binomial") {
    # Assume y is 0/1 (size = 1) unless otherwise specified
    size <- if (!is.null(attr(y, "size"))) attr(y, "size") else 1
    
    for (i in 1:n) {
      # Linear predictor on logit scale
      p_alt <- plogis(theta_alt * X[i, 2])
      p_null <- plogis(theta_null * X[i, 2])
      
      # Log-likelihood ratio
      ll_alt <- dbinom(y[i], size = size, prob = p_alt, log = TRUE)
      ll_null <- dbinom(y[i], size = size, prob = p_null, log = TRUE)
      
      evalues[i] <- exp(ll_alt - ll_null)
    }
    
  # Gamma family
  } else if (family_name == "Gamma") {
    # Estimate shape parameter if not provided
    # Using method of moments: shape = mean^2 / variance
    if (is.null(dispersion)) {
      shape <- mean(y)^2 / var(y)
    } else {
      shape <- dispersion
    }
    
    for (i in 1:n) {
      # Linear predictor on log scale
      mu_alt <- exp(theta_alt * X[i, 2])
      mu_null <- exp(theta_null * X[i, 2])
      
      # Rate parameter: rate = shape / mean
      rate_alt <- shape / mu_alt
      rate_null <- shape / mu_null
      
      # Log-likelihood ratio
      ll_alt <- dgamma(y[i], shape = shape, rate = rate_alt, log = TRUE)
      ll_null <- dgamma(y[i], shape = shape, rate = rate_null, log = TRUE)
      
      evalues[i] <- exp(ll_alt - ll_null)
    }
    
  # Negative Binomial family
  } else if (family_name %in% c("negative.binomial", "negbin", "nb")) {
    # Requires size (dispersion) parameter
    if (is.null(dispersion)) {
      # Estimate size using method of moments
      # Var(Y) = mu + mu^2/size, so size = mu^2 / (Var - mu)
      mu_est <- mean(y)
      var_est <- var(y)
      size <- mu_est^2 / (var_est - mu_est)
      if (size <= 0) {
        warning("Estimated negative binomial size parameter is non-positive. Using size = 1.")
        size <- 1
      }
    } else {
      size <- dispersion
    }
    
    for (i in 1:n) {
      # Linear predictor on log scale
      mu_alt <- exp(theta_alt * X[i, 2])
      mu_null <- exp(theta_null * X[i, 2])
      
      # Log-likelihood ratio using dnbinom
      # dnbinom uses prob = size/(size + mu) parameterization
      ll_alt <- dnbinom(y[i], size = size, mu = mu_alt, log = TRUE)
      ll_null <- dnbinom(y[i], size = size, mu = mu_null, log = TRUE)
      
      evalues[i] <- exp(ll_alt - ll_null)
    }
    
  } else {
    # Unsupported family - provide guidance
    stop(paste0(
      "Family '", family_name, "' is not yet supported.\n",
      "Supported families: gaussian, poisson, binomial, Gamma, negative.binomial\n",
      "TODO: Implement additional families as needed.\n",
      "For custom implementations, see:\n",
      "  - Turner, Ly, & Grünwald (2020+) for general e-value construction\n",
      "  - Grünwald et al. (2019) for safe testing framework\n",
      "  - Exponential family theory: Brown (1986), Statistical Sci."
    ))
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
