test_that("glm_lr_evalues works for Poisson family", {
  set.seed(123)
  n <- 50
  x <- rnorm(n)
  # Generate Poisson data
  lambda <- exp(0.5 * x)
  y <- rpois(n, lambda = lambda)
  X <- cbind(1, x)
  
  # Compute e-values
  evalues <- glm_lr_evalues(y, X, theta_null = 0, theta_alt = 0.5, 
                            family = poisson())
  
  # Check output
  expect_type(evalues, "double")
  expect_length(evalues, n)
  expect_true(all(is.finite(evalues)))
  expect_true(all(evalues > 0))
})

test_that("glm_lr_evalues works for binomial family", {
  set.seed(123)
  n <- 50
  x <- rnorm(n)
  # Generate binomial data
  p <- plogis(0.5 * x)
  y <- rbinom(n, size = 1, prob = p)
  X <- cbind(1, x)
  
  # Compute e-values
  evalues <- glm_lr_evalues(y, X, theta_null = 0, theta_alt = 0.5,
                            family = binomial())
  
  # Check output
  expect_type(evalues, "double")
  expect_length(evalues, n)
  expect_true(all(is.finite(evalues)))
  expect_true(all(evalues > 0))
})

test_that("glm_lr_evalues works for Gamma family", {
  set.seed(123)
  n <- 50
  x <- rnorm(n)
  # Generate Gamma data
  shape <- 2
  mu <- exp(0.5 * x)
  y <- rgamma(n, shape = shape, rate = shape / mu)
  X <- cbind(1, x)
  
  # Compute e-values
  evalues <- glm_lr_evalues(y, X, theta_null = 0, theta_alt = 0.5,
                            family = Gamma())
  
  # Check output
  expect_type(evalues, "double")
  expect_length(evalues, n)
  expect_true(all(is.finite(evalues)))
  expect_true(all(evalues > 0))
})

test_that("glm_lr_evalues works for negative binomial family", {
  set.seed(123)
  n <- 50
  x <- rnorm(n)
  # Generate negative binomial data
  size <- 2
  mu <- exp(0.5 * x)
  y <- rnbinom(n, size = size, mu = mu)
  X <- cbind(1, x)
  
  # Compute e-values with specified dispersion
  evalues <- glm_lr_evalues(y, X, theta_null = 0, theta_alt = 0.5,
                            family = "negative.binomial", dispersion = size)
  
  # Check output
  expect_type(evalues, "double")
  expect_length(evalues, n)
  expect_true(all(is.finite(evalues)))
  expect_true(all(evalues > 0))
})

test_that("glm_lr_evalues estimates dispersion when not provided", {
  set.seed(123)
  n <- 100
  x <- rnorm(n)
  
  # Negative binomial without specified size
  mu <- exp(0.5 * x)
  y <- rnbinom(n, size = 3, mu = mu)
  X <- cbind(1, x)
  
  # Should estimate dispersion
  evalues <- glm_lr_evalues(y, X, theta_null = 0, theta_alt = 0.5,
                            family = "negative.binomial")
  
  expect_type(evalues, "double")
  expect_length(evalues, n)
  expect_true(all(evalues > 0))
})

test_that("glm_lr_evalues gives informative error for unsupported family", {
  n <- 50
  x <- rnorm(n)
  y <- rnorm(n)
  X <- cbind(1, x)
  
  # Test with unsupported family
  expect_error(
    glm_lr_evalues(y, X, theta_null = 0, theta_alt = 0.5,
                   family = "unsupported"),
    "not yet supported"
  )
  
  # Error message should mention references
  expect_error(
    glm_lr_evalues(y, X, theta_null = 0, theta_alt = 0.5,
                   family = "quasi"),
    "Turner, Ly"
  )
})

test_that("glm_lr_evalues E[e] <= 1 under null for different families", {
  set.seed(456)
  n <- 100
  x <- rnorm(n)
  X <- cbind(1, x)
  
  # Poisson under null
  y_pois <- rpois(n, lambda = exp(0 * x + 1))  # theta = 0
  e_pois <- glm_lr_evalues(y_pois, X, theta_null = 0, theta_alt = 0.3,
                           family = poisson())
  expect_true(mean(e_pois) <= 1.5)  # Allow some sampling variation
  
  # Binomial under null
  y_binom <- rbinom(n, size = 1, prob = plogis(0 * x))  # theta = 0
  e_binom <- glm_lr_evalues(y_binom, X, theta_null = 0, theta_alt = 0.3,
                            family = binomial())
  expect_true(mean(e_binom) <= 1.5)
})
