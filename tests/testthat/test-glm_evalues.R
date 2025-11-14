test_that("glm_lr_evalues computes e-values correctly", {
  # Simple test data
  set.seed(123)
  n <- 50
  x <- rnorm(n)
  y <- 0.3 * x + rnorm(n)
  X <- cbind(1, x)
  
  # Compute e-values
  evalues <- glm_lr_evalues(y, X, theta_null = 0, theta_alt = 0.5)
  
  # Check output
  expect_type(evalues, "double")
  expect_length(evalues, n)
  expect_true(all(is.finite(evalues)))
  expect_true(all(evalues > 0))
})

test_that("mix_evalues produces valid output", {
  evalues <- c(0.8, 1.2, 0.9, 1.5, 0.7)
  mixed <- mix_evalues(evalues)
  
  expect_type(mixed, "double")
  expect_length(mixed, length(evalues))
  expect_true(all(is.finite(mixed)))
  expect_true(all(mixed > 0))
})

test_that("mix_evalues with different parameters", {
  evalues <- runif(20, 0.5, 2)
  mixed1 <- mix_evalues(evalues, alpha_mix = 0.5, beta_mix = 0.5)
  mixed2 <- mix_evalues(evalues, alpha_mix = 1, beta_mix = 1)
  
  expect_type(mixed1, "double")
  expect_type(mixed2, "double")
  expect_false(identical(mixed1, mixed2))
})
