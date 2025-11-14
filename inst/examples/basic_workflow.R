## Example: Single Coefficient Test and Multiple Testing with SafeBRMS
## This script demonstrates the core functionality of the SafeBRMS package

# Note: This example requires brms to be installed
# install.packages("brms")

# If SafeBRMS is installed as a package:
# library(SafeBRMS)

# For testing, source the R files directly:
source("R/fit_safe_brms.R")
source("R/glm_evalues.R")
source("R/compute_evalues.R")
source("R/compute_eposterior.R")
source("R/single_coef_test.R")
source("R/multiple_testing.R")

# Set seed for reproducibility
set.seed(42)

# Load example data
data(example_data)
cat("Example data loaded. Summary:\n")
print(summary(example_data))

# The data has true effects: y = 0.5*x1 + 0.3*x2 + epsilon
# x3 has no effect (null is true)

cat("\n=== Example 1: Fit Bayesian Regression Model ===\n")
# Note: This requires brms and rstan to be installed
# Uncomment the following lines to run with brms:

# fit <- fit_safe_brms(
#   formula = y ~ x1 + x2 + x3,
#   data = example_data,
#   chains = 2,        # Use 2 chains for speed
#   iter = 1000,       # 1000 iterations
#   warmup = 500,      # 500 warmup iterations
#   refresh = 0        # Suppress Stan output
# )
# 
# # View model summary
# summary(fit)

cat("\n=== Example 2: Single Coefficient Test ===\n")
# Test whether x1 coefficient is zero
# result_x1 <- single_coef_test(fit, coef_name = "x1", null_value = 0, alpha = 0.05)
# print(result_x1)
# 
# # Expected: Should reject null (x1 has true effect of 0.5)
# 
# # Test whether x3 coefficient is zero
# result_x3 <- single_coef_test(fit, coef_name = "x3", null_value = 0, alpha = 0.05)
# print(result_x3)
# 
# # Expected: Should NOT reject null (x3 has no true effect)

cat("\n=== Example 3: Multiple Testing with BH Correction ===\n")
# Test all coefficients simultaneously with FDR control
# results <- multiple_testing_bh(
#   fit = fit,
#   null_value = 0,
#   fdr_level = 0.05
# )
# 
# print(results)
# 
# # Expected results:
# # - x1 and x2 should be significant (true effects)
# # - x3 should not be significant (no true effect)

cat("\n=== Example 4: E-Posterior Computation ===\n")
# Compute and plot e-posterior for x1
# epost <- compute_eposterior(fit, coef_name = "x1")
# 
# plot(epost$theta, epost$eposterior, type = "l",
#      xlab = "Theta", ylab = "E-Posterior Density",
#      main = "E-Posterior for x1 Coefficient",
#      lwd = 2, col = "blue")
# abline(v = epost$MAP_theta, col = "red", lty = 2, lwd = 2)
# abline(v = 0.5, col = "green", lty = 2, lwd = 2)  # True value
# legend("topright", 
#        legend = c("E-Posterior", "MAP Estimate", "True Value"),
#        col = c("blue", "red", "green"), 
#        lty = c(1, 2, 2),
#        lwd = 2)

cat("\n=== Example 5: GLM E-Values (Direct Computation) ===\n")
# Compute per-observation e-values directly
y <- example_data$y
X <- model.matrix(~ x1 + x2 + x3, data = example_data)

# Test for effect of x1
evalues <- glm_lr_evalues(
  y = y,
  X = X,
  theta_null = 0,
  theta_alt = 0.5,
  family = gaussian()
)

cat("Summary of per-observation e-values:\n")
print(summary(evalues))

# Plot histogram
hist(evalues, breaks = 30, 
     main = "Distribution of Per-Observation E-Values",
     xlab = "E-Value", col = "lightblue", border = "white")
abline(v = 1, col = "red", lty = 2, lwd = 2)

cat("\n=== Example 6: Mixing E-Values ===\n")
# Apply mixing to ensure E <= 1 under null
mixed <- mix_evalues(evalues, alpha_mix = 0.5, beta_mix = 0.5)

cat("Summary of mixed e-values:\n")
print(summary(mixed))

# Compare distributions
par(mfrow = c(1, 2))
hist(evalues, breaks = 30, main = "Original E-Values",
     xlab = "E-Value", col = "lightblue", border = "white")
abline(v = 1, col = "red", lty = 2, lwd = 2)

hist(mixed, breaks = 30, main = "Mixed E-Values",
     xlab = "Mixed E-Value", col = "lightgreen", border = "white")
abline(v = 1, col = "red", lty = 2, lwd = 2)
par(mfrow = c(1, 1))

cat("\n=== All Examples Completed ===\n")
cat("Note: Examples 1-4 require brms/Stan installation.\n")
cat("Examples 5-6 demonstrate core e-value computation without brms.\n")
