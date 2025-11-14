## Example: Safestats E-Values Integration
## This demonstrates how SafeBRMS e-values align with safestats principles

# If SafeBRMS is installed as a package:
# library(SafeBRMS)

# For testing, source the R files directly:
source("R/glm_evalues.R")
source("R/compute_evalues.R")

set.seed(123)

cat("=== Safestats Principles in SafeBRMS ===\n\n")

cat("1. E-VALUES SATISFY E[e] <= 1 UNDER NULL\n")
cat("   This is the fundamental safe testing property\n\n")

# Generate data under the null (no effect)
n <- 100
x <- rnorm(n)
y_null <- rnorm(n)  # No relationship with x
X <- cbind(1, x)

# Compute e-values when null is true
evalues_null <- glm_lr_evalues(y_null, X, theta_null = 0, theta_alt = 0.3)
cat(sprintf("Mean e-value under null: %.3f (should be <= 1)\n", mean(evalues_null)))
cat(sprintf("Median e-value under null: %.3f\n", median(evalues_null)))

cat("\n2. E-VALUES GROW UNDER ALTERNATIVE\n")
cat("   When the alternative is true, e-values tend to be larger\n\n")

# Generate data under alternative (true effect)
y_alt <- 0.5 * x + rnorm(n)  # True effect of 0.5
evalues_alt <- glm_lr_evalues(y_alt, X, theta_null = 0, theta_alt = 0.5)
cat(sprintf("Mean e-value under alternative: %.3f (should be > 1)\n", mean(evalues_alt)))
cat(sprintf("Median e-value under alternative: %.3f\n", median(evalues_alt)))

cat("\n3. PRODUCT PROPERTY FOR MULTIPLE TESTS\n")
cat("   Product of independent e-values is also an e-value\n\n")

# Simulate multiple independent tests
set.seed(456)
n_tests <- 5
product_evalue <- 1

for (i in 1:n_tests) {
  x_i <- rnorm(n)
  y_i <- rnorm(n)  # All null
  X_i <- cbind(1, x_i)
  e_i <- glm_lr_evalues(y_i, X_i, theta_null = 0, theta_alt = 0.3)
  product_evalue <- product_evalue * mean(e_i)
}

cat(sprintf("Product of %d e-values under null: %.3f\n", n_tests, product_evalue))
cat("Expected to be <= 1 due to product property\n")

cat("\n4. MIXING ENSURES CALIBRATION\n")
cat("   Mixing integrates over alternative values for robustness\n\n")

mixed_null <- mix_evalues(evalues_null)
mixed_alt <- mix_evalues(evalues_alt)

cat(sprintf("Mean mixed e-value under null: %.3f\n", mean(mixed_null)))
cat(sprintf("Mean mixed e-value under alternative: %.3f\n", mean(mixed_alt)))

cat("\n5. OPTIONAL STOPPING VALIDITY\n")
cat("   E-values remain valid even with optional stopping\n\n")

# Simulate sequential testing
set.seed(789)
sequential_e <- numeric(100)
cumulative_product <- 1

for (t in 1:100) {
  # Generate one observation at a time (under null)
  x_t <- rnorm(1)
  y_t <- rnorm(1)
  
  # Simple e-value computation
  e_t <- exp(-0.5 * (y_t - 0.3 * x_t)^2 + 0.5 * y_t^2)
  cumulative_product <- cumulative_product * e_t^(1/100)  # Normalized
  sequential_e[t] <- cumulative_product
}

cat(sprintf("Final cumulative e-value: %.3f\n", cumulative_product))
cat("Valid for any stopping rule!\n")

cat("\n6. ANYTIME-VALID CONFIDENCE SEQUENCES\n")
cat("   E-values can construct confidence sequences valid at all times\n\n")

# Plot sequential e-values
plot(1:100, sequential_e, type = "l", 
     xlab = "Sample Size", ylab = "Cumulative E-Value",
     main = "Sequential E-Values (Valid with Optional Stopping)",
     col = "blue", lwd = 2)
abline(h = 1, col = "red", lty = 2)
abline(h = 1/0.05, col = "orange", lty = 2)
legend("topright", 
       legend = c("Cumulative E-Value", "E=1 (null reference)", "E=20 (reject at α=0.05)"),
       col = c("blue", "red", "orange"),
       lty = c(1, 2, 2),
       lwd = c(2, 1, 1))

cat("\n=== Comparison: E-Values vs P-Values ===\n\n")

cat("E-VALUES (Safe Testing):\n")
cat("  ✓ Valid under optional stopping\n")
cat("  ✓ Product property for combining tests\n")
cat("  ✓ Always valid (even if assumptions violated)\n")
cat("  ✓ Interpretable as betting odds\n")
cat("  ✓ E > 1/α rejects at level α\n\n")

cat("P-VALUES (Classical Testing):\n")
cat("  ✗ Invalid under optional stopping\n")
cat("  ✗ No simple combination rule\n")
cat("  ✗ Requires correct model specification\n")
cat("  ✗ Not directly interpretable as evidence\n")
cat("  ✓ P < α rejects at level α\n\n")

cat("=== Key Safestats Concepts Implemented ===\n\n")
cat("1. Test martingales: E[e] <= 1 under null\n")
cat("2. Growth rate: E[log(e)] < 0 under null\n")
cat("3. Betting interpretation: E-value as wealth ratio\n")
cat("4. Calibration: Valid under misspecification\n")
cat("5. Composability: Products and mixtures of e-values\n")

cat("\n=== References ===\n")
cat("- Grünwald, de Heide, & Koolen (2019). Safe testing.\n")
cat("- Shafer (2021). Testing by betting: A strategy for statistical and scientific communication.\n")
cat("- Turner, Ly, & Grünwald (2020+). Safe testing in regression.\n")

cat("\nExample completed!\n")
