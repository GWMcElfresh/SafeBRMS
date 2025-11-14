# SafeBRMS Usage Examples

## Installation

```r
# Install dependencies first
install.packages("brms")

# Install SafeBRMS (from source or GitHub)
# R CMD INSTALL SafeBRMS
# or
# devtools::install_github("GWMcElfresh/SafeBRMS")

library(SafeBRMS)
```

## Example 1: Basic Workflow

```r
# Load example data
data(example_data)
summary(example_data)

# Fit Bayesian regression model
fit <- fit_safe_brms(
  formula = y ~ x1 + x2 + x3,
  data = example_data,
  chains = 2,
  iter = 1000
)

# Test single coefficient
result <- single_coef_test(fit, coef_name = "x1", null_value = 0)
print(result)
# E-value > 20 indicates rejection at α = 0.05
```

## Example 2: Multiple Testing

```r
# Test all coefficients with FDR control
results <- multiple_testing_bh(
  fit = fit,
  null_value = 0,
  fdr_level = 0.05
)

print(results)
# Shows which coefficients are significant after BH correction
```

## Example 3: E-Posterior

```r
# Compute e-posterior for a coefficient
epost <- compute_eposterior(fit, coef_name = "x1")

# Plot
plot(epost$theta, epost$eposterior, type = "l",
     xlab = "Theta", ylab = "E-Posterior Density",
     main = "E-Posterior for x1")
abline(v = epost$MAP_theta, col = "red", lty = 2)
```

## Example 4: Direct E-Value Computation

```r
# Prepare data
y <- example_data$y
X <- model.matrix(~ x1 + x2 + x3, data = example_data)

# Compute per-observation e-values
evalues <- glm_lr_evalues(
  y = y,
  X = X,
  theta_null = 0,
  theta_alt = 0.5,
  family = gaussian()
)

# Summarize
summary(evalues)
hist(evalues, main = "Distribution of E-Values")
```

## Example 5: Mixing E-Values

```r
# Apply mixing to ensure calibration
mixed <- mix_evalues(evalues, alpha_mix = 0.5, beta_mix = 0.5)

# Compare
par(mfrow = c(1, 2))
hist(evalues, main = "Original E-Values")
hist(mixed, main = "Mixed E-Values")
```

## Example 6: Simulation Study

```r
# Simulate data under null
set.seed(123)
n <- 100
x <- rnorm(n)
y_null <- rnorm(n)  # No relationship
X <- cbind(1, x)

# Compute e-values
e_null <- glm_lr_evalues(y_null, X, theta_null = 0, theta_alt = 0.3)
mean(e_null)  # Should be ≤ 1

# Simulate data under alternative
y_alt <- 0.5 * x + rnorm(n)  # True effect
e_alt <- glm_lr_evalues(y_alt, X, theta_null = 0, theta_alt = 0.5)
mean(e_alt)  # Should be > 1
```

## Example 7: Sequential Testing

```r
# Optional stopping is valid with e-values
cumulative_e <- 1
for (i in 1:100) {
  # Generate new observation
  x_i <- rnorm(1)
  y_i <- 0.5 * x_i + rnorm(1)
  
  # Update e-value
  e_i <- exp(-0.5 * (y_i - 0.5 * x_i)^2 + 0.5 * y_i^2)
  cumulative_e <- cumulative_e * e_i
  
  # Can stop anytime if cumulative_e > 20
  if (cumulative_e > 20) {
    cat("Reject null at observation", i, "\n")
    break
  }
}
```

## Interpretation Guide

### E-Values
- **E-value = 1**: No evidence against null
- **E-value > 1**: Evidence against null
- **E-value > 1/α**: Reject null at level α
- **E-value > 20**: Reject at α = 0.05

### E-Posterior
- Distribution over parameter values
- MAP estimate = mode of e-posterior
- Wider distribution = more uncertainty
- Can be used for credible intervals

### Multiple Testing
- `reject = TRUE`: Significant after BH correction
- `adjusted_pvalue`: BH-adjusted significance level
- Results sorted by significance
- Controls FDR at specified level

## Tips

1. **Use more MCMC iterations** for better posterior estimates in real analyses
2. **Check convergence** of brms models (Rhat, effective sample size)
3. **Visualize e-posteriors** to understand parameter uncertainty
4. **Consider mixing parameters** based on prior knowledge
5. **Use product of e-values** for joint tests of multiple parameters

## References

See the package vignette for detailed theory and methodology:
```r
vignette("safebrms_workflow", package = "SafeBRMS")
```

For examples demonstrating safestats principles:
```r
source(system.file("examples/safestats_principles.R", package = "SafeBRMS"))
```
