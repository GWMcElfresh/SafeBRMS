# SafeBRMS

Safe Bayesian Regression Model Selection with E-Values

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

SafeBRMS is an R package that implements safe Bayesian regression model selection using e-values. It combines the power of [brms](https://paul-buerkner.github.io/brms/) for Bayesian modeling with rigorous safe testing procedures based on the Turner/Ly/Grunwald e-variable construction.

### Key Features

- **Fit Bayesian regression models** using brms with safe testing in mind
- **Compute e-values** for hypothesis testing with type-I error guarantees
- **Per-observation GLM e-variables** based on likelihood ratios for exponential families
- **Mixing to ensure E ≤ 1** under the null hypothesis
- **E-posteriors** over parameter grids for robust inference
- **Single and multiple coefficient tests** with proper error control
- **Benjamini-Hochberg correction** for multiple testing using e-values

## Installation

```r
# Install from GitHub
# devtools::install_github("GWMcElfresh/SafeBRMS")

# Load the package
library(SafeBRMS)
```

## Quick Start

```r
library(SafeBRMS)

# Load example data
data(example_data)

# Fit a Bayesian regression model
fit <- fit_safe_brms(
  formula = y ~ x1 + x2 + x3,
  data = example_data,
  chains = 2,
  iter = 1000
)

# Test a single coefficient
result <- single_coef_test(fit, coef_name = "x1")
print(result)

# Multiple testing with FDR control
results <- multiple_testing_bh(fit, fdr_level = 0.05)
print(results)

# Compute e-posterior
epost <- compute_eposterior(fit, coef_name = "x1")
plot(epost$theta, epost$eposterior, type = "l")
```

## What are E-Values?

E-values (evidence values) are alternatives to p-values that provide:

1. **Type-I error guarantees**: E[e] ≤ 1 under the null hypothesis
2. **Product property**: Can multiply e-values while maintaining validity
3. **Optional stopping**: Valid for sequential testing
4. **Model robustness**: Remain valid under model misspecification

At significance level α, we reject the null when e-value > 1/α. For example, at α = 0.05, we reject when e-value > 20.

## Documentation

See the package vignette for a complete workflow:

```r
vignette("safebrms_workflow", package = "SafeBRMS")
```

## Functions

- `fit_safe_brms()`: Fit a Bayesian regression model using brms
- `compute_evalues()`: Compute e-values for coefficients
- `compute_eposterior()`: Compute e-posterior over parameter grid
- `glm_lr_evalues()`: Per-observation GLM likelihood ratio e-values
- `mix_evalues()`: Mix e-values to ensure E ≤ 1 under null
- `single_coef_test()`: Test a single coefficient
- `multiple_testing_bh()`: Multiple testing with Benjamini-Hochberg correction

## Example Dataset

The package includes `example_data`, a simulated dataset with:
- 100 observations
- Response `y` and three predictors `x1`, `x2`, `x3`
- True effects: y = 0.5 * x1 + 0.3 * x2 + ε

## References

- Grünwald, P., de Heide, R., & Koolen, W. (2019). Safe testing. arXiv preprint.
- Turner, D., Ly, A., & Grünwald, P. (2020+). Safe testing and e-values in regression.
- Bürkner, P. C. (2017). brms: An R package for Bayesian multilevel models using Stan.

## License

GPL (>= 3)
