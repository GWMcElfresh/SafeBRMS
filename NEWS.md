# SafeBRMS 0.1.0

## Initial Release

This is the first release of SafeBRMS, providing safe Bayesian regression model selection with e-values.

### Features

* `fit_safe_brms()`: Wrapper for brms model fitting optimized for safe testing
* `glm_lr_evalues()`: Per-observation GLM likelihood ratio e-values for exponential families
* `mix_evalues()`: Mixing procedure to ensure E[e] ≤ 1 under the null hypothesis
* `compute_evalues()`: Coefficient e-values using Turner/Ly/Grunwald construction
* `compute_eposterior()`: E-posterior distributions over parameter grids
* `single_coef_test()`: Safe hypothesis testing for single coefficients
* `multiple_testing_bh()`: Multiple testing with Benjamini-Hochberg FDR control

### Documentation

* Complete package vignette demonstrating workflow and reproducibility
* Roxygen2 documentation for all exported functions
* Example dataset (`example_data`) with 100 observations
* Two example scripts in `inst/examples/`:
  - `basic_workflow.R`: Complete analysis workflow
  - `safestats_principles.R`: Demonstration of safe testing principles

### Tests

* Unit tests for core e-value computation functions
* Data validation tests
* All tests use testthat framework

### References

Implementation based on:
- Grünwald, de Heide, & Koolen (2019). Safe testing.
- Turner, Ly, & Grünwald (2020+). Safe testing and e-values in regression.
- Bürkner (2017). brms: Bayesian regression models using Stan.
