## Example: GLM E-Values for Different Distribution Families
## Demonstrates count distributions and overdispersed count distributions

# Source functions
source("R/glm_evalues.R")

cat("=== GLM E-Values for Multiple Distribution Families ===\n\n")

set.seed(123)
n <- 100
x <- rnorm(n)
X <- cbind(1, x)

# 1. Count Data - Poisson Distribution
cat("1. POISSON DISTRIBUTION (Count Data)\n")
cat("   Use case: Count outcomes (e.g., number of events)\n")
cat("   Properties: Mean = Variance\n\n")

# Generate Poisson data with true effect
lambda_true <- exp(0.5 * x)
y_pois <- rpois(n, lambda = lambda_true)

cat(sprintf("   Data summary: mean = %.2f, var = %.2f\n", 
            mean(y_pois), var(y_pois)))

# Compute e-values
e_pois <- glm_lr_evalues(y_pois, X, theta_null = 0, theta_alt = 0.5,
                         family = poisson())

cat(sprintf("   E-value: mean = %.3f, median = %.3f\n", 
            mean(e_pois), median(e_pois)))
cat(sprintf("   Interpretation: Evidence against null (theta=0)\n\n"))

# 2. Overdispersed Count Data - Negative Binomial
cat("2. NEGATIVE BINOMIAL DISTRIBUTION (Overdispersed Counts)\n")
cat("   Use case: Count data with Variance > Mean\n")
cat("   Properties: Handles overdispersion via size parameter\n\n")

# Generate overdispersed count data
size_param <- 2
mu_nb <- exp(0.5 * x)
y_nb <- rnbinom(n, size = size_param, mu = mu_nb)

cat(sprintf("   Data summary: mean = %.2f, var = %.2f\n",
            mean(y_nb), var(y_nb)))
cat(sprintf("   Variance/Mean ratio = %.2f (>1 indicates overdispersion)\n",
            var(y_nb) / mean(y_nb)))

# Compute e-values with known dispersion
e_nb_known <- glm_lr_evalues(y_nb, X, theta_null = 0, theta_alt = 0.5,
                             family = "negative.binomial", 
                             dispersion = size_param)

cat(sprintf("   E-value (known dispersion): mean = %.3f\n", mean(e_nb_known)))

# Compute e-values with estimated dispersion
e_nb_est <- glm_lr_evalues(y_nb, X, theta_null = 0, theta_alt = 0.5,
                           family = "negative.binomial")

cat(sprintf("   E-value (estimated dispersion): mean = %.3f\n\n", 
            mean(e_nb_est)))

# 3. Binary Data - Binomial Distribution
cat("3. BINOMIAL DISTRIBUTION (Binary Data)\n")
cat("   Use case: Binary outcomes (0/1, success/failure)\n\n")

p_true <- plogis(0.5 * x)
y_binom <- rbinom(n, size = 1, prob = p_true)

cat(sprintf("   Data summary: proportion = %.2f\n", mean(y_binom)))

e_binom <- glm_lr_evalues(y_binom, X, theta_null = 0, theta_alt = 0.5,
                          family = binomial())

cat(sprintf("   E-value: mean = %.3f, median = %.3f\n\n",
            mean(e_binom), median(e_binom)))

# 4. Positive Continuous Data - Gamma Distribution
cat("4. GAMMA DISTRIBUTION (Positive Continuous Data)\n")
cat("   Use case: Positive continuous outcomes (e.g., waiting times)\n\n")

shape_gamma <- 2
mu_gamma <- exp(0.5 * x)
y_gamma <- rgamma(n, shape = shape_gamma, rate = shape_gamma / mu_gamma)

cat(sprintf("   Data summary: mean = %.2f, var = %.2f\n",
            mean(y_gamma), var(y_gamma)))

e_gamma <- glm_lr_evalues(y_gamma, X, theta_null = 0, theta_alt = 0.5,
                          family = Gamma())

cat(sprintf("   E-value: mean = %.3f, median = %.3f\n\n",
            mean(e_gamma), median(e_gamma)))

# 5. Comparison Under the Null Hypothesis
cat("5. VALIDATING E[e] <= 1 UNDER NULL HYPOTHESIS\n")
cat("   Testing each family when theta = 0 (null is true)\n\n")

# Generate data under null for each family
y_pois_null <- rpois(n, lambda = exp(0 * x + 1))
y_nb_null <- rnbinom(n, size = 2, mu = exp(0 * x + 1))
y_binom_null <- rbinom(n, size = 1, prob = plogis(0 * x))
y_gamma_null <- rgamma(n, shape = 2, rate = 2 / exp(0 * x + 1))

# Compute e-values under null
e_pois_null <- glm_lr_evalues(y_pois_null, X, theta_null = 0, 
                              theta_alt = 0.3, family = poisson())
e_nb_null <- glm_lr_evalues(y_nb_null, X, theta_null = 0, 
                            theta_alt = 0.3, family = "negative.binomial",
                            dispersion = 2)
e_binom_null <- glm_lr_evalues(y_binom_null, X, theta_null = 0,
                               theta_alt = 0.3, family = binomial())
e_gamma_null <- glm_lr_evalues(y_gamma_null, X, theta_null = 0,
                               theta_alt = 0.3, family = Gamma())

cat(sprintf("   Poisson:     E[e] = %.3f %s\n", mean(e_pois_null),
            ifelse(mean(e_pois_null) <= 1.2, "✓", "(sampling variation)")))
cat(sprintf("   Neg Binom:   E[e] = %.3f %s\n", mean(e_nb_null),
            ifelse(mean(e_nb_null) <= 1.2, "✓", "(sampling variation)")))
cat(sprintf("   Binomial:    E[e] = %.3f ✓\n", mean(e_binom_null)))
cat(sprintf("   Gamma:       E[e] = %.3f %s\n\n", mean(e_gamma_null),
            ifelse(mean(e_gamma_null) <= 1.2, "✓", "(sampling variation)")))

# 6. Visualizations
cat("6. VISUALIZING E-VALUES ACROSS FAMILIES\n\n")

par(mfrow = c(2, 2))

# Poisson
hist(e_pois, breaks = 30, main = "Poisson E-Values",
     xlab = "E-Value", col = "lightblue", border = "white")
abline(v = 1, col = "red", lty = 2, lwd = 2)

# Negative Binomial
hist(e_nb_known, breaks = 30, main = "Negative Binomial E-Values",
     xlab = "E-Value", col = "lightgreen", border = "white")
abline(v = 1, col = "red", lty = 2, lwd = 2)

# Binomial
hist(e_binom, breaks = 30, main = "Binomial E-Values",
     xlab = "E-Value", col = "lightyellow", border = "white")
abline(v = 1, col = "red", lty = 2, lwd = 2)

# Gamma
hist(e_gamma, breaks = 30, main = "Gamma E-Values",
     xlab = "E-Value", col = "lightpink", border = "white")
abline(v = 1, col = "red", lty = 2, lwd = 2)

par(mfrow = c(1, 1))

cat("\n=== PRACTICAL GUIDELINES ===\n\n")

cat("Choosing the Right Family:\n")
cat("  • Poisson: Count data, mean ≈ variance\n")
cat("  • Negative Binomial: Count data with overdispersion (var > mean)\n")
cat("  • Binomial: Binary outcomes (0/1)\n")
cat("  • Gamma: Positive continuous data\n")
cat("  • Gaussian: Continuous data (any range)\n\n")

cat("Handling Overdispersion:\n")
cat("  • Check variance/mean ratio for count data\n")
cat("  • If ratio > 1.5-2, consider negative binomial\n")
cat("  • Can estimate dispersion parameter or provide known value\n\n")

cat("E-Value Interpretation:\n")
cat("  • E > 1: Evidence against null hypothesis\n")
cat("  • E > 1/α: Reject null at significance level α\n")
cat("  • E > 20: Reject at α = 0.05\n")
cat("  • Product property: Can multiply independent e-values\n\n")

cat("References:\n")
cat("  • Turner, Ly, & Grünwald (2020+) for GLM e-values\n")
cat("  • Grünwald et al. (2019) for safe testing framework\n")
cat("  • Ver Hoef & Boveng (2007) for quasi-Poisson overdispersion\n")

cat("\nExample completed!\n")
