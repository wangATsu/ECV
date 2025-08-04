# ECV: Entrywise Splitting Cross-Validation for Generalized Factor Models

[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/wangATsu/ECV/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/wangATsu/ECV/actions/workflows/R-CMD-check.yml)

## Overview

The **ECV** (Entrywise Splitting Cross-Validation) package provides data-driven methods for determining the number of factors in generalized factor models. It implements both standard ECV and penalized ECV (pECV) approaches.

The package also handles missing data through `pECV.miss()` and provides data-driven estimation of the constraint constant C through `estimate_C()` functions.

## Key Features

- **Automatic data type detection**: The package automatically identifies whether your data is continuous, count, or binary
- **Missing data support**: Handles datasets with missing values using specialized algorithms
- **Data-driven parameter estimation**: Provides methods to estimate the constraint constant C from data
- **Multiple penalty methods**: Implements various penalized versions (p1ECV, p2ECV, p3ECV, p4ECV) for improved selection
- **Efficient implementation**: Core algorithms implemented in C++ for speed

## Installation

You can install the development version of ECV from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install ECV from GitHub
devtools::install_github("wangATsu/ECV")

# Load the package
library(ECV)
```

## Basic Usage

### Main Function: `pECV()`

The primary function for factor number selection:

```r
pECV(resp, C = 5, qmax = 8, fold = 5, tol_val = 0.01, 
     theta0 = NULL, A0 = NULL, seed = 1, data_type = NULL)
```

**Parameters:**
- `resp`: Observation data matrix (n × p)
- `C`: Constraint constant (default: 5)
- `qmax`: Maximum number of factors to consider (default: 8)
- `fold`: Number of folds in cross-validation (default: 5)
- `tol_val`: Convergence tolerance (default: 0.01)
- `theta0`, `A0`: Initial values for factors and loadings (optional)
- `seed`: Random seed for reproducibility
- `data_type`: "continuous", "count", or "binary" (auto-detected if NULL)

**Returns:** A list containing:
- `ECV`: Number of factors selected by standard ECV
- `p1ECV`, `p2ECV`, `p3ECV`, `p4ECV`: Factors number selected by penalized versions
- `ECV_loss`: Cross-validation loss values for each factor number
- `data_type`: Detected/specified data type

## Examples

### Example 1

```r
library(ECV)
set.seed(123)

# Generate count data from a Poisson factor model
n <- 50  # number of observations
p <- 50  # number of variables
q <- 2   # true number of factors

# Generate true parameters
theta_true <- cbind(1, matrix(runif(n * q, -2, 2), n, q))
A_true <- matrix(runif(p * (q + 1), -2, 2), p, (q + 1))

# Generate Poisson observations
lambda <- exp(theta_true %*% t(A_true))
resp <- matrix(rpois(length(lambda), lambda = as.vector(lambda)),
               nrow = nrow(lambda), ncol = ncol(lambda))

# Apply pECV to select the number of factors
result <- pECV(resp, C = 4, qmax = 4, fold = 5)

# View results
print(result$ECV)    # Selected number of factors
print(result$p1ECV)  # Penalized version result
print(result$ECV_loss)  # CV losses for each factor number
```

### Example 2
```r
set.seed(123)
# Generate count data with missing values
n <- 50; p <- 50; q <- 2
theta_true <- cbind(1,matrix(runif(n * q,-2,2), n, q))
A_true <- matrix(runif(p * (q+1),-2,2), p, (q+1))
lambda <- exp(theta_true %*% t(A_true))
resp <- matrix(
rpois(length(lambda), lambda = as.vector(lambda)),
nrow = nrow(lambda), ncol = ncol(lambda))
# Introduce 5% missing values
miss_idx <- sample(1:(n*p), size = 0.05*n*p)
resp[miss_idx] <- NA
result <- pECV.miss(resp, C = 4, qmax = 4, fold = 5)
print(result)
```

### Data-Driven Estimation of C

```r
# For count data (transform first)
count_data <- generate_count_data(n = 50, p = 50, q = 2)
count_transformed <- log(count_data$resp + 1)
C_est_count <- estimate_C(count_transformed, qmax = 4)
```

## Additional Functions

### Data Generation Functions
- `generate_continuous_data()`: Generate Gaussian factor model data
- `generate_count_data()`: Generate Poisson factor model data
- `generate_binary_data()`: Generate binary factor model data
- `generate_*_data_miss()`: Versions with missing values

### Parameter Estimation
- `estimate_C()`: Estimate constraint constant for continuous/count data
- `estimate_C_binary()`: Estimate constraint constant for binary data

### Missing Data Support
- `pECV.miss()`: Main function for data with missing values
- `introduce_missing()`: Add missing values to complete data

## Methodological Details

The ECV method works by:
1. Randomly splitting each data entry across K folds
2. Training models with different numbers of factors on K-1 folds
3. Evaluating on the held-out fold
4. Selecting the number of factors that minimizes cross-validation loss

The penalized versions (p1ECV through p4ECV) add different penalty terms to favor more parsimonious models (consistent estimation).


## Bug Reports and Issues

If you encounter any bugs or issues, please report them on the [GitHub Issues page](https://github.com/wangATsu/ECV/issues). When reporting, please include:

1. A clear description of the issue
2. Reproducible code example (preferably using the [reprex](https://reprex.tidyverse.org/) package)
3. Your session information (`sessionInfo()`)
4. Any relevant data files or error messages

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License

This package is licensed under GPL-3. See the [LICENSE](LICENSE) file for details.
