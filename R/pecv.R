#' @title Entrywise Splitting Cross-Validation to determine the number of factors in generalized factor models
#' @description Uses (Penalized) Entrywise Splitting Cross-Validation method to estimate the number of factors in generalized factor models
#' @param resp Observation data matrix (n x p), can be continuous data, count data, or binary data
#' @param C Constraint constant, default is 5
#' @param qmax Maximum number of factors, default is 8
#' @param fold Number of folds in cross-validation, default is 5
#' @param tol_val Convergence tolerance, default is 0.01 (represents 0.01/number of estimated elements)
#' @param theta0 Initial values matrix for factors, generated from uniform distribution if not provided
#' @param A0 Initial values matrix for loadings, generated from uniform distribution if not provided
#' @param seed Random seed, default is 1
#' @param data_type Data type, options are "continuous", "count", "binary". Function will auto-detect if not specified
#' @return Returns a list containing:
#'   \item{ECV}{Number of factors selected by standard ECV}
#'   \item{p1ECV}{Number of factors selected by ECV with penalty 1}
#'   \item{p2ECV}{Number of factors selected by ECV with penalty 2}
#'   \item{p3ECV}{Number of factors selected by ECV with penalty 3}
#'   \item{p4ECV}{Number of factors selected by ECV with penalty 4}
#'   \item{ECV_loss}{CV loss values for each number of factors}
#'   \item{data_type}{Identified data type}
#' @export
#' @examples
#' set.seed(123)
#' # Generate count data
#' n <- 50; p <- 50; q <- 2
#' theta_true <- cbind(1,matrix(runif(n * q,-2,2), n, q))
#' A_true <- matrix(runif(p * (q+1),-2,2), p, (q+1))
#' lambda <- exp(theta_true %*% t(A_true))
#' resp <- matrix(
#' rpois(length(lambda), lambda = as.vector(lambda)),
#' nrow = nrow(lambda), ncol = ncol(lambda))
#' result <- pECV(resp, C = 4, qmax = 4, fold = 5)
#' print(result)

pECV <- function(resp, C = 5, qmax = 8, fold = 5, tol_val = 0.01,
                 theta0 = NULL, A0 = NULL, seed = 1, data_type = NULL) {

  # Set random seed
  set.seed(seed)

  # Get data dimensions
  n <- nrow(resp)
  p <- ncol(resp)

  # Auto-detect data type (if not specified)
  if (is.null(data_type)) {
    data_type <- detect_data_type(resp)
    message(paste("Detected data type:", data_type))
  }

  # Validate data type
  if (!data_type %in% c("continuous", "count", "binary")) {
    stop("data_type must be one of 'continuous', 'count', or 'binary'")
  }

  # Call appropriate function based on data type
  if (data_type == "count") {
    result <- pECV_poisson(resp, n, p, C, qmax, fold, tol_val, theta0, A0)
  } else if (data_type == "continuous") {
    result <- pECV_gaussian(resp, n, p, C, qmax, fold, tol_val, theta0, A0)
  } else if (data_type == "binary") {
    result <- pECV_binary(resp, n, p, C, qmax, fold, tol_val, theta0, A0)
  }

  # Add data type information
  result$data_type <- data_type

  return(result)
}

#' @title Detect data type
#' @description Automatically detect the type of input data
#' @param resp Data matrix
#' @return Data type string
#' @noRd
detect_data_type <- function(resp) {
  # Check if binary data
  unique_vals <- unique(as.vector(resp))
  unique_vals <- unique_vals[!is.na(unique_vals)]

  if (all(unique_vals %in% c(0, 1))) {
    return("binary")
  }

  # Check if count data (non-negative integers)
  if (all(resp >= 0, na.rm = TRUE) &&
      all(resp == round(resp), na.rm = TRUE)) {
    return("count")
  }

  # Otherwise treat as continuous data
  return("continuous")
}

#' @title Matrix splitting function
#' @description Randomly split matrix into K sub-matrices for cross-validation
#' @param n Number of rows
#' @param p Number of columns
#' @param K Number of splits
#' @return List of K logical matrices
#' @noRd
split_matrix <- function(n, p, K) {
  # Create n x p matrix of all 1s
  original_matrix <- matrix(1, nrow = n, ncol = p)

  # Initialize list of K n x p matrices
  split_matrices <- replicate(K, matrix(0, nrow = n, ncol = p), simplify = FALSE)

  # Sample from multinomial distribution for each element and allocate to K matrices
  for (i in 1:n) {
    for (j in 1:p) {
      counts <- rmultinom(1, size = original_matrix[i, j], prob = rep(1/K, K))
      for (k in 1:K) {
        split_matrices[[k]][i, j] <- counts[k]
      }
    }
  }

  # Convert 1s to TRUE and 0s to FALSE in matrices
  split_matrices <- lapply(split_matrices, function(x) x == 1)

  return(split_matrices)
}

#' @title pECV for Poisson factor model
#' @description pECV implementation for count data
#' @noRd
pECV_poisson <- function(resp, n, p, C, qmax, fold, tol_val, theta0, A0) {

  # Initialize factor loading matrices
  if (is.null(theta0)) {
    # First column is intercept (all 1s), remaining columns from uniform distribution [-1,1]
    theta0 <- cbind(rep(1, n), matrix(runif(n * qmax, -1, 1), n, qmax))
  }

  if (is.null(A0)) {
    A0 <- matrix(runif(p * (qmax + 1), -1, 1), p, qmax + 1)
  }

  # Perform cross-validation
  K <- fold
  split_matrices <- split_matrix(n, p, K)

  # Initialize loss storage matrix
  loss_values <- matrix(0, nrow = K, ncol = qmax)

  # Perform K-fold cross-validation
  for (k in 1:K) {
    # Create train/test masks
    test_mask <- split_matrices[[k]]
    train_mask <- Reduce("+", split_matrices[-k]) > 0

    # Create training set
    train_set <- matrix(0, nrow = nrow(resp), ncol = ncol(resp))
    train_set[train_mask] <- resp[train_mask]

    num_train <- sum(train_mask == 1)
    tol_current <- tol_val / num_train

    # Fit models with different numbers of factors
    for (q_est in 1:qmax) {
      # Note: confirm_CJMLE_poisson_cpp function is required here

      # If actual C++ function is available, use:
      jml_res <- confirm_CJMLE_poisson_cpp(
        train_set, train_mask,
        theta0[, 1:(q_est + 1)], A0[, 1:(q_est + 1)],
        matrix(TRUE, p, q_est + 1),
        C = C, tol = tol_current
      )

      jml_result <- JL_poisson(resp, jml_res$theta, jml_res$A, test_mask)
      loss_values[k, q_est] <- -jml_result$lik
    }
  }

  # Calculate different penalty methods
  pen1 <- max(n, p) * log(min(n, p))
  pen2 <- (n + p) * log((n * p) / (n + p))

  # Calculate CV loss
  cv_loss <- colSums(loss_values)

  # Return results
  result <- list(
    ECV = which.min(cv_loss),
    p1ECV = which.min(cv_loss + ((K * K) / (K - 1)) * (1:qmax) * pen1),
    p2ECV = which.min(cv_loss + ((K * K) / (K - 1)) * (1:qmax) * pen2),
    p3ECV = which.min(cv_loss + (K / (K - 1)) * (1:qmax) * pen1),
    p4ECV = which.min(cv_loss + (K / (K - 1)) * (1:qmax) * pen2),
    ECV_loss = cv_loss
  )

  return(result)
}

#' @title pECV for Gaussian factor model
#' @description pECV implementation for continuous data
#' @noRd
pECV_gaussian <- function(resp, n, p, C, qmax, fold, tol_val, theta0, A0) {

  # Initialize factor loading matrices
  if (is.null(theta0)) {
    theta0 <- cbind(rep(1, n), matrix(runif(n * qmax, -1, 1), n, qmax))
  }

  if (is.null(A0)) {
    A0 <- matrix(runif(p * (qmax + 1), -1, 1), p, qmax + 1)
  }

  # Perform cross-validation
  K <- fold
  split_matrices <- split_matrix(n, p, K)

  # Initialize loss storage matrix
  loss_values <- matrix(0, nrow = K, ncol = qmax)

  # Perform K-fold cross-validation
  for (k in 1:K) {
    # Create train/test masks
    test_mask <- split_matrices[[k]]
    train_mask <- Reduce("+", split_matrices[-k]) > 0

    # Create training set
    train_set <- matrix(0, nrow = nrow(resp), ncol = ncol(resp))
    train_set[train_mask] <- resp[train_mask]

    num_train <- sum(train_mask == 1)
    tol_current <- tol_val / num_train

    # Fit models with different numbers of factors
    for (q_est in 1:qmax) {
      # Note: CJMLE_linear function is required here
      # Since this is an external C++ function, a placeholder is provided

      # If actual C++ function is available, use:
      jml_res <- CJMLE_linear(
        train_set, train_mask,
        theta0[, 1:(q_est + 1)], A0[, 1:(q_est + 1)],
        matrix(TRUE, p, q_est + 1), 1,
        C = C, tol = tol_current, F
      )

      jml_result <- JL_gaussian(resp, jml_res$theta, jml_res$A, test_mask)
      loss_values[k, q_est] <- -jml_result$lik
    }
  }

  # Calculate different penalty methods
  pen1 <- max(n, p) * log(min(n, p))
  pen2 <- (n + p) * log((n * p) / (n + p))

  # Calculate CV loss
  cv_loss <- colSums(loss_values)

  # Return results
  result <- list(
    ECV = which.min(cv_loss),
    p1ECV = which.min(cv_loss + ((K * K) / (K - 1)) * (1:qmax) * pen1),
    p2ECV = which.min(cv_loss + ((K * K) / (K - 1)) * (1:qmax) * pen2),
    p3ECV = which.min(cv_loss + (K / (K - 1)) * (1:qmax) * pen1),
    p4ECV = which.min(cv_loss + (K / (K - 1)) * (1:qmax) * pen2),
    ECV_loss = cv_loss
  )

  return(result)
}

#' @title pECV for Binary factor model
#' @description pECV implementation for binary data
#' @noRd
pECV_binary <- function(resp, n, p, C, qmax, fold, tol_val, theta0, A0) {

  # Check if mirtjml package is installed
  if (!requireNamespace("mirtjml", quietly = TRUE)) {
    stop("mirtjml package is required. Please run: install.packages('mirtjml')")
  }

  # Initialize factor loading matrices (note: binary model has different dimensions)
  if (is.null(theta0)) {
    theta0 <- matrix(runif(n * qmax, -1, 1), n, qmax)
  }

  if (is.null(A0)) {
    A0 <- matrix(runif(p * qmax, -1, 1), p, qmax)
  }

  # Initialize intercept terms
  d0 <- runif(p, -1, 1)

  # Perform cross-validation
  K <- fold
  split_matrices <- split_matrix(n, p, K)

  # Initialize loss storage matrix
  loss_values <- matrix(0, nrow = K, ncol = qmax)

  # Perform K-fold cross-validation
  for (k in 1:K) {
    # Create train/test masks
    test_mask <- split_matrices[[k]]
    train_mask <- Reduce("+", split_matrices[-k]) > 0

    # Create training and test sets
    train_set <- matrix(NA, nrow = nrow(resp), ncol = ncol(resp))
    train_set[train_mask] <- resp[train_mask]
    test_set <- matrix(NA, nrow = nrow(resp), ncol = ncol(resp))
    test_set[test_mask] <- resp[test_mask]

    # Fit models with different numbers of factors
    for (q_est in 1:qmax) {
      tryCatch({
        # Note: mirtjml_expr function from mirtjml package is required here
        # If package is installed, use:
        jml_res <- mirtjml::mirtjml_expr(
          train_set, q_est,
          theta0 = theta0[, 1:q_est],
          A0 = A0[, 1:q_est], d0 = d0,
          tol = 0.01, cc = C,
          print_proc = FALSE
        )


        jml_result <- JL_binary(test_set, jml_res$A_hat, jml_res$theta_hat, jml_res$d_hat)
        loss_values[k, q_est] <- -jml_result$lik
      }, error = function(e) {
        loss_values[k, q_est] <- NA
      })
    }
  }

  # Calculate different penalty methods
  pen1 <- max(n, p) * log(min(n, p))
  pen2 <- (n + p) * log((n * p) / (n + p))

  # Calculate CV loss (handle NA values)
  cv_loss <- apply(loss_values, 2, function(x) {
    non_na <- x[!is.na(x)]
    if (length(non_na) == 0) Inf else sum(non_na)
  })

  # Return results
  result <- list(
    ECV = which.min(cv_loss),
    p1ECV = which.min(cv_loss + (K / (K - 1)) * (1:qmax) * pen1),
    p2ECV = which.min(cv_loss + (K / (K - 1)) * (1:qmax) * pen2),
    p3ECV = which.min(cv_loss + (1 / (K - 1)) * (1:qmax) * pen1),
    p4ECV = which.min(cv_loss + (1 / (K - 1)) * (1:qmax) * pen2),
    ECV_loss = cv_loss
  )

  return(result)
}

#' @title Joint likelihood function for Poisson model
#' @description Calculate log-likelihood for Poisson factor model
#' @noRd
JL_poisson <- function(resp, theta, A, nonmis_ind) {
  temp <- theta %*% t(A)
  temp1 <- resp * temp - exp(temp)
  lik <- sum(temp1[nonmis_ind])
  list(lik = lik, M = temp)
}

#' @title Joint likelihood function for Gaussian model
#' @description Calculate log-likelihood for Gaussian factor model
#' @noRd
JL_gaussian <- function(resp, theta, A, nonmis_ind) {
  temp <- theta %*% t(A)
  lik <- sum((log(dnorm(resp, mean = temp, sd = 1)))[nonmis_ind])
  list(lik = lik, M = temp)
}

#' @title Joint likelihood function for Binary model
#' @description Calculate log-likelihood for Binary factor model
#' @noRd
JL_binary <- function(data, A, Theta, d) {
  N <- nrow(data)
  temp <- Theta %*% t(A) + rep(1, N) %*% t(d)
  M <- temp
  prob <- 1/(1 + exp(-temp))
  temp <- data * log(prob) + (1 - data) * log(1 - prob)
  temp[is.na(temp)] <- 0
  list(lik = sum(temp), M = M)
}

# ===== Auxiliary functions for generating example data =====

#' @title Generate continuous data example
#' @description Generate simulated data based on Gaussian factor model
#' @param n Number of observations
#' @param p Number of variables
#' @param q True number of factors
#' @param noise_sd The variance of noise
#' @export
generate_continuous_data <- function(n = 100, p = 50, q = 3, noise_sd = 1) {
  # Generate true factors and loadings
  theta_true <- cbind(rep(1, n), matrix(runif(n * q, -4, 4), n, q))
  A_true <- matrix(runif(p * (q + 1), -4, 4), p, q + 1)
  # Generate observed data
  signal <- theta_true %*% t(A_true)
  noise <- matrix(rnorm(n * p, sd = noise_sd), n, p)
  resp <- signal + noise

  return(list(resp = resp, true_q = q, theta_true = theta_true, A_true = A_true))
}

#' @title Generate count data example
#' @description Generate simulated data based on Poisson factor model
#' @param n Number of observations
#' @param p Number of variables
#' @param q True number of factors
#' @export
generate_count_data <- function(n = 100, p = 50, q = 3) {
  # Generate true factors and loadings
  theta_true <- cbind(rep(1, n), matrix(runif(n * q, -2, 2), n, q))
  A_true <- matrix(runif(p * (q + 1), -2, 2), p, q + 1)

  # Generate intensity parameters
  lambda <- exp(theta_true %*% t(A_true))

  # Generate Poisson observations
  resp <- matrix(rpois(n * p, lambda = as.vector(lambda)), n, p)

  return(list(resp = resp, true_q = q, theta_true = theta_true, A_true = A_true))
}

#' @title Generate binary data example
#' @description Generate simulated data based on Binary factor model
#' @param n Number of observations
#' @param p Number of variables
#' @param q True number of factors
#' @export
generate_binary_data <- function(n = 100, p = 50, q = 3) {
  # Generate true factors and loadings
  theta_true <- matrix(runif(n * q, -8, 8), n, q)
  A_true <- matrix(runif(p * (q), -8, 8), p, q)
  d_true <- runif(p, -8, 8)

  # Generate probabilities
  logit <- theta_true %*% t(A_true) + rep(1, n) %*% t(d_true)
  prob <- 1 / (1 + exp(-logit))

  # Generate binary observations
  resp <- matrix(rbinom(n * p, 1, prob = as.vector(prob)), n, p)

  return(list(resp = resp, true_q = q, theta_true = theta_true,
              A_true = A_true, d_true = d_true))
}




















