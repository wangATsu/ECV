
#' @title Entrywise Splitting Cross-Validation to Determine the Number of Factors with Missing Data
#' @description Uses (Penalized) Entrywise Splitting Cross-Validation method to estimate the number of factors
#'              in generalized factor models when data contains missing values
#' @param resp Observation data matrix (n x p) with missing values as NA, can be continuous data, count data, or binary data
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
#'   \item{miss_percent}{Percentage of missing data}
#' @export
#' @examples
#' set.seed(123)
#' # Generate count data with missing values
#' n <- 50; p <- 50; q <- 2
#' theta_true <- cbind(1,matrix(runif(n * q,-2,2), n, q))
#' A_true <- matrix(runif(p * (q+1),-2,2), p, (q+1))
#' lambda <- exp(theta_true %*% t(A_true))
#' resp <- matrix(
#' rpois(length(lambda), lambda = as.vector(lambda)),
#' nrow = nrow(lambda), ncol = ncol(lambda))
#' # Introduce 5% missing values
#' miss_idx <- sample(1:(n*p), size = 0.05*n*p)
#' resp[miss_idx] <- NA
#' result <- pECV.miss(resp, C = 4, qmax = 4, fold = 5)
#' print(result)

pECV.miss <- function(resp, C = 5, qmax = 8, fold = 5, tol_val = 0.01,
                      theta0 = NULL, A0 = NULL, seed = 1, data_type = NULL) {

  # Set random seed
  set.seed(seed)

  # Get data dimensions
  n <- nrow(resp)
  p <- ncol(resp)

  # Create missing indicator matrix (1 = observed, 0 = missing)
  miss_ind <- !is.na(resp)
  miss_ind_numeric <- matrix(as.numeric(miss_ind), n, p)

  # Calculate missing percentage
  miss_percent <- (1 - sum(miss_ind) / (n * p)) * 100

  # Replace NA with 0 for computation (will be masked by miss_ind)
  resp[!miss_ind] <- 0

  # Auto-detect data type (if not specified) - only check non-missing values
  if (is.null(data_type)) {
    resp_non_missing <- resp[miss_ind]
    data_type <- detect_data_type_miss(resp_non_missing)
    message(paste("Detected data type:", data_type))
    message(paste("Missing data percentage:", round(miss_percent, 2), "%"))
  }

  # Validate data type
  if (!data_type %in% c("continuous", "count", "binary")) {
    stop("data_type must be one of 'continuous', 'count', or 'binary'")
  }

  # Call appropriate function based on data type
  if (data_type == "count") {
    result <- pECV.miss_poisson(resp, miss_ind_numeric, n, p, C, qmax, fold, tol_val, theta0, A0)
  } else if (data_type == "continuous") {
    result <- pECV.miss_gaussian(resp, miss_ind_numeric, n, p, C, qmax, fold, tol_val, theta0, A0)
  } else if (data_type == "binary") {
    result <- pECV.miss_binary(resp, miss_ind_numeric, n, p, C, qmax, fold, tol_val, theta0, A0)
  }

  # Add data type and missing percentage information
  result$data_type <- data_type
  result$miss_percent <- miss_percent

  return(result)
}

#' @noRd
detect_data_type_miss <- function(resp_non_missing) {
  # Check if binary data
  unique_vals <- unique(resp_non_missing)

  if (all(unique_vals %in% c(0, 1))) {
    return("binary")
  }

  # Check if count data (non-negative integers)
  if (all(resp_non_missing >= 0) &&
      all(resp_non_missing == round(resp_non_missing))) {
    return("count")
  }

  # Otherwise treat as continuous data
  return("continuous")
}

#' @title pECV for Poisson factor model with missing data
#' @description pECV implementation for count data with missing values
#' @noRd
pECV.miss_poisson <- function(resp, miss_ind, n, p, C, qmax, fold, tol_val, theta0, A0) {

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
    # Create train/test masks considering missing data
    test_mask <- split_matrices[[k]]
    train_mask <- Reduce("+", split_matrices[-k]) > 0

    # Apply missing data mask
    test_mask <- test_mask * miss_ind == 1
    train_mask <- train_mask * miss_ind == 1

    # Create training set
    train_set <- matrix(0, nrow = nrow(resp), ncol = ncol(resp))
    train_set[train_mask] <- resp[train_mask]

    num_train <- sum(train_mask == 1)
    if (num_train == 0) {
      loss_values[k, ] <- NA
      next
    }
    tol_current <- tol_val / num_train

    # Fit models with different numbers of factors
    for (q_est in 1:qmax) {
      tryCatch({
        jml_res <- confirm_CJMLE_poisson_cpp(
          train_set, train_mask,
          theta0[, 1:(q_est + 1)], A0[, 1:(q_est + 1)],
          matrix(TRUE, p, q_est + 1),
          C = C, tol = tol_current
        )

        jml_result <- JL_poisson(resp, jml_res$theta, jml_res$A, test_mask)
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
    p1ECV = which.min(cv_loss + ((K * K) / (K - 1)) * (1:qmax) * pen1),
    p2ECV = which.min(cv_loss + ((K * K) / (K - 1)) * (1:qmax) * pen2),
    p3ECV = which.min(cv_loss + (K / (K - 1)) * (1:qmax) * pen1),
    p4ECV = which.min(cv_loss + (K / (K - 1)) * (1:qmax) * pen2),
    ECV_loss = cv_loss
  )

  return(result)
}

#' @title pECV for Gaussian factor model with missing data
#' @description pECV implementation for continuous data with missing values
#' @noRd
pECV.miss_gaussian <- function(resp, miss_ind, n, p, C, qmax, fold, tol_val, theta0, A0) {

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
    # Create train/test masks considering missing data
    test_mask <- split_matrices[[k]]
    train_mask <- Reduce("+", split_matrices[-k]) > 0

    # Apply missing data mask
    test_mask <- test_mask * miss_ind == 1
    train_mask <- train_mask * miss_ind == 1

    # Create training set
    train_set <- matrix(0, nrow = nrow(resp), ncol = ncol(resp))
    train_set[train_mask] <- resp[train_mask]

    num_train <- sum(train_mask == 1)
    if (num_train == 0) {
      loss_values[k, ] <- NA
      next
    }
    tol_current <- tol_val / num_train

    # Fit models with different numbers of factors
    for (q_est in 1:qmax) {
      tryCatch({
        jml_res <- CJMLE_linear(
          train_set, train_mask,
          theta0[, 1:(q_est + 1)], A0[, 1:(q_est + 1)],
          matrix(TRUE, p, q_est + 1), 1,
          C = C, tol = tol_current, FALSE
        )

        jml_result <- JL_gaussian(resp, jml_res$theta, jml_res$A, test_mask)
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
    p1ECV = which.min(cv_loss + ((K * K) / (K - 1)) * (1:qmax) * pen1),
    p2ECV = which.min(cv_loss + ((K * K) / (K - 1)) * (1:qmax) * pen2),
    p3ECV = which.min(cv_loss + (K / (K - 1)) * (1:qmax) * pen1),
    p4ECV = which.min(cv_loss + (K / (K - 1)) * (1:qmax) * pen2),
    ECV_loss = cv_loss
  )

  return(result)
}

#' @title pECV for Binary factor model with missing data
#' @description pECV implementation for binary data with missing values
#' @noRd
pECV.miss_binary <- function(resp, miss_ind, n, p, C, qmax, fold, tol_val, theta0, A0) {

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
    # Create train/test masks considering missing data
    test_mask <- split_matrices[[k]]
    train_mask <- Reduce("+", split_matrices[-k]) > 0

    # Apply missing data mask
    test_mask <- test_mask * miss_ind == 1
    train_mask <- train_mask * miss_ind == 1

    # Create training and test sets with NA for missing values
    train_set <- matrix(NA, nrow = nrow(resp), ncol = ncol(resp))
    train_set[train_mask] <- resp[train_mask]
    test_set <- matrix(NA, nrow = nrow(resp), ncol = ncol(resp))
    test_set[test_mask] <- resp[test_mask]

    # Check if there's enough data
    if (sum(train_mask) < q_est * 2) {
      loss_values[k, ] <- NA
      next
    }

    # Fit models with different numbers of factors
    for (q_est in 1:qmax) {
      tryCatch({
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

# ===== Data generation functions with missing values =====


#' @title Introduce missing values to a matrix
#' @description Helper function to introduce random missing values
#' @param mat Data matrix
#' @param miss_prop Proportion of missing values
#' @noRd
introduce_missing <- function(mat, miss_prop) {
  n <- nrow(mat)
  p <- ncol(mat)

  # Completely random missing
  miss_idx <- sample(1:(n*p), size = floor(miss_prop * n * p))
  mat[miss_idx] <- NA

  return(mat)
}

#' @title Generate continuous data with missing values
#' @description Generate simulated data based on Gaussian factor model with missing values
#' @param n Number of observations
#' @param p Number of variables
#' @param q True number of factors
#' @param noise_sd Standard deviation of noise
#' @param miss_prop Proportion of missing values (default 0.05)
#' @export
generate_continuous_data_miss <- function(n = 100, p = 50, q = 3, noise_sd = 1,
                                          miss_prop = 0.05) {
  # Generate complete data
  data_complete <- generate_continuous_data(n, p, q, noise_sd)
  resp <- data_complete$resp

  # Introduce missing values
  resp_miss <- introduce_missing(resp, miss_prop)

  return(list(resp = resp_miss,
              resp_complete = resp,
              true_q = q,
              theta_true = data_complete$theta_true,
              A_true = data_complete$A_true,
              miss_prop = miss_prop))
}

#' @title Generate count data with missing values
#' @description Generate simulated data based on Poisson factor model with missing values
#' @param n Number of observations
#' @param p Number of variables
#' @param q True number of factors
#' @param miss_prop Proportion of missing values (default 0.05)
#' @export
generate_count_data_miss <- function(n = 100, p = 50, q = 3,
                                     miss_prop = 0.05) {
  # Generate complete data
  data_complete <- generate_count_data(n, p, q)
  resp <- data_complete$resp

  # Introduce missing values
  resp_miss <- introduce_missing(resp, miss_prop)

  return(list(resp = resp_miss,
              resp_complete = resp,
              true_q = q,
              theta_true = data_complete$theta_true,
              A_true = data_complete$A_true,
              miss_prop = miss_prop))
}

#' @title Generate binary data with missing values
#' @description Generate simulated data based on Binary factor model with missing values
#' @param n Number of observations
#' @param p Number of variables
#' @param q True number of factors
#' @param miss_prop Proportion of missing values (default 0.05)
#' @export
generate_binary_data_miss <- function(n = 100, p = 50, q = 3,
                                      miss_prop = 0.05) {
  # Generate complete data
  data_complete <- generate_binary_data(n, p, q)
  resp <- data_complete$resp

  # Introduce missing values
  resp_miss <- introduce_missing(resp, miss_prop)

  return(list(resp = resp_miss,
              resp_complete = resp,
              true_q = q,
              theta_true = data_complete$theta_true,
              A_true = data_complete$A_true,
              d_true = data_complete$d_true,
              miss_prop = miss_prop))
}
