#' @title Estimate constraint constant C for continuous data
#' @description Data-driven estimation of the constraint constant C in alternating maximization algorithm for continuous data
#'              using truncated SVD approach. This function decomposes the data matrix
#'              and estimates C based on the maximum row norms.
#' @param X n x p continuous data matrix
#' @param qmax Rank for truncated SVD (default 8)
#' @param safety Safety parameter for conservative estimation (default 1.2)
#' @return A list containing:
#'   \item{qmax}{Truncation rank used}
#'   \item{safety}{Safety parameter applied}
#'   \item{C_norm_hat}{Original maximum row norm}
#'   \item{C_est}{Final conservative estimate of C}
#'   \item{a_norms}{Row norms of factor matrix A}
#'   \item{b_norms}{Row norms of factor matrix B}
#' @details
#' The function performs the following steps:
#' 1. Computes truncated SVD of X with rank qmax
#' 2. Constructs factor matrices A = U * sqrt(D) and B = V * sqrt(D)
#' 3. Calculates row 2-norms for matrices A and B
#' 4. Takes the maximum norm and multiplies by safety parameter
#'
#' For count data, it is recommended to transform the data using log(X + 1) before
#' applying this function.
#' @export
#' @examples
#' # Example 1: Continuous data
#' set.seed(123)
#' n <- 100; p <- 50; q <- 3
#' theta_true <- matrix(runif(n * q), n, q)
#' A_true <- matrix(runif(p * q), p, q)
#' X <- theta_true %*% t(A_true) + matrix(rnorm(n * p, sd = 0.5), n, p)
#'
#' # Estimate C
#' C_result <- estimate_C(X, qmax = 5)
#' print(C_result$C_est)
#'
#' # Example 2: Count data (apply log transformation)
#' lambda <- exp(theta_true %*% t(A_true))
#' X_count <- matrix(rpois(n * p, lambda = as.vector(lambda)), n, p)
#' X_transformed <- log(X_count + 1)
#' C_count <- estimate_C(X_transformed, qmax = 5)
#' print(C_count$C_est)

estimate_C <- function(X, qmax = 8, safety = 1.2) {
  # X      : n×p continuous data matrix
  # qmax   : rank for truncated SVD (default 8)
  # safety : safety parameter (default 1.2)

  # Check for irlba package
  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Please install 'irlba' package: install.packages('irlba')")
  }

  # Check input
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (any(is.na(X))) {
    warning("Data contains NA values. Consider using pECV.miss() for missing data.")
  }

  n <- nrow(X)
  p <- ncol(X)

  # Ensure qmax is not larger than matrix dimensions
  qmax <- min(qmax, n - 1, p - 1)

  # 1) Truncated SVD
  svd_res <- irlba::irlba(X, nv = qmax, nu = qmax)
  U <- svd_res$u           # n×qmax
  V <- svd_res$v           # p×qmax
  d <- svd_res$d[1:qmax]   # first qmax singular values

  # 2) Construct factor matrices A_hat, B_hat
  D_sqrt <- diag(sqrt(d), nrow = qmax, ncol = qmax)
  A_hat  <- U %*% D_sqrt   # n×qmax
  B_hat  <- V %*% D_sqrt   # p×qmax

  # 3) Calculate row 2-norms
  a_norms <- sqrt(rowSums(A_hat^2))  # length n
  b_norms <- sqrt(rowSums(B_hat^2))  # length p

  # 4) Take maximum norm and multiply by safety factor
  C_norm_hat <- max(c(a_norms, b_norms))
  C_est      <- C_norm_hat * safety

  # 5) Return results
  list(
    qmax        = qmax,
    safety      = safety,
    C_norm_hat  = C_norm_hat,  # original max row-norm
    C_est       = C_est,        # final conservative estimate
    a_norms     = a_norms,
    b_norms     = b_norms
  )
}

#' @title Estimate constraint constant C for binary data
#' @description Data-driven estimation of the constraint constant C for binary data
#'              using cross-window smoothing and empirical logit transformation.
#' @param X n x p binary data matrix (0/1 values)
#' @param qmax Rank for truncated SVD (default 8)
#' @param safety Safety parameter for conservative estimation (default 1.5)
#' @param eps Small constant to avoid logit divergence when p=0 or p=1 (default 1e-12)
#' @param radius Radius for cross-window smoothing (default 1)
#' @return A list containing:
#'   \item{radius}{Cross-window radius used}
#'   \item{qmax}{Truncation rank used}
#'   \item{safety}{Safety parameter applied}
#'   \item{C0}{Original maximum row norm}
#'   \item{C_est}{Final conservative estimate of C}
#'   \item{a_norms}{Row norms of factor matrix A}
#'   \item{b_norms}{Row norms of factor matrix B}
#'   \item{Mhat}{Logit-transformed matrix}
#'   \item{P_smooth}{Smoothed probability matrix}
#'   \item{N_counts}{Count of values in each smoothing window}
#' @details
#' The function performs the following steps:
#' 1. Applies cross-window smoothing to estimate probabilities
#' 2. Performs empirical logit transformation with smoothing
#' 3. Computes truncated SVD of the transformed matrix
#' 4. Constructs matrices A and B and calculates row norms
#' 5. Estimates C as the maximum norm times safety parameter
#'
#' The cross-window smoothing helps stabilize probability estimates,
#' especially for sparse binary data.
#' @export

estimate_C_binary <- function(X,
                              qmax   = 8,      # truncation rank
                              safety = 1.5,    # safety amplification factor
                              eps    = 1e-12,  # avoid p=0/1 causing logit divergence
                              radius = 1       # cross-window smoothing radius
) {

  # Check for irlba package
  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Please install 'irlba' package: install.packages('irlba')")
  }

  # Check input
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  # Verify binary data
  unique_vals <- unique(as.vector(X[!is.na(X)]))
  if (!all(unique_vals %in% c(0, 1))) {
    warning("Data contains non-binary values. Function expects 0/1 values.")
  }

  n <- nrow(X)
  p <- ncol(X)

  # Ensure qmax is reasonable
  qmax <- min(qmax, n - 1, p - 1)

  # 1) Cross-window smoothing
  P_hat    <- matrix(NA_real_, n, p)
  N_counts <- matrix(0L,       n, p)

  for (i in seq_len(n)) {
    for (j in seq_len(p)) {
      # Collect all values in cross window
      window_values <- c(X[i, j])

      # Vertical direction
      for (d in seq_len(radius)) {
        if (i - d >= 1)       window_values <- c(window_values, X[i - d, j])
        if (i + d <= n)       window_values <- c(window_values, X[i + d, j])
      }

      # Horizontal direction
      for (d in seq_len(radius)) {
        if (j - d >= 1)       window_values <- c(window_values, X[i, j - d])
        if (j + d <= p)       window_values <- c(window_values, X[i, j + d])
      }

      # Smoothed proportion & effective count
      vals <- window_values[!is.na(window_values)]
      P_hat[i, j]    <- mean(vals)
      N_counts[i, j] <- length(vals)
    }
  }

  # 2) Empirical logit transformation (with smoothing eps)
  P_adj <- (P_hat * N_counts + eps) / (N_counts + 5*eps)
  Mhat  <- log(P_adj / (1 - P_adj))

  # 3) Truncated SVD → factor decomposition
  svd_res <- irlba::irlba(Mhat, nv = qmax, nu = qmax)
  U <- svd_res$u     # n×qmax
  V <- svd_res$v     # p×qmax
  d <- svd_res$d[1:qmax]
  Dsqrt <- diag(sqrt(d), qmax, qmax)
  A_hat <- U %*% Dsqrt  # n×qmax
  B_hat <- V %*% Dsqrt  # p×qmax

  # 4) Calculate row norms
  a_norms <- sqrt(rowSums(A_hat^2))
  b_norms <- sqrt(rowSums(B_hat^2))

  # 5) Take maximum norm, then multiply by safety factor
  C0    <- max(c(a_norms, b_norms))
  C_est <- C0 * safety

  list(
    radius   = radius,
    qmax     = qmax,
    safety   = safety,
    C0       = C0,
    C_est    = C_est,
    a_norms  = a_norms,
    b_norms  = b_norms,
    Mhat     = Mhat,
    P_smooth = P_hat,
    N_counts = N_counts
  )
}

# ===== Comprehensive Examples =====

#' @examples
#' \dontrun{
#' # Load ECV package
#' library(ECV)
#'
#' # ===== Example 1: Complete workflow for continuous data =====
#' set.seed(123)
#'
#' # Generate data
#' continuous_data <- generate_continuous_data(n = 100, p = 50, q = 3)
#'
#' # Estimate C
#' C_est_result <- estimate_C(continuous_data$resp, qmax = 6)
#' cat("Estimated C for continuous data:", C_est_result$C_est, "\n")
#'
#' # Use estimated C in pECV
#' result_cont <- pECV(
#'   resp = continuous_data$resp,
#'   C = C_est_result$C_est,  # Use estimated C
#'   qmax = 6,
#'   fold = 5,
#'   seed = 1
#' )
#' cat("Selected factors with estimated C:", result_cont$ECV, "\n")
#'
#' # ===== Example 2: Complete workflow for count data =====
#' count_data <- generate_count_data(n = 50, p = 50, q = 2)
#'
#' # Transform count data
#' count_transformed <- log(count_data$resp + 1)
#'
#' # Estimate C
#' C_est_count <- estimate_C(count_transformed, qmax = 4)
#' cat("Estimated C for count data:", C_est_count$C_est, "\n")
#'
#' # Use in pECV
#' result_count <- pECV(
#'   resp = count_data$resp,
#'   C = C_est_count$C_est,
#'   qmax = 4,
#'   fold = 5,
#'   seed = 1,
#'   data_type = "count"
#' )
#' cat("Selected factors with estimated C:", result_count$ECV, "\n")
#'
