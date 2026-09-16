# Gaussian-likelihood shared-plus-zone covariance fit used in the manuscript.

.ispat3d_positive_factors <- function(x, rank) {
  eig <- eigen((x + t(x)) / 2, symmetric = TRUE)
  i <- seq_len(min(rank, ncol(x)))
  sweep(eig$vectors[, i, drop = FALSE], 2, sqrt(pmax(eig$values[i], 0)), "*")
}

.ispat3d_softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

#' Convert a full covariance matrix to partial correlations
#' @param covariance Symmetric positive-definite covariance matrix.
#' @param ridge Small numerical diagonal ridge (default 1e-6).
#' @return A signed partial-correlation matrix.
#' @export
ispat3d_partial_correlation <- function(covariance, ridge = 1e-6) {
  covariance <- as.matrix(covariance)
  if (nrow(covariance) != ncol(covariance) || any(!is.finite(covariance)))
    stop("covariance must be a finite square matrix")
  theta <- solve(covariance + diag(ridge, nrow(covariance)))
  out <- -theta / sqrt(outer(diag(theta), diag(theta)))
  diag(out) <- 1
  dimnames(out) <- dimnames(covariance)
  out
}

#' Fit shared and zone-specific factor covariance from sufficient statistics
#'
#' Fits Sigma_q = Phi Phi' + Lambda_q Lambda_q' + diag(psi_q) by a weighted
#' Gaussian covariance likelihood. This is the current estimator used after
#' Vecchia GP adjustment, not a variational factor posterior.
#' @param covariances Named list of complete sample covariance matrices.
#' @param counts Number of adjusted cells in each zone, in list order.
#' @param rank Shared and zone-specific loading rank (default 5).
#' @param maxit Maximum L-BFGS-B iterations.
#' @return Shared covariance, full zone covariances, partial correlations,
#' loadings, uniquenesses, and optimization diagnostics.
#' @export
ispat3d_fit_covariance <- function(covariances, counts, rank = 5L, maxit = 350L) {
  if (!is.list(covariances) || length(covariances) < 2L)
    stop("covariances must contain at least two zones")
  qn <- length(covariances)
  counts <- as.numeric(counts)
  g <- ncol(as.matrix(covariances[[1L]]))
  if (length(counts) != qn || any(!is.finite(counts)) || any(counts <= g))
    stop("counts must match zones and exceed the number of variables")
  rank <- as.integer(rank)
  if (length(rank) != 1L || is.na(rank) || rank < 1L || rank >= g)
    stop("rank must be between 1 and G-1")
  mats <- lapply(covariances, as.matrix)
  labels <- colnames(mats[[1L]])
  for (x in mats) {
    if (!identical(dim(x), c(g, g)) || any(!is.finite(x)))
      stop("all covariances must be finite G-by-G matrices")
  }
  scales <- sqrt(pmax(Reduce("+", lapply(mats, diag)) / qn, 1e-8))
  corrs <- lapply(mats, function(x) sweep(sweep(x, 1L, scales, "/"), 2L, scales, "/"))
  weights <- counts / sum(counts)
  pooled <- Reduce("+", Map(function(x, w) x * w, corrs, weights))
  phi0 <- .ispat3d_positive_factors(pooled * 0.65, rank)
  lam0 <- lapply(corrs, function(x) .ispat3d_positive_factors(
    (x - tcrossprod(phi0)) * 0.55, rank))
  psi0 <- Map(function(x, lam) pmax(diag(x) - rowSums(phi0^2) -
                                    rowSums(lam^2), 0.05), corrs, lam0)
  start <- c(as.vector(phi0), unlist(lapply(lam0, as.vector)),
             unlist(lapply(psi0, log)))
  n_phi <- g * rank
  n_lam <- qn * n_phi
  unpack <- function(par) {
    phi <- matrix(par[seq_len(n_phi)], g, rank)
    lam <- lapply(seq_len(qn), function(q) {
      first <- n_phi + (q - 1L) * n_phi + 1L
      matrix(par[first:(first + n_phi - 1L)], g, rank)
    })
    eta <- matrix(utils::tail(par, qn * g), g, qn)
    list(phi = phi, lam = lam, eta = eta)
  }
  evaluate <- function(par) {
    b <- unpack(par)
    shared <- tcrossprod(b$phi)
    value <- 0
    grad_phi <- matrix(0, g, rank)
    grad_lam <- vector("list", qn)
    grad_eta <- matrix(0, g, qn)
    for (q in seq_len(qn)) {
      psi <- .ispat3d_softplus(b$eta[, q]) + 1e-5
      sigma <- shared + tcrossprod(b$lam[[q]]) + diag(psi)
      ch <- tryCatch(chol(sigma), error = function(e) NULL)
      if (is.null(ch)) stop("non-positive covariance during optimization")
      inv <- chol2inv(ch)
      value <- value + weights[q] * (2 * sum(log(diag(ch))) +
                                     sum(diag(inv %*% corrs[[q]]))) / 2
      d <- weights[q] * (inv - inv %*% corrs[[q]] %*% inv) / 2
      grad_phi <- grad_phi + 2 * d %*% b$phi
      grad_lam[[q]] <- 2 * d %*% b$lam[[q]]
      grad_eta[, q] <- diag(d) * stats::plogis(b$eta[, q])
    }
    value <- value + 0.004 * mean(b$phi^2) +
      0.012 * mean(unlist(b$lam)^2)
    grad_phi <- grad_phi + 0.008 * b$phi / n_phi
    grad_lam <- lapply(grad_lam, function(x) x)
    for (q in seq_len(qn))
      grad_lam[[q]] <- grad_lam[[q]] + 0.024 * b$lam[[q]] / n_lam
    list(value = value, gradient = c(as.vector(grad_phi),
      unlist(lapply(grad_lam, as.vector)), as.vector(grad_eta)))
  }
  cache <- new.env(parent = emptyenv())
  cached <- function(par) {
    if (is.null(cache$par) || !identical(par, cache$par)) {
      cache$par <- par
      cache$result <- evaluate(par)
    }
    cache$result
  }
  fit <- stats::optim(start, fn = function(par) cached(par)$value,
               gr = function(par) cached(par)$gradient,
               method = "L-BFGS-B", control = list(maxit = as.integer(maxit),
               factr = 1e7))
  b <- unpack(fit$par)
  dimnames(b$phi) <- list(labels, paste0("shared_", seq_len(rank)))
  for (q in seq_len(qn)) rownames(b$lam[[q]]) <- labels
  zone_names <- names(covariances)
  if (is.null(zone_names)) zone_names <- paste0("Zone_", seq_len(qn))
  names(b$lam) <- zone_names
  shared <- tcrossprod(b$phi)
  full <- lapply(seq_len(qn), function(q) {
    psi <- .ispat3d_softplus(b$eta[, q]) + 1e-5
    x <- shared + tcrossprod(b$lam[[q]]) + diag(psi)
    x * outer(scales, scales)
  })
  shared <- shared * outer(scales, scales)
  dimnames(shared) <- list(labels, labels)
  full <- lapply(full, function(x) {dimnames(x) <- list(labels, labels); x})
  names(full) <- zone_names
  partial <- lapply(full, ispat3d_partial_correlation)
  uniqueness <- lapply(seq_len(qn), function(q)
    (.ispat3d_softplus(b$eta[, q]) + 1e-5) * scales^2)
  names(uniqueness) <- zone_names
  list(shared = shared, full = full, partial = partial,
       loadings_shared = sweep(b$phi, 1L, scales, "*"),
       loadings_zone = lapply(b$lam, function(x) sweep(x, 1L, scales, "*")),
       uniqueness = uniqueness, counts = stats::setNames(counts, zone_names),
       objective = fit$value, convergence = fit$convergence,
       message = fit$message, iterations = fit$counts)
}

