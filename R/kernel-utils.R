#' Anisotropic 3D spatial covariance kernel
#'
#' Internal helper building an N x N spatial covariance/correlation matrix
#' over 3D coordinates (x, y, z), using separate length scales for the
#' in-plane (x, y) and depth (z) directions, so that spatial smoothness can
#' differ across tissue sections versus within a section.
#'
#' @param S_loc_c N x 3 matrix of (x, y, z) coordinates.
#' @param lS In-plane (x, y) length scale.
#' @param lZ Depth (z) length scale.
#' @param Kernel \code{"Matern"} (Matern 3/2) or any other value for a
#'   squared-exponential kernel.
#'
#' @return An N x N kernel matrix with unit diagonal.
#' @keywords internal
#' @noRd
.build_kernel <- function(S_loc_c, lS, lZ, Kernel) {
  N <- nrow(S_loc_c)
  KS <- matrix(0, N, N)
  for (i in 1:N) {
    for (j in 1:i) {
      dx <- S_loc_c[i, 1] - S_loc_c[j, 1]
      dy <- S_loc_c[i, 2] - S_loc_c[j, 2]
      dz <- S_loc_c[i, 3] - S_loc_c[j, 3]
      if (Kernel == "Matern") {
        r <- sqrt((dx / lS)^2 + (dy / lS)^2 + (dz / lZ)^2)
        # Guard against Inf*0 = NaN when r is very large (short lengthscale)
        KS[i, j] <- if (r > 500) 0 else (1 + sqrt(3) * r) * exp(-sqrt(3) * r)
      } else {
        r2 <- (dx / lS)^2 + (dy / lS)^2 + (dz / lZ)^2
        KS[i, j] <- if (r2 > 1e6) 0 else exp(-0.5 * r2)
      }
      KS[j, i] <- KS[i, j]
    }
    KS[i, i] <- 1
  }
  KS
}

#' Negative log marginal likelihood for the 3D spatial GP
#'
#' Internal helper: negative log marginal likelihood of a single marker's
#' expression under a Gaussian process with the anisotropic 3D kernel from
#' \code{\link{.build_kernel}}, used by \code{\link[stats]{optim}} (L-BFGS-B)
#' to estimate GP hyperparameters per marker per zone via type-II marginal
#' likelihood.
#'
#' @param theta Numeric vector
#'   \code{c(log_lS, log_lZ, log_sigmaS, log_sigmaEPS, beta)}.
#' @param y Numeric vector of marker expression at each location.
#' @param S_loc_c N x 3 matrix of (x, y, z) coordinates.
#' @param Kernel \code{"Matern"} or \code{"RBF"}, passed to \code{\link{.build_kernel}}.
#'
#' @return The negative log marginal likelihood (a single number); a large
#'   penalty value (1e10) if the covariance is not numerically invertible.
#' @keywords internal
#' @noRd
.gp_neglml <- function(theta, y, S_loc_c, Kernel) {
  lS <- exp(theta[1])
  lZ <- exp(theta[2])
  sigmaS <- exp(theta[3])
  sigmaEPS <- exp(theta[4])
  beta <- theta[5]
  N <- length(y)

  KS <- .build_kernel(S_loc_c, lS, lZ, Kernel)
  V <- sigmaS^2 * KS + sigmaEPS^2 * diag(N)

  cholV <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(cholV)) {
    V <- V + 1e-6 * diag(N)
    cholV <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(cholV)) return(1e10)
  }

  resid <- y - beta
  alpha <- backsolve(cholV, forwardsolve(t(cholV), resid))
  0.5 * (sum(resid * alpha) + 2 * sum(log(diag(cholV))) + N * log(2 * pi))
}

#' Fit the 3D spatial GP for one marker within one zone
#'
#' Internal worker called (in parallel) once per marker by
#' \code{\link{run_ispat_spot}}. Estimates GP hyperparameters
#' (length scales, spatial/error variances, intercept) for one marker via
#' type-II marginal likelihood optimization, then returns the spatially
#' adjusted residual \code{z_hat = (y - beta_hat) - B_hat} to be passed on to
#' the multi-study factor analysis step.
#'
#' @param count Marker/gene row index into \code{Y}.
#' @param N_c Number of locations in the current zone.
#' @param pos Indices (into \code{Y}'s columns) of the locations in this zone.
#' @param Y Marker/gene expression matrix.
#' @param S_loc_c N_c x 3 matrix of (x, y, z) coordinates for this zone.
#' @param Kernel \code{"Matern"} or \code{"RBF"}.
#'
#' @return A one-column data frame of the spatially adjusted residual for
#'   this marker.
#' @keywords internal
#' @noRd
ext_loop_3d <- function(count, N_c, pos, Y, S_loc_c, Kernel) {

  y <- as.vector(as.matrix(Y)[count, pos])

  xy_dists <- as.vector(dist(S_loc_c[, 1:2])); xy_dists <- xy_dists[xy_dists > 0]
  z_dists <- as.vector(dist(S_loc_c[, 3, drop = FALSE])); z_dists <- z_dists[z_dists > 0]
  lS0 <- if (length(xy_dists) > 0) median(xy_dists) else 0.1
  lZ0 <- if (length(z_dists) > 0) median(z_dists) else lS0

  opt <- tryCatch(
    optim(
      par = c(log(lS0), log(lZ0), log(1), log(0.5), 0),
      fn = .gp_neglml,
      y = y, S_loc_c = S_loc_c, Kernel = Kernel,
      method = "L-BFGS-B",
      lower = c(log(lS0 / 10), log(lZ0 / 10), log(0.01), log(0.01), -10),
      upper = c(log(lS0 * 10), log(lZ0 * 10), log(10), log(10), 10),
      control = list(maxit = 200, factr = 1e9)
    ),
    error = function(e) {
      # Fallback to median-heuristic values if optimizer fails
      list(par = c(log(lS0), log(lZ0), log(1), log(0.5), mean(y)))
    }
  )

  lS_est <- exp(opt$par[1])
  lZ_est <- exp(opt$par[2])
  sigmaS_est <- exp(opt$par[3])
  sigmaEPS_est <- exp(opt$par[4])
  beta_est <- opt$par[5]

  KS_est <- .build_kernel(S_loc_c, lS_est, lZ_est, Kernel)
  if (!matrixcalc::is.positive.definite(KS_est)) {
    KS_est <- as.matrix(Matrix::nearPD(KS_est)$mat)
  }

  sS <- max(sigmaS_est, 0.011)
  sEP <- max(sigmaEPS_est, 0.011)

  # Posterior mean of spatial nuisance: B_hat = sigmaS^2 * K * V^{-1} * (y - beta)
  # Spatially adjusted residual fed to MSFA: z_hat = (y - beta) - B_hat
  #   = sigmaEPS^2 * V^{-1} * (y - beta)
  V <- sS^2 * KS_est + sEP^2 * diag(N_c)
  cholV <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(cholV)) cholV <- chol(V + 1e-6 * diag(N_c))
  alpha <- backsolve(cholV, forwardsolve(t(cholV), y - beta_est))
  B_hat <- sS^2 * KS_est %*% alpha
  storeZ <- as.vector((y - beta_est) - B_hat)

  data.frame(storeZ)
}
