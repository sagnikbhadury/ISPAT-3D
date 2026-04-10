#------------------------------------------------------------------------------------
# ISPat-3D: 3D Spatial Interaction Pattern Analysis
#
# GP hyperparameters (lS, lZ, sigmaS, sigmaEPS, beta) are estimated per marker
# per zone via type-II marginal likelihood optimization (L-BFGS-B).
# De-spatialised residuals are passed to MSFA to recover zone-specific networks.
#------------------------------------------------------------------------------------

pacman::p_load(doParallel, foreach, Matrix, matrixcalc, parallel)

source("svi_msfa.R")
source("cavi_msfa.R")

#------------------------------------------------------------------------------------
# .build_kernel: anisotropic 3D kernel matrix
#------------------------------------------------------------------------------------

.build_kernel <- function(S_loc_c, lS, lZ, Kernel) {
  N  <- nrow(S_loc_c)
  KS <- matrix(0, N, N)
  for (i in 1:N) {
    for (j in 1:i) {
      dx <- S_loc_c[i,1] - S_loc_c[j,1]
      dy <- S_loc_c[i,2] - S_loc_c[j,2]
      dz <- S_loc_c[i,3] - S_loc_c[j,3]
      if (Kernel == "Matern") {
        r       <- sqrt((dx/lS)^2 + (dy/lS)^2 + (dz/lZ)^2)
        # Guard against Inf*0 = NaN when r is very large (short lengthscale)
        KS[i,j] <- if (r > 500) 0 else (1 + sqrt(3)*r) * exp(-sqrt(3)*r)
      } else {
        r2      <- (dx/lS)^2 + (dy/lS)^2 + (dz/lZ)^2
        KS[i,j] <- if (r2 > 1e6) 0 else exp(-0.5 * r2)
      }
      KS[j,i] <- KS[i,j]
    }
    KS[i,i] <- 1
  }
  KS
}

#------------------------------------------------------------------------------------
# .gp_neglml: negative log marginal likelihood
# theta = c(log_lS, log_lZ, log_sigmaS, log_sigmaEPS, beta)
#------------------------------------------------------------------------------------

.gp_neglml <- function(theta, y, S_loc_c, Kernel) {
  lS       <- exp(theta[1])
  lZ       <- exp(theta[2])
  sigmaS   <- exp(theta[3])
  sigmaEPS <- exp(theta[4])
  beta     <- theta[5]
  N        <- length(y)

  KS <- .build_kernel(S_loc_c, lS, lZ, Kernel)
  V  <- sigmaS^2 * KS + sigmaEPS^2 * diag(N)

  cholV <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(cholV)) {
    V     <- V + 1e-6 * diag(N)
    cholV <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(cholV)) return(1e10)
  }

  resid <- y - beta
  alpha <- backsolve(cholV, forwardsolve(t(cholV), resid))
  0.5 * (sum(resid * alpha) + 2*sum(log(diag(cholV))) + N*log(2*pi))
}

#------------------------------------------------------------------------------------
# ext_loop_3d: fit GP for one marker within one zone; return spatially adjusted
# residual z_hat = (y - beta_hat) - B_hat for input to MSFA
#------------------------------------------------------------------------------------

ext_loop_3d <- function(count, N_c, pos, Y, S_loc_c, Kernel) {

  y <- as.vector(as.matrix(Y)[count, pos])

  xy_dists <- as.vector(dist(S_loc_c[, 1:2])); xy_dists <- xy_dists[xy_dists > 0]
  z_dists  <- as.vector(dist(S_loc_c[, 3, drop=FALSE])); z_dists <- z_dists[z_dists > 0]
  lS0 <- if (length(xy_dists) > 0) median(xy_dists) else 0.1
  lZ0 <- if (length(z_dists)  > 0) median(z_dists)  else lS0

  opt <- tryCatch(
    optim(
      par     = c(log(lS0), log(lZ0), log(1), log(0.5), 0),
      fn      = .gp_neglml,
      y       = y, S_loc_c = S_loc_c, Kernel = Kernel,
      method  = "L-BFGS-B",
      lower   = c(log(lS0/10), log(lZ0/10), log(0.01), log(0.01), -10),
      upper   = c(log(lS0*10), log(lZ0*10), log(10),   log(10),    10),
      control = list(maxit = 200, factr = 1e9)
    ),
    error = function(e) {
      # Fallback to median-heuristic values if optimizer fails
      list(par = c(log(lS0), log(lZ0), log(1), log(0.5), mean(y)))
    }
  )

  lS_est       <- exp(opt$par[1])
  lZ_est       <- exp(opt$par[2])
  sigmaS_est   <- exp(opt$par[3])
  sigmaEPS_est <- exp(opt$par[4])
  beta_est     <- opt$par[5]

  KS_est <- .build_kernel(S_loc_c, lS_est, lZ_est, Kernel)
  if (!matrixcalc::is.positive.definite(KS_est)) {
    KS_est <- as.matrix(Matrix::nearPD(KS_est)$mat)
  }

  sS  <- max(sigmaS_est,   0.011)
  sEP <- max(sigmaEPS_est, 0.011)

  # Posterior mean of spatial nuisance: B_hat = sigmaS^2 * K * V^{-1} * (y - beta)
  # Spatially adjusted residual fed to MSFA: z_hat = (y - beta) - B_hat
  #   = sigmaEPS^2 * V^{-1} * (y - beta)
  V     <- sS^2 * KS_est + sEP^2 * diag(N_c)
  cholV <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(cholV)) cholV <- chol(V + 1e-6 * diag(N_c))
  alpha  <- backsolve(cholV, forwardsolve(t(cholV), y - beta_est))
  B_hat  <- sS^2 * KS_est %*% alpha
  storeZ <- as.vector((y - beta_est) - B_hat)

  return(data.frame(storeZ))
}

#------------------------------------------------------------------------------------
# run_ispat_spot: run ISPat-3D for a single spot
#------------------------------------------------------------------------------------

run_ispat_spot <- function(Y_sp, S_sp, ncores,
                           Kernel      = "Matern",
                           MSFA_method = "CAVI") {

  G          <- nrow(Y_sp)
  Clusters   <- sort(unique(S_sp[, 4]))
  C          <- length(Clusters)
  zone_names <- c("Very Low", "Low", "Intermediate", "High", "Very High")

  N_c   <- numeric(C)
  Z_est <- vector("list", C)

  for (c in seq_along(Clusters)) {
    pos_c   <- which(S_sp[, 4] == Clusters[c])
    N_c[c]  <- length(pos_c)
    S_loc_c <- S_sp[pos_c, 1:3]

    cl_par <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl_par)

    Z_raw <- suppressWarnings({
      foreach(my_count = 1:G,
              .packages = c("Matrix", "matrixcalc"),
              .combine  = cbind,
              .export   = c("ext_loop_3d", ".build_kernel", ".gp_neglml")) %dopar% {
        ext_loop_3d(my_count, N_c[c], pos_c, Y_sp, S_loc_c, Kernel)
      }
    })

    parallel::stopCluster(cl_par)

    # Ensure MSFA receives a proper numeric matrix (N_c x G)
    Z_est[[c]] <- as.matrix(Z_raw)
  }

  # Factor counts: ceiling(2*log(G)) — robust across the G=15,20,25 range
  K_shared <- ceiling(2 * log(G))
  K_zone   <- rep(ceiling(2 * log(G)), C)

  if (MSFA_method == "CAVI") {
    VBfit <- cavi_msfa(Z_est, K_shared, K_zone, scale = FALSE)
  } else {
    VBfit <- svi_msfa(Z_est, K_shared, K_zone, scale = FALSE)
  }

  Shared_Net <- tcrossprod(VBfit$mean_phi)
  Zone_Nets  <- lapply(1:C, function(c) {
    tcrossprod(VBfit$mean_lambda_s[[c]]) + tcrossprod(VBfit$mean_phi)
  })
  names(Zone_Nets) <- zone_names[1:C]

  return(list(
    Shared_Net = Shared_Net,
    Zone_Nets  = Zone_Nets,
    Z_est      = Z_est,
    VBfit      = VBfit
  ))
}

#------------------------------------------------------------------------------------
# ISPAT_3D: outer loop over all spots
#------------------------------------------------------------------------------------

ISPAT_3D <- function(Y, S, spots_vec,
                     ncores      = max(1, parallel::detectCores() - 2),
                     Kernel      = c("Matern", "RBF"),
                     MSFA_method = c("CAVI", "SVI")) {

  Kernel      <- match.arg(Kernel)
  MSFA_method <- match.arg(MSFA_method)

  spots_uniq <- unique(spots_vec)
  n_spots    <- length(spots_uniq)
  G          <- nrow(Y)

  cat("ISPat-3D\n")
  cat("PCs (G):", G, "| Total cells:", ncol(Y), "| Spots:", n_spots, "\n")
  cat("Kernel:", Kernel, "| MSFA:", MSFA_method, "\n\n")

  results        <- vector("list", n_spots)
  names(results) <- spots_uniq

  for (sp_idx in seq_along(spots_uniq)) {
    sp       <- spots_uniq[sp_idx]
    sp_cells <- which(spots_vec == sp)

    cat(sprintf("[%d/%d] Spot: %s | Cells: %d\n",
                sp_idx, n_spots, sp, length(sp_cells)))

    results[[sp]] <- run_ispat_spot(
      Y_sp        = Y[, sp_cells, drop = FALSE],
      S_sp        = S[sp_cells, , drop = FALSE],
      ncores      = ncores,
      Kernel      = Kernel,
      MSFA_method = MSFA_method
    )

    cat(sprintf("  Spot %s complete.\n\n", sp))
  }

  cat("ISPat-3D complete. Results available for", n_spots, "spots.\n")
  return(results)
}
