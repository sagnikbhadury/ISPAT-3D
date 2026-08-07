#' Per-marker GP worker reading from a shared bigmemory matrix
#'
#' Internal worker used by \code{\link{run_ispat_spot_bigmem}}. Identical to
#' \code{\link{ext_loop_3d}}, except that it attaches to a memory-mapped
#' \pkg{bigmemory} matrix and reads only its own marker row, instead of
#' receiving a copy of the full expression matrix in each parallel worker.
#'
#' @param count Marker/gene row index.
#' @param N_c Number of locations in the current zone.
#' @param pos Indices (into the bigmemory matrix's columns) of the locations
#'   in this zone.
#' @param bm_desc A \pkg{bigmemory} descriptor (from
#'   \code{bigmemory::describe()}) for the shared expression matrix.
#' @param S_loc_c N_c x 3 matrix of (x, y, z) coordinates for this zone.
#' @param Kernel \code{"Matern"} or \code{"RBF"}.
#'
#' @return A one-column data frame of the spatially adjusted residual for
#'   this marker.
#' @keywords internal
#' @noRd
ext_loop_3d_bigmem <- function(count, N_c, pos, bm_desc, S_loc_c, Kernel) {

  # Attach to shared memory -- zero-copy read of one row
  Y_bm <- bigmemory::attach.big.matrix(bm_desc)
  y <- as.vector(Y_bm[count, pos])

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

  V <- sS^2 * KS_est + sEP^2 * diag(N_c)
  cholV <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(cholV)) cholV <- chol(V + 1e-6 * diag(N_c))
  alpha <- backsolve(cholV, forwardsolve(t(cholV), y - beta_est))
  B_hat <- sS^2 * KS_est %*% alpha
  storeZ <- as.vector((y - beta_est) - B_hat)

  data.frame(storeZ)
}

#' Memory-efficient variant of \code{run_ispat_spot} for large spots
#'
#' Drop-in replacement for \code{\link{run_ispat_spot}} that avoids
#' serializing the full marker expression matrix to every parallel worker.
#' Instead, the expression matrix is written once to a memory-mapped
#' \pkg{bigmemory} backing file, and each worker attaches to the shared file
#' and reads only its own marker row. This cuts per-worker data transfer from
#' roughly gigabytes to kilobytes, which matters when a spot has many cells
#' and/or many markers.
#'
#' Requires the \pkg{bigmemory} package, which is not installed automatically
#' with ISPAT3D; install it separately (\code{install.packages("bigmemory")})
#' before calling this function.
#'
#' @inheritParams run_ispat_spot
#' @param bm_dir Directory to use for the temporary bigmemory backing file;
#'   defaults to \code{tempdir()}. The backing file is removed when the
#'   function finishes.
#'
#' @return Same structure as \code{\link{run_ispat_spot}}.
#' @seealso \code{\link{run_ispat_spot}}
#' @export
run_ispat_spot_bigmem <- function(Y_sp, S_sp, ncores,
                                   Kernel = "Matern",
                                   MSFA_method = "CAVI",
                                   bm_dir = tempdir()) {

  if (!requireNamespace("bigmemory", quietly = TRUE)) {
    stop("Package 'bigmemory' is required for run_ispat_spot_bigmem(). ",
         "Install it with install.packages(\"bigmemory\").", call. = FALSE)
  }

  G <- nrow(Y_sp)
  Clusters <- sort(unique(S_sp[, 4]))
  C <- length(Clusters)
  zone_names <- c("Very Low", "Low", "Intermediate", "High", "Very High")

  # --- Write Y_sp once to a bigmemory backing file ---
  # Workers attach via descriptor -- no serialization of Y to each worker
  bm_file <- file.path(bm_dir, "Y_bigmem.bin")
  bm_desc_file <- file.path(bm_dir, "Y_bigmem.desc")

  message("  Writing Y to bigmemory backing file ...")
  Y_bm <- bigmemory::as.big.matrix(
    x = Y_sp,
    type = "double",
    backingfile = "Y_bigmem.bin",
    backingpath = bm_dir,
    descriptorfile = "Y_bigmem.desc"
  )
  bm_desc <- bigmemory::describe(Y_bm)
  message(sprintf("  Y written: %d x %d", nrow(Y_bm), ncol(Y_bm)))

  N_c <- numeric(C)
  Z_est <- vector("list", C)

  for (c in seq_along(Clusters)) {

    pos_c <- which(S_sp[, 4] == Clusters[c])
    N_c[c] <- length(pos_c)
    S_loc_c <- S_sp[pos_c, 1:3, drop = FALSE]

    message(sprintf("  Zone %d (%s): %d cells, %d markers ...",
                    Clusters[c], zone_names[Clusters[c]], N_c[c], G))

    cl_par <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl_par)

    Z_raw <- suppressWarnings(
      foreach::foreach(my_count = 1:G,
                        .packages = c("bigmemory", "Matrix", "matrixcalc"),
                        .combine = cbind,
                        .export = c("ext_loop_3d_bigmem", ".build_kernel", ".gp_neglml"),
                        .noexport = c("Y_sp", "Y_bm")) %dopar% {
        ext_loop_3d_bigmem(my_count, N_c[c], pos_c, bm_desc, S_loc_c, Kernel)
      }
    )

    parallel::stopCluster(cl_par)

    Z_est[[c]] <- as.matrix(Z_raw)
    message(sprintf("    Zone %d GP stage done. Z_est dim: %d x %d",
                    Clusters[c], nrow(Z_est[[c]]), ncol(Z_est[[c]])))
  }

  # --- Clean up backing files ---
  if (file.exists(bm_file)) file.remove(bm_file)
  if (file.exists(bm_desc_file)) file.remove(bm_desc_file)

  K_shared <- ceiling(2 * log(G))
  K_zone <- rep(ceiling(2 * log(G)), C)

  message("  Running MSFA ...")
  if (MSFA_method == "CAVI") {
    VBfit <- cavi_msfa(Z_est, K_shared, K_zone, scale = FALSE)
  } else {
    VBfit <- svi_msfa(Z_est, K_shared, K_zone, scale = FALSE)
  }

  Shared_Net <- tcrossprod(VBfit$mean_phi)
  Zone_Nets <- lapply(1:C, function(c) {
    tcrossprod(VBfit$mean_phi) +
      tcrossprod(VBfit$mean_lambda_s[[c]]) +
      diag(as.numeric(VBfit$mean_psi_s[[c]]))
  })
  names(Zone_Nets) <- zone_names[1:C]

  list(
    Shared_Net = Shared_Net,
    Zone_Nets = Zone_Nets,
    Z_est = Z_est,
    VBfit = VBfit
  )
}
