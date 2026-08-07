#' Run ISPat-3D for a single spot
#'
#' Fits the 3D spatial GP (via \code{\link{ext_loop_3d}}) to every marker
#' within every zone of a single spot/field of view, then runs a multi-study
#' factor analysis (\code{\link{cavi_msfa}} or \code{\link{svi_msfa}}) across
#' zones on the resulting de-spatialised residuals to recover a shared
#' cell-type interaction network and zone-specific networks.
#'
#' @param Y_sp G x N marker/gene expression matrix for this spot (G markers,
#'   N cells).
#' @param S_sp N x 4 matrix of (x, y, z, zone) for this spot; the 4th column
#'   is an integer zone/density annotation.
#' @param ncores Number of cores for the parallel per-marker GP fit.
#' @param Kernel Spatial covariance kernel: \code{"Matern"} or \code{"RBF"}.
#' @param MSFA_method Multi-study factor analysis algorithm: \code{"CAVI"} or
#'   \code{"SVI"}.
#'
#' @return A list with elements \code{Shared_Net} (G x G shared factor
#'   covariance component), \code{Zone_Nets} (named list of G x G full
#'   zone-specific covariance matrices),
#'   \code{Z_est} (de-spatialised residuals per zone), and \code{VBfit} (raw
#'   MSFA fit).
#'
#' @seealso \code{\link{ISPAT_3D}}, which calls this function once per spot;
#'   \code{\link{run_ispat_spot_bigmem}} for a memory-efficient variant.
#' @export
run_ispat_spot <- function(Y_sp, S_sp, ncores,
                            Kernel = "Matern",
                            MSFA_method = "CAVI") {

  G <- nrow(Y_sp)
  Clusters <- sort(unique(S_sp[, 4]))
  C <- length(Clusters)
  zone_names <- c("Very Low", "Low", "Intermediate", "High", "Very High")

  N_c <- numeric(C)
  Z_est <- vector("list", C)

  for (c in seq_along(Clusters)) {
    pos_c <- which(S_sp[, 4] == Clusters[c])
    N_c[c] <- length(pos_c)
    S_loc_c <- S_sp[pos_c, 1:3]

    cl_par <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl_par)

    Z_raw <- suppressWarnings({
      foreach::foreach(my_count = 1:G,
                        .packages = c("Matrix", "matrixcalc"),
                        .combine = cbind,
                        .export = c("ext_loop_3d", ".build_kernel", ".gp_neglml")) %dopar% {
        ext_loop_3d(my_count, N_c[c], pos_c, Y_sp, S_loc_c, Kernel)
      }
    })

    parallel::stopCluster(cl_par)

    # Ensure MSFA receives a proper numeric matrix (N_c x G)
    Z_est[[c]] <- as.matrix(Z_raw)
  }

  # Factor counts: ceiling(2*log(G)) -- robust across the G=15,20,25 range
  K_shared <- ceiling(2 * log(G))
  K_zone <- rep(ceiling(2 * log(G)), C)

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

#' ISPat-3D: spatial network estimation for 3D multiplexed tissue imaging
#'
#' Runs \code{\link{run_ispat_spot}} across every spot/field of view in a 3D
#' (or serial-section) multiplexed imaging dataset, estimating a shared and a
#' zone-specific cell-type interaction network for each spot. GP
#' hyperparameters (length scales, spatial/error variances, intercept) are
#' estimated per marker per zone via type-II marginal likelihood optimization
#' (L-BFGS-B); de-spatialised residuals are then passed to a multi-study
#' factor analysis to recover the interaction networks.
#'
#' @param Y G x M marker/gene expression matrix across all spots (G markers,
#'   M total cells).
#' @param S M x 4 matrix of (x, y, z, zone) coordinates/annotations across
#'   all spots.
#' @param spots_vec Length-M vector giving the spot/field-of-view identifier
#'   for each cell/location in \code{Y} and \code{S}.
#' @param ncores Number of cores for the parallel per-marker GP fit. Defaults
#'   to two fewer than the number of available cores.
#' @param Kernel Spatial covariance kernel: \code{"Matern"} or \code{"RBF"}.
#' @param MSFA_method Multi-study factor analysis algorithm: \code{"CAVI"} or
#'   \code{"SVI"}.
#'
#' @return A named list (one element per spot) of the result from
#'   \code{\link{run_ispat_spot}}.
#'
#' @examples
#' \dontrun{
#' result <- ISPAT_3D(Y, S, spots_vec, ncores = 8, Kernel = "Matern", MSFA_method = "CAVI")
#' result[["spot_1"]]$Shared_Net
#' }
#'
#' @export
ISPAT_3D <- function(Y, S, spots_vec,
                      ncores = max(1, parallel::detectCores() - 2),
                      Kernel = c("Matern", "RBF"),
                      MSFA_method = c("CAVI", "SVI")) {

  Kernel <- match.arg(Kernel)
  MSFA_method <- match.arg(MSFA_method)

  spots_uniq <- unique(spots_vec)
  n_spots <- length(spots_uniq)
  G <- nrow(Y)

  message("ISPat-3D")
  message("PCs (G): ", G, " | Total cells: ", ncol(Y), " | Spots: ", n_spots)
  message("Kernel: ", Kernel, " | MSFA: ", MSFA_method)

  results <- vector("list", n_spots)
  names(results) <- spots_uniq

  for (sp_idx in seq_along(spots_uniq)) {
    sp <- spots_uniq[sp_idx]
    sp_cells <- which(spots_vec == sp)

    message(sprintf("[%d/%d] Spot: %s | Cells: %d", sp_idx, n_spots, sp, length(sp_cells)))

    results[[sp]] <- run_ispat_spot(
      Y_sp = Y[, sp_cells, drop = FALSE],
      S_sp = S[sp_cells, , drop = FALSE],
      ncores = ncores,
      Kernel = Kernel,
      MSFA_method = MSFA_method
    )

    message(sprintf("  Spot %s complete.", sp))
  }

  message("ISPat-3D complete. Results available for ", n_spots, " spots.")
  results
}
