# Current volumetric and matched section-wise planar ISPAT workflow.

.ispat3d_zone_names <- function(zones) {
  values <- unique(as.character(zones))
  if (is.factor(zones)) return(levels(zones)[levels(zones) %in% values])
  if (is.numeric(zones)) return(as.character(sort(unique(zones))))
  conventional <- c("Very Low", "Low", "Intermediate", "High", "Very High")
  if (all(values %in% conventional))
    return(conventional[conventional %in% values])
  values
}

.ispat3d_pipeline <- function(Y, coords, zones, sections, rank, anchor_fraction,
                              anchor_min, anchor_max, neighbors, gp_maxit,
                              factor_maxit, seed, threads, return_residuals,
                              planar) {
  Y <- as.matrix(Y)
  coords <- as.matrix(coords)
  storage.mode(Y) <- storage.mode(coords) <- "double"
  n <- nrow(Y); g <- ncol(Y)
  if (n < 2L || g < 2L || nrow(coords) != n ||
      ncol(coords) != 3L || length(zones) != n ||
      length(sections) != n || anyNA(zones) || anyNA(sections) ||
      any(!is.finite(Y)) || any(!is.finite(coords)))
    stop("Y (cells by variables), coords (cells by x,y,z), zones, and sections must have matching finite rows")
  if (is.null(colnames(Y))) colnames(Y) <- paste0("Variable_", seq_len(g))
  if (anyDuplicated(colnames(Y))) stop("Y column labels must be unique")
  if (!is.finite(anchor_fraction) || anchor_fraction <= 0 ||
      !is.finite(anchor_min) || anchor_min < 1 ||
      !is.finite(anchor_max) || anchor_max < anchor_min)
    stop("invalid anchor selection parameters")
  zone_names <- .ispat3d_zone_names(zones)
  if (length(zone_names) < 2L) stop("at least two zones are required")
  if (planar && anyNA(sections)) stop("sections are required for the planar fit")
  set.seed(seed)
  covariances <- residuals <- vector("list", length(zone_names))
  names(covariances) <- names(residuals) <- zone_names
  counts <- numeric(length(zone_names))
  log <- list()
  entry <- 0L
  for (q in seq_along(zone_names)) {
    idx <- which(as.character(zones) == zone_names[q])
    counts[q] <- length(idx)
    if (length(idx) <= g) stop("zone ", zone_names[q], " has too few cells")
    R <- matrix(0, length(idx), g, dimnames = list(NULL, colnames(Y)))
    if (planar) {
      groups <- split(seq_along(idx), sections[idx])
    } else {
      groups <- list(all = seq_along(idx))
    }
    for (h in seq_along(groups)) {
      loc <- groups[[h]]
      rows <- idx[loc]
      xyz <- if (planar) coords[rows, 1:2, drop = FALSE] else coords[rows, , drop = FALSE]
      anchor_n <- min(length(rows), as.integer(min(anchor_max,
        max(anchor_min, ceiling(anchor_fraction * length(rows))))))
      anchor <- if (length(rows) <= anchor_n) seq_along(rows) else
        .ispat3d_spatial_pick(seq_along(rows), xyz[, 1], xyz[, 2],
          if (planar) rep(1L, length(rows)) else sections[rows], anchor_n)
      for (j in seq_len(g)) {
        y <- Y[rows, j]
        fitted <- if (planar) {
          tryCatch(.ispat3d_gp_adjust(y, xyz, anchor,
            seed = seed + 1000L * q + 100L * h + j,
            neighbors = neighbors, maxit = gp_maxit, threads = threads),
            error = function(e) list(residual = y - mean(y),
              status = paste0("fallback:", conditionMessage(e)),
              fit_seconds = NA_real_, predict_seconds = NA_real_,
              parameters = rep(NA_real_, ncol(xyz) + 2L)))
        } else .ispat3d_gp_adjust(y, xyz, anchor,
          seed = seed + 1000L * q + j, neighbors = neighbors,
          maxit = gp_maxit, threads = threads)
        R[loc, j] <- fitted$residual
        entry <- entry + 1L
        log[[entry]] <- data.frame(zone = zone_names[q],
          section = if (planar) as.character(names(groups)[h]) else "pooled",
          variable = colnames(Y)[j], cells = length(rows),
          anchors = length(anchor), status = fitted$status,
          fit_seconds = fitted$fit_seconds,
          predict_seconds = fitted$predict_seconds)
      }
    }
    covariances[[q]] <- stats::cov(R)
    if (return_residuals) residuals[[q]] <- R
  }
  active <- which(vapply(seq_len(g), function(j)
    max(vapply(covariances, function(x) x[j, j], numeric(1))) > 1e-10,
    logical(1)))
  if (length(active) < 2L) stop("fewer than two nonconstant variables")
  reduced <- lapply(covariances, function(x) x[active, active, drop = FALSE])
  fitted <- ispat3d_fit_covariance(reduced, counts,
    rank = min(as.integer(rank), length(active) - 1L), maxit = factor_maxit)
  embed <- function(x, diagonal = FALSE) {
    out <- matrix(0, g, g, dimnames = list(colnames(Y), colnames(Y)))
    out[active, active] <- x
    if (diagonal) diag(out)[-active] <- 1
    out
  }
  fitted$shared <- embed(fitted$shared)
  fitted$full <- lapply(fitted$full, embed)
  fitted$partial <- lapply(fitted$partial, embed, diagonal = TRUE)
  fitted$active_variables <- colnames(Y)[active]
  fitted$inactive_variables <- colnames(Y)[-active]
  fitted$gp_log <- do.call(rbind, log)
  fitted$covariances <- covariances
  fitted$method <- if (planar) "section-wise planar Vecchia GP + covariance likelihood" else
    "pooled volumetric Vecchia GP + covariance likelihood"
  if (return_residuals) fitted$residuals <- residuals
  fitted
}

#' Fit current ISPAT-3D on selected cells
#'
#' Fits a MatÃ©rn-3/2 anisotropic GP per variable and zone using a 15-neighbor
#' Vecchia likelihood on spatially balanced anchors; predicts at all selected
#' cells; then fits shared-plus-zone covariance from complete residual
#' covariance summaries by Gaussian maximum likelihood.
#'
#' @param Y Numeric cells-by-variables matrix, usually log1p(1e9 * KDE).
#' @param coords Numeric cells-by-3 matrix in registered x,y,z coordinates.
#' @param zones Zone label for each row.
#' @param sections Section identifier for balanced anchor selection.
#' @param rank Shared and zone-specific factor rank (default 5).
#' @param anchor_fraction Fraction of cells selected as GP anchors.
#' @param anchor_min,anchor_max Minimum and maximum number of anchors per fit.
#' @param neighbors Number of Vecchia neighbors.
#' @param gp_maxit Maximum GP likelihood iterations.
#' @param factor_maxit Maximum covariance likelihood iterations.
#' @param seed Random seed for anchor selection and GP fits.
#' @param threads GPBoost threads per GP fit.
#' @param return_residuals Whether to return full adjusted residual matrices.
#' @return List containing shared covariance, full zone covariances, zone
#' partial correlations, diagnostics, and optionally residuals.
#' @export
ispat3d_fit <- function(Y, coords, zones, sections = rep(1L, nrow(Y)),
                        rank = 5L, anchor_fraction = 0.10,
                        anchor_min = 300L, anchor_max = 5000L,
                        neighbors = 15L, gp_maxit = 30L,
                        factor_maxit = 350L, seed = 2026L, threads = 4L,
                        return_residuals = FALSE) {
  .ispat3d_pipeline(Y, coords, zones, sections, rank, anchor_fraction,
    anchor_min, anchor_max, neighbors, gp_maxit, factor_maxit, seed,
    threads, return_residuals, planar = FALSE)
}

#' Fit the matched section-wise ISPAT-2D comparator
#'
#' Uses the same selected rows and covariance estimator as ispat3d_fit(), but
#' fits independent planar GPs within each zone-section group. Groups with
#' fewer than 10 cells are mean-centered; failed planar GP fits are recorded
#' and mean-centered. This changes both spatial dimension and section pooling.
#' @inheritParams ispat3d_fit
#' @param sections Required section identifier for every selected cell.
#' @return Same structure as ispat3d_fit().
#' @export
ispat3d_fit_2d <- function(Y, coords, zones, sections,
                           rank = 5L, anchor_fraction = 0.10,
                           anchor_min = 300L, anchor_max = 5000L,
                           neighbors = 15L, gp_maxit = 30L,
                           factor_maxit = 350L, seed = 2026L, threads = 2L,
                           return_residuals = FALSE) {
  if (missing(sections)) stop("sections must be supplied for the planar fit")
  .ispat3d_pipeline(Y, coords, zones, sections, rank, anchor_fraction,
    anchor_min, anchor_max, neighbors, gp_maxit, factor_maxit, seed,
    threads, return_residuals, planar = TRUE)
}

