# Spatially balanced row selection shared by volumetric and planar workflows.
.ispat3d_allocation <- function(counts, budget) {
  if (budget >= sum(counts)) return(as.integer(counts))
  raw <- as.numeric(budget) * counts / sum(counts)
  out <- floor(raw)
  missing <- budget - sum(out)
  if (missing > 0L) out[order(raw - out, decreasing = TRUE)[seq_len(missing)]] <-
    out[order(raw - out, decreasing = TRUE)[seq_len(missing)]] + 1L
  as.integer(out)
}

.ispat3d_spatial_pick <- function(index, x, y, section, budget) {
  if (length(index) <= budget) return(index)
  groups <- split(index, section[index])
  nsec <- .ispat3d_allocation(lengths(groups), budget)
  out <- integer()
  for (k in seq_along(groups)) {
    members <- groups[[k]]
    take <- nsec[k]
    if (!take) next
    pool <- if (length(members) > 3L * take) sample(members, 3L * take) else members
    xx <- x[pool]; yy <- y[pool]
    bx <- pmin(10L, 1L + floor(10 * (xx - min(xx)) / (max(xx) - min(xx) + 1e-9)))
    by <- pmin(10L, 1L + floor(10 * (yy - min(yy)) / (max(yy) - min(yy) + 1e-9)))
    bins <- split(pool, paste(bx, by, sep = "_"))
    first <- vapply(bins, function(v) sample(v, 1L), integer(1))
    if (length(first) > take) first <- sample(first, take)
    if (length(first) < take) {
      left <- setdiff(pool, first)
      first <- c(first, sample(left, min(take - length(first), length(left))))
    }
    out <- c(out, first)
  }
  if (length(out) > budget) out <- sample(out, budget)
  if (length(out) < budget) out <- c(out, sample(setdiff(index, out), budget - length(out)))
  out
}

#' Select spatially balanced cells by zone and section
#' @param coords N-by-3 registered coordinate matrix.
#' @param zones Zone label per row.
#' @param sections Section identifier per row.
#' @param budget Number per zone when budget_kind="per_zone"; total otherwise.
#' @param budget_kind One of "per_zone" or "total".
#' @param seed Sampling seed.
#' @return Named list of source row indices for each zone.
#' @export
ispat3d_sample <- function(coords, zones, sections, budget,
                           budget_kind = c("per_zone", "total"), seed = 2026L) {
  budget_kind <- match.arg(budget_kind)
  coords <- as.matrix(coords)
  if (ncol(coords) != 3L || any(!is.finite(coords)) ||
      length(zones) != nrow(coords) || length(sections) != nrow(coords) ||
      anyNA(zones) || anyNA(sections))
    stop("coords, zones, and sections must have matching complete rows")
  budget <- as.integer(budget)
  if (length(budget) != 1L || is.na(budget) || budget < 1L)
    stop("budget must be a positive integer")
  zone_names <- unique(as.character(zones))
  counts <- vapply(zone_names, function(z) sum(as.character(zones) == z), integer(1))
  target <- if (budget_kind == "total") .ispat3d_allocation(counts, budget) else
    rep(budget, length(counts))
  if (any(target > counts)) stop("budget exceeds available cells in a zone")
  set.seed(seed)
  out <- lapply(seq_along(zone_names), function(q) {
    idx <- which(as.character(zones) == zone_names[q])
    .ispat3d_spatial_pick(idx, coords[, 1], coords[, 2], sections, target[q])
  })
  names(out) <- zone_names
  out
}

