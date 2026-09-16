#' Simulate a small registered cell map with section-wise KDE inputs
#'
#' Generates a didactic cell map with three annotated cell types, section-wise
#' Gaussian kernel-density estimates evaluated at every cell, and two relative
#' tumor-density zones per section. It is designed for package examples and
#' checks, not for biological simulation studies.
#' @param n_per_section Cells per section; an even multiple of six, at least 24.
#' @param n_sections Number of serial sections, at least two.
#' @param bandwidth Gaussian KDE bandwidth in the simulated coordinate units.
#' @param seed Reproducible simulation seed.
#' @return A list containing coordinates, sections, source cell types, KDE
#'   values, transformed model matrix Y, zones, and tumor-density score.
#' @export
ispat3d_example_image <- function(n_per_section = 24L, n_sections = 3L,
                                  bandwidth = 0.12, seed = 2026L) {
  n_per_section <- as.integer(n_per_section)
  n_sections <- as.integer(n_sections)
  if (length(n_per_section) != 1L || is.na(n_per_section) ||
      n_per_section < 24L || n_per_section %% 6L != 0L)
    stop("n_per_section must be an even multiple of six, at least 24")
  if (length(n_sections) != 1L || is.na(n_sections) || n_sections < 2L)
    stop("n_sections must be at least two")
  if (length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0)
    stop("bandwidth must be positive")
  set.seed(seed)
  n <- n_per_section * n_sections
  section <- rep(seq_len(n_sections), each = n_per_section)
  types <- c("Tumor", "T_cell", "Macrophage")
  source_type <- unlist(lapply(seq_len(n_sections), function(i)
    sample(rep(types, each = n_per_section / 3L))), use.names = FALSE)
  centers <- c(Tumor = 0.70, T_cell = 0.30, Macrophage = 0.50)
  x <- pmin(0.99, pmax(0.01,
    stats::rnorm(n, centers[source_type] + 0.03 * (section - 1L), 0.18)))
  y <- stats::runif(n)
  z <- 0.15 * (section - 1L)
  coords <- cbind(x = x, y = y, z = z)
  kde <- matrix(0, n, length(types), dimnames = list(NULL, types))
  for (s in seq_len(n_sections)) {
    loc <- which(section == s)
    for (g in seq_along(types)) {
      source <- loc[source_type[loc] == types[g]]
      dx <- outer(x[loc], x[source], "-")
      dy <- outer(y[loc], y[source], "-")
      kde[loc, g] <- rowMeans(exp(-(dx^2 + dy^2) /
        (2 * bandwidth^2))) / (2 * pi * bandwidth^2)
    }
  }
  tumor_score <- kde[, "Tumor"]
  high <- logical(n)
  for (s in seq_len(n_sections)) {
    loc <- which(section == s)
    high[loc] <- rank(tumor_score[loc], ties.method = "first") >
      n_per_section / 2L
  }
  zones <- factor(ifelse(high, "High", "Low"), levels = c("Low", "High"))
  list(coords = coords, sections = section, cell_type = source_type,
       kde = kde, Y = log1p(1e9 * kde), zones = zones,
       tumor_score = tumor_score)
}
