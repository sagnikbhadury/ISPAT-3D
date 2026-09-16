#' Extract signed conditional-association edges
#'
#' Uses partial correlations from a full fitted zone covariance. A low-rank
#' shared factor covariance alone is not invertible and is not a network input.
#' @param x An ISPAT3D fit returned by ispat3d_fit() or a partial-correlation
#'   matrix.
#' @param zone Zone name when x is a fit.
#' @param threshold Minimum absolute partial correlation to retain.
#' @return A data frame with source, target, effect, sign, and magnitude,
#'   sorted by decreasing absolute effect.
#' @export
ispat3d_edge_table <- function(x, zone = NULL, threshold = 0) {
  if (is.list(x) && !is.null(x$partial)) {
    if (is.null(zone) || !zone %in% names(x$partial))
      stop("supply a zone present in fit$partial")
    x <- x$partial[[zone]]
  }
  x <- as.matrix(x)
  if (nrow(x) != ncol(x) || nrow(x) < 2L || any(!is.finite(x)))
    stop("x must be a finite square partial-correlation matrix")
  if (!isTRUE(all.equal(x, t(x), tolerance = 1e-7)))
    stop("x must be symmetric")
  if (length(threshold) != 1L || !is.finite(threshold) ||
      threshold < 0 || threshold > 1)
    stop("threshold must be between 0 and 1")
  labels <- colnames(x)
  if (is.null(labels)) labels <- paste0("Variable_", seq_len(ncol(x)))
  loc <- which(upper.tri(x) & abs(x) > threshold, arr.ind = TRUE)
  if (!nrow(loc)) return(data.frame(from = character(), to = character(),
    partial_correlation = numeric(), sign = character(),
    absolute_effect = numeric()))
  effects <- x[loc]
  out <- data.frame(from = labels[loc[, 1]], to = labels[loc[, 2]],
    partial_correlation = effects,
    sign = ifelse(effects >= 0, "positive", "negative"),
    absolute_effect = abs(effects))
  out[order(-out$absolute_effect, out$from, out$to), , drop = FALSE]
}

.ispat3d_partial_for_plot <- function(x, zone) {
  if (is.list(x) && !is.null(x$full)) {
    if (is.null(zone) || !zone %in% names(x$full))
      stop("supply a zone present in fit$full")
    return(ispat3d_partial_correlation(x$full[[zone]]))
  }
  if (!is.null(zone)) stop("zone is only used when x is an ISPAT3D fit")
  as.matrix(x)
}

#' Plot a circular partial-correlation network using base graphics
#'
#' The fitted zone covariance is transformed to partial correlations using
#' the full shared, zone-specific, and uniqueness terms. Red and blue edges
#' show positive and negative conditional density associations.
#' @param x An ISPAT3D fit or partial-correlation matrix.
#' @param zone Zone name when x is a fit.
#' @param threshold Minimum absolute partial correlation to display.
#' @param main Plot title.
#' @param positive,negative Edge colors for positive and negative associations.
#' @param vertex_cex Node size multiplier.
#' @param label_cex Label size multiplier.
#' @param edge_scale Controls edge widths.
#' @return Invisibly returns the displayed edge table.
#' @export
ispat3d_plot_network <- function(x, zone = NULL, threshold = 0.05,
                                 main = zone, positive = "#B2182B",
                                 negative = "#2166AC", vertex_cex = 1.1,
                                 label_cex = 0.75, edge_scale = 4) {
  pcor <- .ispat3d_partial_for_plot(x, zone)
  edges <- ispat3d_edge_table(pcor, threshold = threshold)
  g <- ncol(pcor)
  labels <- colnames(pcor)
  if (is.null(labels)) labels <- paste0("Variable_", seq_len(g))
  angle <- pi / 2 - 2 * pi * (seq_len(g) - 1) / g
  xx <- cos(angle); yy <- sin(angle)
  graphics::plot(NA_real_, NA_real_, type = "n", xlim = c(-1.55, 1.55),
    ylim = c(-1.4, 1.4), xlab = "", ylab = "", axes = FALSE, asp = 1,
    main = main)
  if (nrow(edges)) {
    max_effect <- max(edges$absolute_effect)
    for (k in seq_len(nrow(edges))) {
      a <- match(edges$from[k], labels)
      b <- match(edges$to[k], labels)
      graphics::segments(xx[a], yy[a], xx[b], yy[b],
        col = if (edges$partial_correlation[k] >= 0) positive else negative,
        lwd = 0.5 + edge_scale * edges$absolute_effect[k] / max_effect)
    }
  }
  graphics::points(xx, yy, pch = 21, bg = "white", col = "#333333",
    cex = vertex_cex)
  graphics::text(1.19 * xx, 1.19 * yy, labels = labels,
    cex = label_cex, adj = ifelse(xx > 0.15, 0, ifelse(xx < -0.15, 1, 0.5)),
    xpd = NA)
  invisible(edges)
}

#' Facet all fitted zone networks in one base-R figure
#' @param fit An ISPAT3D fit from ispat3d_fit() or ispat3d_fit_2d().
#' @param zones Zone names to plot, in order.
#' @param columns Number of panel columns.
#' @param ... Further arguments passed to ispat3d_plot_network().
#' @return Invisibly returns a named list of displayed edge tables.
#' @export
ispat3d_plot_zones <- function(fit, zones = names(fit$full),
                               columns = 2L, ...) {
  if (!is.list(fit) || is.null(fit$full) ||
      !length(zones) || any(!zones %in% names(fit$full)))
    stop("fit must have full zone covariances for every requested zone")
  columns <- as.integer(columns)
  if (length(columns) != 1L || is.na(columns) || columns < 1L)
    stop("columns must be a positive integer")
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(ceiling(length(zones) / columns), columns),
                mar = c(1, 1, 3, 1))
  out <- lapply(zones, function(z)
    ispat3d_plot_network(fit, zone = z, main = z, ...))
  names(out) <- zones
  invisible(out)
}
