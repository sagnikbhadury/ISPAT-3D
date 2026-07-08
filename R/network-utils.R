#' Convert a covariance/network matrix to partial correlations
#'
#' Utility for turning a shared or zone-specific network returned by
#' \code{\link{run_ispat_spot}} / \code{\link{ISPAT_3D}} (e.g.
#' \code{result$Shared_Net}) into a partial correlation matrix, which is
#' often easier to interpret and threshold for visualization: entries close
#' to zero indicate cell types that are conditionally independent given all
#' other cell types.
#'
#' @param Sigma A symmetric G x G covariance/network matrix.
#' @param cell_types Optional character vector of length G giving row/column
#'   names for the result; if omitted, the row/column names of \code{Sigma}
#'   are used (if present).
#'
#' @return A symmetric G x G partial correlation matrix with unit diagonal.
#'
#' @examples
#' \dontrun{
#' pcor <- cov_to_pcor(result$Shared_Net, cell_types = rownames(Y))
#' }
#'
#' @export
cov_to_pcor <- function(Sigma, cell_types = rownames(Sigma)) {
  Theta <- tryCatch(solve(Sigma + diag(1e-6, nrow(Sigma))),
                     error = function(e) MASS::ginv(Sigma))
  D <- sqrt(diag(Theta))
  Pcor <- -Theta / outer(D, D)
  diag(Pcor) <- 1
  if (!is.null(cell_types)) {
    rownames(Pcor) <- colnames(Pcor) <- cell_types
  }
  Pcor
}

#' Chord diagram of a partial correlation network
#'
#' Draws a circular ("chord") diagram of a partial correlation matrix (see
#' \code{\link{cov_to_pcor}}), with nodes arranged on a circle, edges colored
#' by the sign of the partial correlation (positive = conditionally
#' co-localised, negative = conditional spatial exclusion), and edge
#' width/opacity proportional to edge strength.
#'
#' This function requires the \pkg{dplyr}, \pkg{igraph}, \pkg{tidygraph},
#' \pkg{ggraph}, \pkg{ggplot2}, and \pkg{scales} packages, none of which are
#' installed automatically with ISPAT3D; install them separately before
#' calling this function.
#'
#' @param pcor_mat A symmetric partial correlation matrix with row/column
#'   names, as returned by \code{\link{cov_to_pcor}}.
#' @param title Plot title.
#' @param threshold Only edges with \code{abs(partial correlation) >= threshold}
#'   are drawn.
#' @param node_colors Optional named character vector of node fill colors,
#'   keyed by the row/column names of \code{pcor_mat}. If omitted, a default
#'   qualitative palette is used.
#'
#' @return A \pkg{ggplot2}/\pkg{ggraph} plot object, or \code{NULL} (invisibly,
#'   with a message) if no edges exceed \code{threshold}.
#'
#' @examples
#' \dontrun{
#' pcor <- cov_to_pcor(result$Shared_Net, cell_types = rownames(Y))
#' plot_pcor_chord(pcor, title = "Shared network", threshold = 0.05)
#' }
#'
#' @export
plot_pcor_chord <- function(pcor_mat, title, threshold, node_colors = NULL) {

  needed <- c("dplyr", "igraph", "tidygraph", "ggraph", "ggplot2", "scales")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("plot_pcor_chord() requires the following packages: ",
         paste(missing, collapse = ", "), ". Install them and try again.", call. = FALSE)
  }

  cnames <- rownames(pcor_mat)
  n_nodes <- length(cnames)

  if (is.null(node_colors)) {
    node_colors <- stats::setNames(
      grDevices::hcl.colors(n_nodes, palette = "Dark 3"), cnames
    )
  }

  edge_df <- expand.grid(from = cnames, to = cnames, stringsAsFactors = FALSE)
  edge_df <- dplyr::filter(edge_df, from < to)
  edge_df <- dplyr::mutate(edge_df,
    pcor = mapply(function(f, t) pcor_mat[f, t], from, to),
    abs_pcor = abs(pcor),
    direction = ifelse(pcor >= 0, "Positive", "Negative")
  )
  edge_df <- dplyr::filter(edge_df, abs_pcor >= threshold)
  edge_df <- dplyr::arrange(edge_df, dplyr::desc(abs_pcor))

  if (nrow(edge_df) == 0) {
    message(sprintf("[%s] No edges above threshold %.3f", title, threshold))
    return(invisible(NULL))
  }

  node_df <- data.frame(name = cnames, stringsAsFactors = FALSE)

  g <- igraph::graph_from_data_frame(
    d = edge_df[, c("from", "to", "pcor", "abs_pcor", "direction")],
    directed = FALSE,
    vertices = node_df
  )

  node_order <- igraph::V(g)$name
  angles <- seq(0, 2 * pi, length.out = n_nodes + 1)[-(n_nodes + 1)]
  names(angles) <- cnames
  node_angles <- angles[node_order]

  tg <- tidygraph::as_tbl_graph(g)

  ggraph::ggraph(tg, layout = "linear", circular = TRUE) +
    ggraph::geom_edge_arc(
      ggplot2::aes(edge_colour = direction, edge_width = abs_pcor, edge_alpha = abs_pcor),
      strength = 0.45
    ) +
    ggraph::geom_node_point(
      size = 9, shape = 21, fill = node_colors[node_order], color = "white", stroke = 1.4
    ) +
    ggraph::geom_node_label(
      ggplot2::aes(label = name),
      size = 3.0, fontface = "bold",
      fill = scales::alpha("white", 0.85), label.size = 0,
      nudge_x = cos(node_angles) * 0.28,
      nudge_y = sin(node_angles) * 0.28
    ) +
    ggraph::scale_edge_color_manual(
      values = c("Positive" = "#1565C0", "Negative" = "#B71C1C"),
      name = "Partial correlation"
    ) +
    ggraph::scale_edge_width_continuous(range = c(0.4, 3.5), guide = "none") +
    ggraph::scale_edge_alpha_continuous(range = c(0.3, 0.95), guide = "none") +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("%d edges  |  threshold |pcor| >= %.3f", nrow(edge_df), threshold),
      caption = "Edge width proportional to |partial correlation|"
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5,
                                          margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5, color = "grey40"),
      plot.caption = ggplot2::element_text(size = 8, hjust = 0.5, color = "grey55"),
      plot.margin = ggplot2::margin(20, 20, 20, 20),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 9)
    )
}
