# GPBoost Vecchia spatial adjustment used by the current analysis.

.ispat3d_gp_adjust <- function(y, coords, anchor, seed, neighbors = 15L,
                               maxit = 30L, threads = 2L) {
  n <- length(y)
  if (n < 10L || stats::sd(y) < 1e-9)
    return(list(residual = y - mean(y), status = if (n < 10L) "small_section" else
                  "constant", fit_seconds = 0, predict_seconds = 0,
                parameters = rep(NA_real_, ncol(coords) + 2L)))
  if (stats::sd(y[anchor]) < 1e-9)
    return(list(residual = y - mean(y), status = "constant_anchor",
                fit_seconds = 0, predict_seconds = 0,
                parameters = rep(NA_real_, ncol(coords) + 2L)))
  if (!requireNamespace("gpboost", quietly = TRUE))
    stop("The gpboost package is required for Vecchia GP adjustment")
  status <- "ok"
  acoords <- coords[anchor, , drop = FALSE]
  ya <- y[anchor]
  nn <- min(as.integer(neighbors), length(anchor) - 1L)
  fit_model <- function(params) gpboost::fitGPModel(
    gp_coords = acoords, cov_function = "matern_ard", cov_fct_shape = 1.5,
    gp_approx = "vecchia", num_neighbors = nn,
    num_parallel_threads = as.integer(threads), seed = as.integer(seed),
    y = ya, X = matrix(1, length(ya), 1), params = params)
  tfit <- system.time({
    model <- tryCatch(fit_model(list(maxit = as.integer(maxit),
                                  delta_rel_conv = 1e-3)),
      error = function(e) {
        status <<- "retry_init"
        v <- max(stats::var(ya), 1e-4)
        ranges <- apply(acoords, 2, function(x) max(diff(range(x)) / 3, 1e-3))
        fit_model(list(maxit = as.integer(maxit), delta_rel_conv = 1e-3,
                       init_cov_pars = c(max(0.2 * v, 1e-4),
                                         max(0.8 * v, 1e-4), ranges)))
      })
  })
  tpred <- system.time({
    pred <- stats::predict(model, gp_coords_pred = coords, X_pred = matrix(1, n, 1),
                    predict_response = FALSE)[["mu"]]
  })
  if (length(pred) != n || any(!is.finite(pred)))
    stop("Vecchia GP prediction failed")
  list(residual = as.numeric(y - pred), status = status,
       fit_seconds = unname(tfit[["elapsed"]]),
       predict_seconds = unname(tpred[["elapsed"]]),
       parameters = as.numeric(model$get_cov_pars()))
}

