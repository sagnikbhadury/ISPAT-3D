test_that("current covariance fitter returns full positive-definite zone matrices", {
  set.seed(19)
  x <- matrix(rnorm(80 * 4), 80, 4)
  y <- matrix(rnorm(75 * 4), 75, 4)
  colnames(x) <- colnames(y) <- c("A","B","C","D")
  fit <- ispat3d_fit_covariance(list(Low = cov(x), High = cov(y)),
                                c(80, 75), rank = 2, maxit = 30)
  expect_equal(dim(fit$shared), c(4L, 4L))
  expect_true(all(vapply(fit$full, function(z)
    min(eigen(z, symmetric = TRUE, only.values = TRUE)$values) > 0,
    logical(1))))
  expect_true(all(vapply(fit$partial, function(z)
    all(diag(z) == 1), logical(1))))
})

test_that("current Vecchia 3D and section-wise 2D entry points run", {
  skip_if_not_installed("gpboost")
  set.seed(71)
  n <- 48L
  xyz <- cbind(runif(n), runif(n), rep(0:3, each = 12) / 2)
  zone <- rep(c("Low","High"), each = 24)
  section <- rep(1:4, each = 12)
  Y <- cbind(A = sin(xyz[, 1] * 4) + rnorm(n, 0, .2),
             B = cos(xyz[, 2] * 3) + rnorm(n, 0, .2),
             C = rnorm(n))
  args <- list(Y = Y, coords = xyz, zones = zone, sections = section,
               rank = 2, anchor_min = 10, anchor_max = 20,
               neighbors = 5, gp_maxit = 2, factor_maxit = 25,
               threads = 2)
  f3 <- do.call(ispat3d_fit, args)
  f2 <- do.call(ispat3d_fit_2d, args)
  expect_identical(names(f3$full), names(f2$full))
  expect_length(f3$gp_log$status, 6L)
  expect_length(f2$gp_log$status, 12L)
  expect_true(all(vapply(f3$partial, function(z) all(is.finite(z)), logical(1))))
  expect_true(all(vapply(f2$partial, function(z) all(is.finite(z)), logical(1))))
})

test_that("example image yields matched section zones and plot-ready networks", {
  image <- ispat3d_example_image(n_per_section = 30L, n_sections = 3L)
  expect_equal(dim(image$coords), c(90L, 3L))
  expect_equal(dim(image$kde), c(90L, 3L))
  expect_true(all(table(image$zones, image$sections) == 15L))
  expect_true(all(is.finite(image$Y)))

  pcor <- matrix(c(1, 0.3, -0.2, 0.3, 1, 0,
                   -0.2, 0, 1), 3L, 3L,
                 dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  edges <- ispat3d_edge_table(pcor, threshold = 0.1)
  expect_equal(nrow(edges), 2L)
  expect_equal(edges$sign, c("positive", "negative"))
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(ispat3d_plot_network(pcor, threshold = 0.1))
  fit <- list(full = list(Low = diag(3L), High = diag(3L)))
  expect_invisible(ispat3d_plot_zones(fit, threshold = 0.1))
})
