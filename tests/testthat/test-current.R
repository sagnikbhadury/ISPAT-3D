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

