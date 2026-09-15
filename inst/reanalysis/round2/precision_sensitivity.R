# Sensitivity analysis requested during peer review.
# Fits sparse precision matrices directly to saved GP-adjusted residuals,
# bypassing the MSFA diagonal-uniqueness covariance decomposition.

rm(list = ls())
suppressPackageStartupMessages(library(glasso))

set.seed(2026)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Usage: precision_sensitivity.R analysis_root")
out_dir <- file.path(args[1],"precision_sensitivity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

datasets <- list(
  CRC = list(base = file.path(args[1],"crc"), ns = c(50000, 100000, 150000)),
  BC  = list(base = file.path(args[1],"bc"), ns = c(25000, 45000, 65000))
)
zones <- c("Very Low", "Low", "Intermediate", "High", "Very High")

read_mat <- function(path) as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
upper <- function(x) x[upper.tri(x)]

theta_to_pcor <- function(theta) {
  ans <- -theta / outer(sqrt(diag(theta)), sqrt(diag(theta)))
  diag(ans) <- 1
  ans
}

# EBIC-selected graphical lasso; gamma=0.5 favors reproducible sparse graphs.
fit_direct_graph <- function(z, gamma = 0.5) {
  keep <- vapply(as.data.frame(z), function(x) {
    s <- sd(x, na.rm = TRUE)
    is.finite(s) && s > 1e-10 && all(is.finite(x))
  }, logical(1))
  z <- z[, keep, drop = FALSE]
  z <- scale(z)
  s <- cov(z)
  n <- nrow(z)
  p <- ncol(z)
  grid <- exp(seq(log(0.005), log(0.5), length.out = 40))
  fits <- lapply(grid, function(rho) {
    fit <- glasso(s, rho = rho)
    theta <- fit$wi
    logdet <- as.numeric(determinant(theta, logarithm = TRUE)$modulus)
    loglik <- n * (logdet - sum(s * theta)) / 2
    edges <- sum(abs(theta[upper.tri(theta)]) > 1e-8)
    ebic <- -2 * loglik + edges * log(n) + 4 * gamma * edges * log(p)
    list(rho = rho, ebic = ebic, theta = theta)
  })
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "ebic"))]]
  pcor <- theta_to_pcor(best$theta)
  rownames(pcor) <- colnames(pcor) <- colnames(z)
  list(pcor = pcor, rho = best$rho, ebic = best$ebic)
}

fisher_mean <- function(mats) {
  zs <- lapply(mats, function(x) {
    clipped <- x
    clipped[clipped > 0.999999] <- 0.999999
    clipped[clipped < -0.999999] <- -0.999999
    atanh(clipped)
  })
  ans <- tanh(Reduce(`+`, zs) / length(zs))
  diag(ans) <- 1
  rownames(ans) <- rownames(mats[[1]])
  colnames(ans) <- colnames(mats[[1]])
  ans
}

comparison <- list()
edge_output <- list()
penalties <- list()

for (dataset in names(datasets)) {
  info <- datasets[[dataset]]
  replicate_graphs <- vector("list", length(info$ns))

  for (i in seq_along(info$ns)) {
    n <- info$ns[i]
    residuals <- readRDS(file.path(info$base, paste0("tier_", format(n, scientific = FALSE, trim = TRUE)),
                                   "gp_residuals.rds"))
    replicate_graphs[[i]] <- list()
    for (zone in zones) {
      fit <- fit_direct_graph(residuals[[zone]])
      replicate_graphs[[i]][[zone]] <- fit$pcor
      penalties[[length(penalties) + 1]] <- data.frame(
        Dataset = dataset, Zone = zone, Subsample_size = n,
        Selected_rho = fit$rho, EBIC = fit$ebic)
    }
  }

  for (zone in zones) {
    direct_mats <- lapply(replicate_graphs, `[[`, zone)
    retained <- Reduce(intersect, lapply(direct_mats, rownames))
    direct_mats <- lapply(direct_mats, function(x) x[retained, retained, drop = FALSE])
    direct <- fisher_mean(direct_mats)
    msfa <- read_mat(file.path(info$base, "combined",
                               paste0("pcor_combined_", gsub(" ", "_", zone), ".csv")))
    common <- intersect(rownames(direct), rownames(msfa))
    direct <- direct[common, common, drop = FALSE]
    msfa <- msfa[common, common, drop = FALSE]
    x <- upper(msfa)
    y <- upper(direct)
    comparison[[length(comparison) + 1]] <- data.frame(
      Dataset = dataset,
      Zone = zone,
      Spearman = suppressWarnings(cor(x, y, method = "spearman")),
      Sign_agreement_all = mean(sign(x) == sign(y)),
      Sign_agreement_MSFA_abs_ge_0.1 = if (any(abs(x) >= 0.1))
        mean(sign(x[abs(x) >= 0.1]) == sign(y[abs(x) >= 0.1])) else NA_real_,
      MSFA_edges_abs_ge_0.1 = sum(abs(x) >= 0.1),
      Direct_edges_abs_ge_0.1 = sum(abs(y) >= 0.1)
    )
    ij <- which(upper.tri(direct), arr.ind = TRUE)
    edge_output[[length(edge_output) + 1]] <- data.frame(
      Dataset = dataset, Zone = zone,
      Cell_1 = rownames(direct)[ij[, 1]], Cell_2 = colnames(direct)[ij[, 2]],
      Direct_partial_correlation = direct[ij], MSFA_partial_correlation = msfa[ij]
    )
  }
}

comparison <- do.call(rbind, comparison)
edge_output <- do.call(rbind, edge_output)
penalties <- do.call(rbind, penalties)
write.csv(comparison, file.path(out_dir, "msfa_vs_direct_residual_graph.csv"), row.names = FALSE)
write.csv(edge_output, file.path(out_dir, "direct_residual_graph_edges.csv"), row.names = FALSE)
write.csv(penalties, file.path(out_dir, "direct_residual_graph_penalties.csv"), row.names = FALSE)

pdf(file.path(out_dir, "msfa_vs_direct_residual_graph.pdf"), width = 9, height = 5)
par(mfrow = c(1, 2), mar = c(7, 4, 3, 1))
for (dataset in names(datasets)) {
  d <- subset(comparison, Dataset == dataset)
  barplot(d$Spearman, names.arg = d$Zone, las = 2, ylim = c(-1, 1),
          ylab = "Spearman correlation", main = paste(dataset, "MSFA vs direct graph"),
          col = "#4477AA")
  abline(h = 0, lty = 2)
}
dev.off()

print(comparison)
