# =============================================================================
# CRC1 ISPat-3D — Sequential Cluster Script
#
# Example analysis script using the ISPAT3D package; not part of the
# installable package itself (see inst/examples/).
#
# SLURM header (add to top when submitting):
#   #SBATCH --job-name=ispat3d_crc1
#   #SBATCH --ntasks=1
#   #SBATCH --cpus-per-task=1
#   #SBATCH --mem=16G
#   #SBATCH --time=02:00:00
#   #SBATCH --output=ispat3d_crc1_%j.log
#
# Run with:
#   Rscript CRC_ISPat3D_cluster.R
#
# Output:
#   ispat3d_crc1_result.rds   — full ISPat-3D result
#   pcor_<zone>.csv           — partial correlation matrices
#   network_<zone>.png/.pdf   — chord diagrams
# =============================================================================

pacman::p_load(data.table, tidyverse, ggplot2, ggraph, igraph,
               tidygraph, Matrix, matrixcalc)

library(ISPAT3D)

# =============================================================================
# PATHS AND CONSTANTS
# =============================================================================

DATA_PATH <- "./lin2023_crc_output/lin2023_crc1_kde_surfaces.csv"
OUT_DIR   <- "./lin2023_crc_output/ispat3d_results/"
FIG_DIR   <- "./lin2023_crc_output/figures/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

CELL_TYPES  <- c("Tumor", "CD8_T", "CD4_T", "Treg", "T_cell",
                  "B_cell", "Macrophage", "Other_immune", "Stroma")
G           <- length(CELL_TYPES)
ZONE_LABELS <- c("Very Low", "Low", "Intermediate", "High", "Very High")
MAX_PER_ZONE <- 150000
set.seed(2026)

# =============================================================================
# STEP 1 — LOAD, SUBSAMPLE, GP REGRESSION (sequential for loop)
#
# For each zone c, for each marker g:
#   - Fit GP via L-BFGS-B marginal likelihood (inherited from ISPAT_3D.R)
#   - Compute spatially adjusted residual:
#       z_hat = (y - beta_hat) - B_hat
#             = sigmaEPS^2 * V^{-1} * (y - beta_hat)
#   - Store z_hat as column g of Z_est[[c]]  (N_c x G)
# =============================================================================

message("\n=== STEP 1: Load data and run GP regression ===\n")

message("Loading KDE surface data ...")
dat <- data.table::fread(DATA_PATH, data.table = FALSE)
message(sprintf("Loaded: %d cells x %d columns", nrow(dat), ncol(dat)))

# --- Build Y and S ---
kde_cols <- paste0("kde_", CELL_TYPES)
Y <- t(as.matrix(dat[, kde_cols]))
Y <- log1p(Y * 1e9)          # stabilise extreme KDE dynamic range
rownames(Y) <- CELL_TYPES

S <- as.matrix(dat[, c("X", "Y", "Z_um", "pathology_zone")])
colnames(S) <- c("X", "Y", "Z", "zone")

rm(dat); gc()                 # free the 4 GB data.frame immediately
message("dat freed from memory.")

# --- Spatially stratified subsample within each zone ---
message("\nSubsampling zones ...")
zones_present <- sort(unique(S[, 4]))
keep_idx      <- c()

for (z in zones_present) {
  z_idx <- which(S[, 4] == z)
  n_z   <- length(z_idx)
  if (n_z <= MAX_PER_ZONE) {
    keep_idx <- c(keep_idx, z_idx)
    message(sprintf("  Zone %d: %d cells — keeping all", z, n_z))
  } else {
    s_z    <- S[z_idx, 1:2]
    x_cut  <- cut(s_z[, 1], breaks = ceiling(sqrt(MAX_PER_ZONE)))
    y_cut  <- cut(s_z[, 2], breaks = ceiling(sqrt(MAX_PER_ZONE)))
    strata <- paste(x_cut, y_cut, sep = "_")
    s_idx  <- tapply(seq_along(z_idx), strata, function(ii) {
      sample(ii, min(max(1L, round(MAX_PER_ZONE * length(ii) / n_z)), length(ii)))
    })
    s_idx  <- z_idx[unlist(s_idx)]
    if (length(s_idx) > MAX_PER_ZONE) s_idx <- sample(s_idx, MAX_PER_ZONE)
    keep_idx <- c(keep_idx, s_idx)
    message(sprintf("  Zone %d: %d cells — subsampled to %d", z, n_z, length(s_idx)))
  }
}

Y_sub <- Y[, keep_idx, drop = FALSE]
S_sub <- S[keep_idx, , drop = FALSE]
rm(Y, S); gc()
message(sprintf("\nY_sub: %d x %d | S_sub: %d x 4", nrow(Y_sub), ncol(Y_sub), nrow(S_sub)))

# --- Sequential GP regression ---
# Z_est[[c]] is N_c x G — spatially adjusted residuals for zone c
Z_est <- vector("list", length(zones_present))
names(Z_est) <- ZONE_LABELS[zones_present]

for (c in seq_along(zones_present)) {

  z       <- zones_present[c]
  pos_c   <- which(S_sub[, 4] == z)
  N_c     <- length(pos_c)
  S_loc_c <- S_sub[pos_c, 1:3, drop = FALSE]

  message(sprintf("\n  Zone %d (%s): %d cells",
                  z, ZONE_LABELS[z], N_c))

  Z_zone <- matrix(0, nrow = N_c, ncol = G,
                   dimnames = list(NULL, CELL_TYPES))

  # Lengthscale initialisation — computed once per zone, shared across markers
  xy_d <- as.vector(dist(S_loc_c[, 1:2])); xy_d <- xy_d[xy_d > 0]
  z_d  <- as.vector(dist(S_loc_c[, 3, drop = FALSE])); z_d <- z_d[z_d > 0]
  lS0  <- if (length(xy_d) > 0) median(xy_d) else 0.1
  lZ0  <- if (length(z_d)  > 0) median(z_d)  else lS0

  for (g in 1:G) {

    y <- as.vector(Y_sub[g, pos_c])

    # L-BFGS-B marginal likelihood optimisation
    opt <- tryCatch(
      optim(
        par     = c(log(lS0), log(lZ0), log(1), log(0.5), 0),
        fn      = .gp_neglml,
        y       = y, S_loc_c = S_loc_c, Kernel = "Matern",
        method  = "L-BFGS-B",
        lower   = c(log(lS0/10), log(lZ0/10), log(0.01), log(0.01), -10),
        upper   = c(log(lS0*10), log(lZ0*10), log(10),   log(10),    10),
        control = list(maxit = 200, factr = 1e9)
      ),
      error = function(e) {
        message(sprintf("    [Zone %d, %s] optim failed: %s — using heuristic",
                        z, CELL_TYPES[g], conditionMessage(e)))
        list(par = c(log(lS0), log(lZ0), log(1), log(0.5), mean(y)))
      }
    )

    lS_est       <- exp(opt$par[1])
    lZ_est       <- exp(opt$par[2])
    sigmaS_est   <- exp(opt$par[3])
    sigmaEPS_est <- exp(opt$par[4])
    beta_est     <- opt$par[5]

    KS_est <- .build_kernel(S_loc_c, lS_est, lZ_est, "Matern")
    if (!matrixcalc::is.positive.definite(KS_est))
      KS_est <- as.matrix(Matrix::nearPD(KS_est)$mat)

    sS  <- max(sigmaS_est,   0.011)
    sEP <- max(sigmaEPS_est, 0.011)

    V     <- sS^2 * KS_est + sEP^2 * diag(N_c)
    cholV <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(cholV))
      cholV <- tryCatch(chol(V + 1e-6 * diag(N_c)),
                        error = function(e) chol(sEP^2 * diag(N_c)))

    alpha  <- backsolve(cholV, forwardsolve(t(cholV), y - beta_est))
    B_hat  <- sS^2 * KS_est %*% alpha

    # Spatially adjusted residual: z_hat = (y - beta) - B_hat
    Z_zone[, g] <- as.vector((y - beta_est) - B_hat)

    # Free large matrices immediately to keep peak memory low
    rm(KS_est, V, cholV, alpha, B_hat); gc()

    message(sprintf("    [%s] lS=%.1f lZ=%.1f sigmaS=%.3f sigmaEPS=%.3f",
                    CELL_TYPES[g], lS_est, lZ_est, sigmaS_est, sigmaEPS_est))
  }

  Z_est[[c]] <- Z_zone
  message(sprintf("  Zone %d GP done. Z_est dim: %d x %d",
                  z, nrow(Z_zone), ncol(Z_zone)))
}

# Checkpoint: save GP residuals before MSFA in case of crash
saveRDS(Z_est, file.path(OUT_DIR, "Z_est_gp_residuals.rds"))
message("\nGP residuals saved to Z_est_gp_residuals.rds")

# =============================================================================
# STEP 2 — COMBINE GP RESIDUALS FOR MSFA
#
# Z_est is already in the correct format for cavi_msfa / svi_msfa:
#   list of C matrices, each N_c x G
# No further transformation needed — just verify dimensions.
# =============================================================================

message("\n=== STEP 2: Verify GP residual structure ===\n")

for (c in seq_along(Z_est)) {
  message(sprintf("  Z_est[[%d]] (%s): %d x %d  |  col means range [%.3f, %.3f]",
                  c, names(Z_est)[c],
                  nrow(Z_est[[c]]), ncol(Z_est[[c]]),
                  min(colMeans(Z_est[[c]])),
                  max(colMeans(Z_est[[c]]))))
}

# =============================================================================
# STEP 3 — MSFA AND NETWORK EXTRACTION
#
# MSFA decomposes Z_est into:
#   Shared factors Phi (G x K_shared) — common across all zones
#   Zone factors Lambda_c (G x K_c)   — zone-specific
#
# Networks:
#   Shared_Net  = Phi Phi^T                        (G x G)
#   Zone_Net_c  = Phi Phi^T + Lambda_c Lambda_c^T + Psi_c  (G x G)
# =============================================================================

message("\n=== STEP 3: MSFA and network extraction ===\n")

K_shared <- ceiling(2 * log(G))
K_zone   <- rep(ceiling(2 * log(G)), length(Z_est))
message(sprintf("K_shared = %d | K_zone = %d (per zone)", K_shared, K_zone[1]))

VBfit <- cavi_msfa(Z_est, K_shared, K_zone, scale = FALSE)

Shared_Net <- tcrossprod(VBfit$mean_phi)
rownames(Shared_Net) <- colnames(Shared_Net) <- CELL_TYPES

Zone_Nets <- lapply(seq_along(Z_est), function(c) {
  net <- tcrossprod(VBfit$mean_phi) +
    tcrossprod(VBfit$mean_lambda_s[[c]]) +
    diag(as.numeric(VBfit$mean_psi_s[[c]]))
  rownames(net) <- colnames(net) <- CELL_TYPES
  net
})
names(Zone_Nets) <- ZONE_LABELS[zones_present]

ispat_result <- list(
  Shared_Net = Shared_Net,
  Zone_Nets  = Zone_Nets,
  Z_est      = Z_est,
  VBfit      = VBfit
)

saveRDS(ispat_result, file.path(OUT_DIR, "ispat3d_crc1_result.rds"))
message("ISPat-3D result saved.")

# =============================================================================
# STEP 4 — PARTIAL CORRELATIONS AND NETWORK FIGURES
#
# Covariance -> partial correlation:
#   Theta    = Sigma^{-1}
#   pcor_ij  = -Theta_ij / sqrt(Theta_ii * Theta_jj)
#
# Chord diagram: circular layout, edges coloured by sign of pcor.
#   Positive (blue): conditionally co-localised cell types
#   Negative (red):  conditional spatial exclusion
# =============================================================================

message("\n=== STEP 4: Partial correlations and figures ===\n")

## cov_to_pcor() is provided by the ISPAT3D package.

all_networks <- c(list(Shared = Shared_Net), Zone_Nets)
pcor_list    <- lapply(all_networks, cov_to_pcor)

# Save partial correlation CSVs
for (nm in names(pcor_list)) {
  write.csv(round(pcor_list[[nm]], 6),
            file.path(OUT_DIR, paste0("pcor_", gsub(" ", "_", nm), ".csv")))
}
message("Partial correlation matrices saved.")

# --- Chord diagram function ---
CELL_TYPE_COLORS <- c(
  Tumor         = "#E63946",
  CD8_T         = "#2196F3",
  CD4_T         = "#4CAF50",
  Treg          = "#FF9800",
  T_cell        = "#9C27B0",
  B_cell        = "#00BCD4",
  Macrophage    = "#795548",
  Other_immune  = "#607D8B",
  Stroma        = "#CDDC39"
)

## plot_pcor_chord() is provided by the ISPAT3D package (pass node_colors =
## CELL_TYPE_COLORS below to keep this study's custom per-cell-type palette).

# Adaptive threshold: 30th percentile of all |pcor| values
all_abs_pcor <- unlist(lapply(pcor_list, function(m) abs(m[upper.tri(m)])))
threshold    <- max(0.03, quantile(all_abs_pcor, 0.30, na.rm = TRUE))
message(sprintf("Edge threshold: |pcor| >= %.4f", threshold))

network_titles <- c(
  Shared        = "Shared Network (All Zones)",
  "Very Low"    = "Zone: Very Low Tumor Density",
  Low           = "Zone: Low Tumor Density",
  Intermediate  = "Zone: Intermediate Tumor Density",
  High          = "Zone: High Tumor Density",
  "Very High"   = "Zone: Very High Tumor Density"
)

for (nm in names(pcor_list)) {
  title_str <- network_titles[nm]
  message(sprintf("  Plotting: %s", nm))

  p <- plot_pcor_chord(pcor_list[[nm]], title = title_str, threshold = threshold,
                       node_colors = CELL_TYPE_COLORS)
  if (is.null(p)) next

  fname <- file.path(FIG_DIR, paste0("network_", gsub("[^A-Za-z0-9]", "_", nm)))
  ggsave(paste0(fname, ".png"), plot = p, width = 9, height = 9,
         dpi = 300, bg = "white")
  ggsave(paste0(fname, ".pdf"), plot = p, width = 9, height = 9,
         device = cairo_pdf)
  message(sprintf("    Saved: %s.png / .pdf", basename(fname)))
}

# Summary of top interactions per network
message("\nTop 5 partial correlations per network:\n")
for (nm in names(pcor_list)) {
  m <- pcor_list[[nm]]
  top <- expand.grid(from = CELL_TYPES, to = CELL_TYPES,
                     stringsAsFactors = FALSE) %>%
    dplyr::filter(from < to) %>%
    dplyr::mutate(pcor = mapply(function(f, t) m[f, t], from, to)) %>%
    dplyr::arrange(dplyr::desc(abs(pcor))) %>%
    dplyr::slice_head(n = 5)
  cat(sprintf("--- %s ---\n", nm)); print(top); cat("\n")
}

message("\n=== Analysis complete ===")
message(sprintf("Results: %s", OUT_DIR))
message(sprintf("Figures: %s", FIG_DIR))
