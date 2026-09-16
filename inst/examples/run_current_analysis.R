# Usage: Rscript run_current_analysis.R crc|bc processed_kde.csv output_dir
suppressPackageStartupMessages({library(ISPAT3D); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || !args[1] %in% c("crc", "bc"))
  stop("Usage: Rscript run_current_analysis.R crc|bc processed_kde.csv output_dir")
dataset <- args[1]; out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
spec <- if (dataset == "crc") list(
  z = "Z_um", section = "section_number", budgets = c(50000L,100000L,150000L),
  kind = "per_zone", labels = c("Tumor","CD8_T","CD4_T","Treg","T_cell",
    "B_cell","Macrophage","Other_immune","Stroma")) else list(
  z = "Z", section = NULL, budgets = c(25000L,45000L,65000L),
  kind = "total", labels = c("Tumor_HER2pos","Tumor_basal","Tumor_luminal",
    "Tumor_other","Endothelial","Macrophage","CAF","Myoepithelial",
    "CD8_T_cell","CD4_T_cell","B_cell","Plasma_cell"))
needed <- unique(c("X","Y",spec$z,spec$section,"pathology_zone",
                   paste0("kde_",spec$labels)))
dat <- data.table::fread(args[2], select = needed)
coords <- as.matrix(dat[, c("X","Y",spec$z), with = FALSE])
sections <- if (is.null(spec$section)) as.integer(floor(dat[[spec$z]] / 2)) else
  dat[[spec$section]]
zones <- dat[["pathology_zone"]]
if (is.numeric(zones)) {
  labels <- c("Very Low","Low","Intermediate","High","Very High")
  zones <- factor(labels[as.integer(zones)], levels = labels)
}
for (tier in seq_along(spec$budgets)) {
  budget <- spec$budgets[tier]
  selection <- ispat3d_sample(coords, zones, sections, budget,
    spec$kind, seed = 2026L + tier)
  rows <- unlist(selection, use.names = FALSE)
  kde <- as.matrix(dat[rows, paste0("kde_",spec$labels), with = FALSE])
  Y <- log1p(1e9 * kde)
  colnames(Y) <- spec$labels
  target <- file.path(out, paste0("tier_",budget))
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  saveRDS(selection, file.path(target,"selected_ids.rds"))
  fit3d <- ispat3d_fit(Y, coords[rows,,drop=FALSE], zones[rows],
    sections[rows], rank = 5L, return_residuals = TRUE)
  saveRDS(fit3d, file.path(target,"fit_3d.rds"))
  if (tier == 2L) {
    fit2d <- ispat3d_fit_2d(Y, coords[rows,,drop=FALSE], zones[rows],
      sections[rows], rank = 5L, return_residuals = TRUE)
    saveRDS(fit2d, file.path(target,"fit_2d.rds"))
  }
  message(dataset," tier ",budget," complete")
}

