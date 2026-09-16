# ISPAT3D

ISPAT3D estimates how **conditional associations among local cell-type densities change across ordered tumor-burden zones** in a registered three-dimensional multiplex image. It accepts cell-level neighborhood-density measurements and registered coordinates; it does not perform image segmentation, cell classification, section registration, or KDE estimation.

The current package implements the analysis used in the manuscript:

1. Fit a separate anisotropic Matérn-3/2 GP to each density variable within each tumor-burden zone. A 15-neighbor Vecchia likelihood is fitted on spatially balanced anchors, then the fitted field is predicted at **every selected cell**.
2. Compute the complete residual covariance in each zone.
3. Fit `Σ_q = ΦΦᵀ + Λ_qΛ_qᵀ + Ψ_q` by Gaussian covariance maximum likelihood, using the actual zone row counts and full residual covariance summaries.
4. Invert each **full** fitted zone covariance to report signed partial correlations. The optional matched 2D analysis fits planar GPs independently by section on the **same selected cells**, then uses the same covariance fitter.

The package API is R-only. The historical exact-GP, CAVI/SVI, and bigmemory APIs have been removed from the current package.

## Requirements and installation

Use R 4.5 or a compatible recent version. Install GPBoost from CRAN, then the GitHub package:

```r
install.packages("gpboost")
install.packages("remotes")
remotes::install_github("sagnikbhadury/ISPAT-3D")
library(ISPAT3D)
packageVersion("ISPAT3D")
```

Version 0.3.0 introduces the current package API. `gpboost` is imported automatically when ISPAT3D is installed. The manuscript applications used R 4.5.2 and GPBoost 1.7.4. For the large CSV example below, install `data.table` separately:

```r
install.packages("data.table")
```

You can also download the repository source and run `R CMD build .` followed by `R CMD check --no-manual ISPAT3D_0.3.0.tar.gz`. An R source install does not require Python.

## Prepare the inputs

Every selected cell is one **row** in all inputs, in exactly the same order:

| Input | Shape and meaning |
| --- | --- |
| `Y` | Numeric **N × G** matrix: cells by cell-type density variables. Give columns unique cell-type names. |
| `coords` | Numeric **N × 3** matrix: registered `x`, `y`, `z` coordinates. Use a common spatial unit across axes. |
| `zones` | Length-N factor, character, or numeric vector of tumor-burden zone labels. The usual order is Very Low, Low, Intermediate, High, Very High. |
| `sections` | Length-N section identifier. Needed for balanced 3D anchor selection and required for the matched 2D fit. |

For the manuscript, `Y` contains `log1p(1e9 * KDE)` for each cell type, evaluated at every selected cell. KDE surfaces were estimated within each section. The tumor-burden variable used to create `zones` is separate from the response matrix: CRC zones come from quintiles of the tumor-cell KDE; breast zones come from quintiles of pan-cytokeratin intensity. A cell may have its own tumor-cell KDE among the `Y` variables and still be assigned a zone using that tumor-burden score. Zone labels are **relative within a specimen**, not universal biological states.

The processed CRC CSV expected by the included example has `X`, `Y`, `Z_um`, `section_number`, `pathology_zone`, and `kde_` columns for Tumor, CD8_T, CD4_T, Treg, T_cell, B_cell, Macrophage, Other_immune, and Stroma. The processed breast CSV has `X`, `Y`, `Z`, `pathology_zone`, and `kde_` columns for Tumor_HER2pos, Tumor_basal, Tumor_luminal, Tumor_other, Endothelial, Macrophage, CAF, Myoepithelial, CD8_T_cell, CD4_T_cell, B_cell, and Plasma_cell. In that specimen, section IDs are computed as `floor(Z/2)`. Raw images and processed tables are not included in this software repository; use the data releases cited in the manuscript.

### Sample the cell roster

`ispat3d_sample()` selects cells proportionally across sections within each zone and spreads picks over a 10 × 10 in-plane grid. It returns **original row indices**, so save them for exact 3D/2D matching.

```r
library(ISPAT3D)

# coords, zones, and sections refer to every eligible source cell.
selected <- ispat3d_sample(
  coords, zones, sections,
  budget = 100000L,
  budget_kind = "per_zone",
  seed = 2027L
)
saveRDS(selected, "selected_ids.rds")
rows <- unlist(selected, use.names = FALSE)
```

Use `budget_kind = "per_zone"` for a CRC-style budget, or `budget_kind = "total"` to allocate one breast-style budget proportionally over zones. The manuscript used **50,000, 100,000, and 150,000 cells per CRC zone** (250,000–750,000 cells total), and **25,000, 45,000, and 65,000 breast cells total across all five zones**. The available breast specimen has 73,936 cells overall. A budget greater than a zone's population is rejected.

### Fit the 3D model

```r
# kde has the same source rows as coords and one column per cell type.
Y <- log1p(1e9 * as.matrix(kde[rows, , drop = FALSE]))
fit3d <- ispat3d_fit(
  Y = Y,
  coords = coords[rows, , drop = FALSE],
  zones = zones[rows],
  sections = sections[rows],
  rank = 5L,
  return_residuals = TRUE
)
saveRDS(fit3d, "fit_3d.rds")
```

`ispat3d_fit()` fits a Matérn-3/2 GP separately for every zone-variable pair using up to 5,000 spatially balanced anchors, approximately 10% of each zone's selected cells and at least 300 when available. It uses 15 Vecchia neighbors and at most 30 GP optimizer iterations by default. Predictions and residual covariance use **all selected cells**, not only anchors. Constant variables are mean-centered. The shared and zone-specific covariance loading ranks are both `rank`, with default 5 for the real-data analyses.

The main return fields are:

| Field | Meaning |
| --- | --- |
| `shared` | Shared factor covariance `ΦΦᵀ`; it is low rank and has no separately interpreted partial-correlation graph. |
| `full` | Named list of full fitted zone covariances `Σ_q`. |
| `partial` | Named list of signed partial-correlation matrices from the full zone covariances. |
| `covariances` | Complete GP-adjusted sample covariance for each zone. |
| `gp_log` | Per-variable GP status, anchor count, fit time, and prediction time. |
| `counts` | Actual selected row count per zone. |
| `active_variables` / `inactive_variables` | Variables included in or excluded from the covariance likelihood. |
| `residuals` | Cell-level GP-adjusted matrices, only when `return_residuals = TRUE`. |

A constant variable is excluded from the active covariance likelihood and re-embedded as a zero covariance row and column for label provenance. Its partial-correlation diagonal is a placeholder 1; do not interpret pairs involving it. This applies to the structurally zero breast luminal-tumor KDE.

```r
fit3d$partial[["Very High"]]
fit3d$gp_log
fit3d$active_variables
```

Positive and negative entries are **residual conditional density associations given the other measured variables**. They are not direct observations of cell contact, ligand–receptor signaling, immune suppression, or causality. No edge-wise p-values or confidence intervals are produced.

### Fit the matched section-wise 2D comparator

Use the **same** `Y`, `coords`, `zones`, `sections`, and row order. The planar function takes the x–y columns of `coords` and fits separate GPs within each zone-section group.

```r
fit2d <- ispat3d_fit_2d(
  Y, coords[rows, , drop = FALSE], zones[rows], sections[rows],
  rank = 5L, return_residuals = TRUE
)
saveRDS(fit2d, "fit_2d.rds")

zone <- "Very High"
delta <- fit3d$partial[[zone]] - fit2d$partial[[zone]]
```

Groups with fewer than 10 cells or constant measurements are mean-centered. A failed planar GP is also mean-centered and recorded with a `fallback:` status in `gp_log`; inspect that log before interpreting contrasts. This comparison holds the cell roster and downstream covariance estimator fixed, but changes **both coordinate dimension and section pooling**. It is not a Z-only causal ablation.

### Fit from precomputed residual covariance

If GP residual covariances have already been saved, the factor stage is independently available:

```r
fit_cov <- ispat3d_fit_covariance(
  covariances = list(
    "Very Low" = cov(residuals_very_low),
    "Low" = cov(residuals_low),
    "Intermediate" = cov(residuals_intermediate),
    "High" = cov(residuals_high),
    "Very High" = cov(residuals_very_high)
  ),
  counts = c(nrow(residuals_very_low), nrow(residuals_low),
             nrow(residuals_intermediate), nrow(residuals_high),
             nrow(residuals_very_high)),
  rank = 5L
)
fit_cov$partial[["High"]]
```

This objective uses the complete sample covariance in every zone, weights zones by their actual row counts, standardizes variables by pooled residual standard deviation, and optimizes the loading matrices and positive uniquenesses by L-BFGS-B with analytic gradients. It is a Gaussian covariance likelihood, **not** the historical CAVI/SVI posterior. `ispat3d_partial_correlation()` converts any full positive-definite covariance to a partial-correlation matrix using a small numerical ridge.

## Run the included CRC or breast example

After preparing a processed KDE CSV with the columns listed above:

```text
Rscript inst/examples/run_current_analysis.R crc path/to/crc_kde.csv path/to/crc_output
Rscript inst/examples/run_current_analysis.R bc path/to/breast_kde.csv path/to/breast_output
```

The example runs all three dataset-specific sampling budgets, saves `selected_ids.rds` and `fit_3d.rds` in each tier directory, and saves `fit_2d.rds` for the middle tier. It uses the packaged R fitter; no Python covariance step is needed. These full application runs can take substantial time and memory, especially the CRC tiers. Start with a smaller budget to check your input columns and GPBoost installation.

## Reproduce the submitted figures

The `inst/reanalysis/round2` directory retains the analysis and figure scripts used for the submitted results, including saved simulation reruns, the focused GP experiment, section-block sign stability, matched 2D fitting, and graphical-lasso sensitivity. Its `README.md` lists the CLI inputs and analysis order. Those archived runs used a Python implementation of the same Gaussian covariance objective; the R package uses an analytic-gradient R implementation. Optimization and floating-point differences mean a fresh package fit may not reproduce every saved coefficient bit-for-bit. The original registered image data, processed KDE tables, fitted outputs, and manuscript are not redistributed in this code repository.

## Interpretation and troubleshooting

- **Matrix orientation:** `Y` must be cells × variables; `coords` must be cells × 3. If your input is variables × cells, transpose it before calling `ispat3d_fit()`.
- **Zone counts:** Each zone must have more selected rows than active variables. Check `table(zones[rows])` before fitting.
- **Coordinate scale:** x, y, and z must use a coherent physical scale. Anisotropic ranges are fitted separately, but incorrect units or registration can still distort GP adjustment.
- **Planar fallback:** Inspect `table(fit2d$gp_log$status)`. A `fallback:` status means that particular zone-section-variable fit was mean-centered.
- **Memory:** Use `return_residuals = FALSE` unless cell-level residuals are needed. The covariance fit itself works from small G × G summaries, but GP prediction still processes every selected cell.
- **Interpretation:** A smoother GP can remove real broad-scale biology as well as nuisance variation. Compare scales or use independent biological evidence before treating a conditional edge as a mechanism.

For package function signatures, use `?ispat3d_fit`, `?ispat3d_fit_2d`, `?ispat3d_sample`, and `?ispat3d_fit_covariance` in R.

