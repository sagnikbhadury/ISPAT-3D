# ISPAT3D: Spatial Network Estimation in 3D Multiplexed Cancer Imaging

An R package extending [ISPAT](https://github.com/sagnikbhadury/ISPAT) to
volumetric (3D) or serial-section multiplexed tissue imaging. Gaussian
process hyperparameters are estimated per marker per tumor-density zone via
type-II marginal likelihood optimization (L-BFGS-B); the resulting
de-spatialised residuals are passed to a multi-study factor analysis to
recover a shared cell-type interaction network and zone-specific interaction
networks for each spot/field of view.

## Installation

```r
# install.packages("remotes")
remotes::install_github("sagnikbhadury/ISPAT-3D")
```

`run_ispat_spot_bigmem()` additionally requires the `bigmemory` package, and
`plot_pcor_chord()` requires `dplyr`, `igraph`, `tidygraph`, `ggraph`,
`ggplot2`, and `scales`; none of these are installed automatically, since
they are only needed for those specific optional functions.

## Usage

```r
library(ISPAT3D)

result <- ISPAT_3D(Y, S, spots_vec, ncores = 8, Kernel = "Matern", MSFA_method = "CAVI")

# Shared network for one spot, as a partial correlation matrix
pcor <- cov_to_pcor(result[["spot_1"]]$Shared_Net, cell_types = rownames(Y))
plot_pcor_chord(pcor, title = "Shared network", threshold = 0.05)
```

For very large spots (many cells and/or many markers), `run_ispat_spot_bigmem()`
is a memory-efficient drop-in replacement for `run_ispat_spot()` that avoids
serializing the full expression matrix to every parallel worker.

## Package layout

* `R/ISPAT3D.R` — `ISPAT_3D()` (loops over spots) and `run_ispat_spot()`
  (single spot).
* `R/bigmem.R` — `run_ispat_spot_bigmem()`, a memory-mapped variant for large
  spots.
* `R/kernel-utils.R` — internal GP kernel and per-marker fitting helpers.
* `R/msfa.R` — the multi-study factor analysis step, `cavi_msfa()` and
  `svi_msfa()` (shared with the ISPAT package).
* `R/network-utils.R` — `cov_to_pcor()` and `plot_pcor_chord()`, utilities
  for turning a network into a partial correlation matrix and chord diagram.
* `inst/examples/` — two case-study analysis scripts (breast cancer IMC data
  from Kuett et al. 2022, and CRC imaging data from Lin et al. 2023) that use
  this package end to end. These are example pipelines, not part of the
  package's public API, and expect their own local input data.

See `?ISPAT_3D` for full documentation of the main function.
