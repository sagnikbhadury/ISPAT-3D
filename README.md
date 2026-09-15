# ISPAT3D

ISPAT3D estimates conditional associations among cell-type density variables in registered three-dimensional tumor images. The R package retains the original exact Gaussian-process and CAVI workflow for smaller inputs. The revision workflow in [inst/reanalysis/round2](inst/reanalysis/round2/README.md) uses a Matérn Gaussian process with a 15-neighbor Vecchia likelihood on spatially balanced anchors, then fits the same shared-plus-zone factor covariance form from the full adjusted covariance summaries. It was developed for the larger colorectal and breast imaging applications.

The original [ISPAT](https://github.com/sagnikbhadury/ISPAT) package provides the methodological reference for planar analysis. The revision scripts include a matched section-wise planar implementation with the same selected cells and downstream factor covariance fit, allowing a practical 2D comparison at these data sizes.

Install the R package with:

```r
remotes::install_github("sagnikbhadury/ISPAT-3D")
```

For a small input, the legacy API remains:

```r
library(ISPAT3D)
result <- ISPAT_3D(Y, S, spots_vec, ncores = 8,
                   Kernel = "Matern", MSFA_method = "CAVI")
zone_cov <- result[["spot_1"]]$Zone_Nets[[1]]
pcor <- cov_to_pcor(zone_cov, cell_types = rownames(Y))
```

The revision workflow is script-based because its full-data covariance likelihood is implemented in Python. Its CLI inputs, analysis order, and output provenance are documented in the linked reanalysis directory. Input imaging tables and submitted manuscript files are not included in this software archive. The colorectal and breast source datasets are described in the manuscript and their original releases.
