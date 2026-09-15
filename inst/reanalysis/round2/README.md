# Scalable ISPAT3D revision workflow

These scripts reproduce the revised analysis stages and figures from processed KDE cell tables. The raw tissue imaging data, processed KDE CSV files, fitted outputs, and manuscript are not part of this software archive.

The colorectal input CSV needs columns `X`, `Y`, `Z_um`, `section_number`, `pathology_zone`, and the nine `kde_<cell_type>` columns named in `run_large_3d.R`. The breast input needs `X`, `Y`, `Z`, `pathology_zone`, and its twelve KDE columns. The breast section identifier is `floor(Z/2)`. Run commands from this directory so R can source `scalable_gp_stage.R`.

The revised main analysis follows the same sequence as the original project: cell-type KDE measurements, spatial adjustment, shared and zone-specific MSFA covariance, partial correlations, the matched serial-section planar comparison, and direct graphical-lasso sensitivity. The computational substitutions are a 15-neighbor Vecchia GP likelihood on spatially balanced anchors and a Gaussian maximum-likelihood fit of the same MSFA covariance form from zone covariance sufficient statistics. These substitutions should be cited as changes in estimation; they are not the original dense GP/CAVI computations.

The 3D budgets are 50,000, 100,000, and 150,000 **per zone** for CRC. Breast budgets are 25,000, 45,000, and 65,000 **total across all five zones**. Each run saves selected source row IDs, adjusted residuals, GP hyperparameters, covariance summaries, and a JSON manifest. The paired 2D comparison uses the 100,000 CRC and 45,000 breast selected row IDs, so its cell roster matches the corresponding middle 3D tier exactly. A section-zone with fewer than ten cells is mean-centered because its GP cannot be estimated reliably; this affected 30 of 45,000 breast cells and none of the CRC cells.

Dependencies: R `data.table`, `gpboost`, `jsonlite`, `glasso`; Python `numpy`, `pandas`, `torch`, `matplotlib`. GPBoost 1.7.4 and R 4.5.2 were used for the revised applications. Python 3.11+ is recommended.

Use paths appropriate to your own copy of the public processed tables:

```text
Rscript run_large_3d.R crc <crc_kde.csv> <analysis_root>
Rscript run_large_3d.R bc <breast_kde.csv> <analysis_root>
python fit_msfa_covariances.py --manifest <analysis_root>/crc/tier_50000/manifest.json --output <analysis_root>/crc/tier_50000/factor_fit --rank 5
python fit_msfa_covariances.py --manifest <analysis_root>/crc/tier_100000/manifest.json --output <analysis_root>/crc/tier_100000/factor_fit --rank 5
python fit_msfa_covariances.py --manifest <analysis_root>/crc/tier_150000/manifest.json --output <analysis_root>/crc/tier_150000/factor_fit --rank 5
python fit_msfa_covariances.py --manifest <analysis_root>/bc/tier_25000/manifest.json --output <analysis_root>/bc/tier_25000/factor_fit --rank 5
python fit_msfa_covariances.py --manifest <analysis_root>/bc/tier_45000/manifest.json --output <analysis_root>/bc/tier_45000/factor_fit --rank 5
python fit_msfa_covariances.py --manifest <analysis_root>/bc/tier_65000/manifest.json --output <analysis_root>/bc/tier_65000/factor_fit --rank 5
python combine_and_plot.py --root <analysis_root>/crc --budgets 50000 100000 150000 --dataset crc
python combine_and_plot.py --root <analysis_root>/bc --budgets 25000 45000 65000 --dataset bc
Rscript run_sectionwise_2d.R crc 100000 <crc_kde.csv> <analysis_root>/crc/tier_100000/selected_ids.rds <analysis_root>/crc/sectionwise_2d_min10_tier_100000
Rscript run_sectionwise_2d.R bc 45000 <breast_kde.csv> <analysis_root>/bc/tier_45000/selected_ids.rds <analysis_root>/bc/sectionwise_2d_min10_tier_45000
python fit_msfa_covariances.py --manifest <analysis_root>/crc/sectionwise_2d_min10_tier_100000/manifest.json --output <analysis_root>/crc/sectionwise_2d_min10_tier_100000/factor_fit --rank 5
python fit_msfa_covariances.py --manifest <analysis_root>/bc/sectionwise_2d_min10_tier_45000/manifest.json --output <analysis_root>/bc/sectionwise_2d_min10_tier_45000/factor_fit --rank 5
python plot_2d_deltas.py --root <analysis_root> --dataset crc --budget 100000
python plot_2d_deltas.py --root <analysis_root> --dataset bc --budget 45000
Rscript precision_sensitivity.R <analysis_root>
```

The `run_saved_simulations.R`, `summarize_simulations.py`, `measure_simulation_runtime.py`, `focused_gp_vecchia_validation.R`, and `plot_focused_gp_validation.py` scripts remake Figures 3–6 from the saved generated simulation inputs or the focused synthetic design. `make_overview_figures.py` remakes Figures 1–2. `make_section_stats.R`, `block_bootstrap_precision.py`, and `plot_fullpair_with_block.py` produce the section-block stability tables and full-pair zone heatmaps. The section-block resampling uses fixed GP residuals; it does not re-estimate the GP or account for image-registration uncertainty.

The precision comparison applies the first-round EBIC graphical-lasso protocol to the new large-fit residuals. A low-rank shared loading product has no unregularized precision matrix, so biological interpretation and 2D deltas use full zone covariances.
