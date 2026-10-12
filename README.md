# HuNoV-LRPF

R code for the manuscript:

**A constrained B-spline mixed-effects workflow for batch-structured analysis of human norovirus RNA trajectories in zebrafish embryos**

HuNoV-LRPF is a reusable analysis workflow based on cubic B-splines and mixed-effects modelling. It fits batch-structured `log10` RNA-load trajectories, extracts prespecified trajectory features, and evaluates batch-level prediction and numerical stability.

## Repository contents

```text
HuNoV-LRPF/
├── README.md
├── LICENSE
├── 01_HuNoV_LRPF_main_analysis.R
├── 02_HuNoV_LRPF_single_locked_application.R
└── R/
    ├── HuNoV_LRPF_core.R
    ├── HuNoV_LRPF_plots.R
    └── HuNoV_LRPF_reporting.R
```

- `01_HuNoV_LRPF_main_analysis.R`: reproduces the manuscript analysis, including the GZJY122 sensitivity analyses.
- `02_HuNoV_LRPF_single_locked_application.R`: applies the locked manuscript specification to one additional inoculum without repeating model or spline selection.
- `R/`: shared analysis, plotting, and reporting functions.

## Requirements

The manuscript analysis was run with R 4.4.1 and `lme4` 1.1-37. Required packages are:

```r
install.packages(c(
  "dplyr", "tidyr", "purrr", "tibble", "readr", "ggplot2",
  "lme4", "MASS", "DescTools", "Matrix", "performance",
  "forcats", "patchwork", "scales"
))
```

The base R package `splines` is also used.

## Data

The analytical data are deposited separately in Zenodo:

**Dataset DOI:** [10.5281/zenodo.21502817](https://doi.org/10.5281/zenodo.21502817)

Required files:

```text
df_315.csv
df_324.csv
df_17.csv
df_19.csv
df_29.csv
df_122.csv
```

Place the CSV files beside the entry script, in a `data/` or `input/` folder, or specify their directory in the script configuration.

Each CSV must contain:

| Column | Description |
|---|---|
| `batch` | Independent experimental-run identifier |
| `dpi` | Days post injection |
| `VL` | Positive RT-qPCR RNA concentration in copies/µL before `log10` transformation |

## Run the manuscript analysis

M2 was prespecified to allow peak timing and local slopes to vary among inocula. BIC selected df = 4 within M2, while M1 had the lowest BIC across the model hierarchy.

From R:

```r
source("01_HuNoV_LRPF_main_analysis.R")
```

or from a terminal:

```bash
Rscript 01_HuNoV_LRPF_main_analysis.R
```

The script generates manuscript figures and tables, supplementary outputs (including Table S17), source data, model objects, checksums, metadata, and session information.

For GZJY122, the script also runs two sensitivity analyses using the locked spline basis: a fixed-effect-only refit and feature extraction from the five four-batch LOBO training fits. These analyses assess the numerical stability of the population-level feature point estimates without repeating model or spline selection.

## Apply the locked specification to one new inoculum

Edit the following values in `02_HuNoV_LRPF_single_locked_application.R`:

```r
inoculum_id <- "NEW_INOCULUM"
genotype_ptype <- "GENOTYPE[P-TYPE]"
input_file <- "new_inoculum.csv"
```

Then run:

```bash
Rscript 02_HuNoV_LRPF_single_locked_application.R
```

The locked specification uses a cubic B-spline with df = 4, an internal knot at 3 dpi, boundary knots at 0 and 7 dpi, a batch random intercept, a 0–5-dpi log-scale AUC, prespecified structural QC, conditional fixed-effect simulation, leave-one-batch-out prediction, residual diagnostics, and numerical grid checks. Model coefficients and variance parameters are estimated from the new data.

## Reproducibility checkpoints

A manuscript-matched run should recover approximately:

- 223 development observations from 28 eligible batches;
- M2 df = 4 BIC: 733.74;
- M2 df = 3 BIC: 761.77;
- M1 df = 4 BIC: 702.91;
- GZJY122 fitted peak: 8.50 `log10` copies/µL at 3.56 dpi;
- GZJY122 0–5-dpi log-scale AUC: 30.97 `log10(copies/µL)·day`;
- GZJY122 LOBO RMSE: 0.546;
- GZJY122 LOBO correlation: 0.967;
- GZJY122 LOBO coverage of nominal 95% prediction intervals: 90.0% (36/40);
- GZJY122 batch random-intercept variance: 0.0000 in the full-data fit, with all five LOBO training fits at the random-effect boundary;
- identical primary mixed-model and fixed-effect-only feature estimates to the reported precision, with all five four-batch sensitivity refits completed successfully.

Small platform-dependent numerical differences may occur.

## Interpretation

The workflow describes RT-qPCR RNA-load trajectories and is not a mechanistic viral-dynamics or infectivity model. Given the daily sampling schedule, fitted peak times are interpreted as a 3–4-dpi window; decimal values are retained above for numerical reproducibility. The GZJY122 sensitivity analyses assess the stability of population-level feature estimates; the five-batch dataset did not resolve a separate batch variance component. Inferences are restricted to the tested inocula and experimental conditions.

## License

The source code is released under the MIT License. See `LICENSE`.

## Citation

Please cite the associated research article and the Zenodo dataset when using these data. If a fixed software release is archived separately, please also cite the corresponding software DOI.
