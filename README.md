# HuNoV-LRPF

R code for the manuscript:

**Quantitative characterization of human norovirus RNA trajectories in zebrafish embryos using a low-complexity B-spline mixed-effects workflow**

HuNoV-LRPF reconstructs batch-clustered `log10` RNA-load trajectories, extracts prespecified trajectory features, and evaluates batch-level prediction and numerical stability.

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

- `01_HuNoV_LRPF_main_analysis.R`: reproduces the manuscript analysis.
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

From R:

```r
source("01_HuNoV_LRPF_main_analysis.R")
```

or from a terminal:

```bash
Rscript 01_HuNoV_LRPF_main_analysis.R
```

The script generates manuscript figures and tables, supplementary outputs, source data, model objects, checksums, metadata, and session information.

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

The locked specification uses a cubic B-spline with df = 4, an internal knot at 3 dpi, boundary knots at 0 and 7 dpi, a batch random intercept, a 0–5-dpi log-scale AUC, prespecified structural QC, conditional fixed-effect simulation, leave-one-batch-out prediction, residual diagnostics, and numerical grid checks.

## Reproducibility checkpoints

A manuscript-matched run should recover approximately:

- 223 development observations from 28 eligible batches;
- M2 df = 4 BIC: 733.74;
- M2 df = 3 BIC: 761.77;
- GZJY122 peak: 8.50 `log10` copies/µL at 3.56 dpi;
- GZJY122 0–5-dpi log-scale AUC: 30.97;
- GZJY122 LOBO RMSE: 0.546;
- GZJY122 LOBO correlation: 0.967;
- GZJY122 95% prediction-interval coverage: 90.0%.

Small platform-dependent numerical differences may occur.

## Interpretation limits

The workflow describes RT-qPCR RNA-load trajectories and is not a mechanistic viral-dynamics or infectivity model. Inferences are restricted to the tested inocula and experimental conditions.

## License

The source code is released under the MIT License. See `LICENSE`.

## Citation

Please cite the associated research article and the Zenodo dataset when using these data. If a fixed software release is archived separately, please also cite the corresponding software DOI.
