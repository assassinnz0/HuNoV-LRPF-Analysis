# HuNoV-LRPF standalone locked-specification analysis (version 10.3)
#
# Purpose:
#   Apply the manuscript specification to one new inoculum without repeating
#   development-stage model or spline selection.
#
# Required input CSV columns:
#   batch, dpi, VL
#
# Locked specification:
#   cubic B-spline, df = 4, internal knot = 3 dpi, boundary knots = 0 and 7 dpi,
#   batch random intercept, 0-5 dpi AUC, structural QC based only on sampling
#   coverage and data integrity, conditional fixed-effect simulation, LOBO
#   prediction, residual diagnostics, and numerical grid checks.

locate_entry_script <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) return(sub("^--file=", "", file_arg[1]))

  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      value <- frame$ofile
      if (is.null(value) || length(value) == 0L) NA_character_ else as.character(value[1])
    },
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files) & nzchar(frame_files)]
  if (length(frame_files) == 0L) return(NA_character_)
  tail(frame_files, 1)
}

entry_script <- locate_entry_script()
project_dir <- if (is.na(entry_script)) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  dirname(normalizePath(entry_script, winslash = "/", mustWork = TRUE))
}

source(file.path(project_dir, "R", "HuNoV_LRPF_core.R"))
source(file.path(project_dir, "R", "HuNoV_LRPF_plots.R"))
source(file.path(project_dir, "R", "HuNoV_LRPF_reporting.R"))
require_hunov_packages()
require_hunov_plot_packages()

# =============================================================================
# 1. User configuration
# =============================================================================

# Edit these values before running the script.
inoculum_id <- "122"
genotype_ptype <- "GII.4"
input_file <- "df_122.csv"

if (
  identical(inoculum_id, "NEW_INOCULUM") ||
  identical(genotype_ptype, "GENOTYPE[P-TYPE]") ||
  identical(input_file, "new_inoculum.csv")
) {
  stop(
    "Edit inoculum_id, genotype_ptype, and input_file before running the script.",
    call. = FALSE
  )
}

# Set data_dir to an explicit folder if the CSV file is not beside the script.
data_dir <- NULL
output_root <- file.path(project_dir, paste0("HuNoV_LRPF_single_", inoculum_id))

# These settings match the locked GZJY122 application in the manuscript.
cfg <- list(
  trajectory_window = c(0, 7),
  auc_window = c(0, 5),
  reference_dpi = 0,
  bs_degree = 3L,
  bs_intercept = FALSE,
  feature_grid_points = 300L,
  grid_points = c(300L, 1001L, 5001L),
  n_sim = 1000L,
  seed = 123L,
  qc = list(
    min_batches = 4L,
    min_timepoints = 6L,
    require_dpi0 = TRUE,
    min_followup_dpi = 6,
    require_all_within_boundary = TRUE,
    allow_duplicate_dpi = FALSE,
    early_window = c(0, 1),
    middle_window = c(2, 4),
    late_window = c(5, 7),
    min_early_timepoints = 1L,
    min_middle_timepoints = 2L,
    min_late_timepoints = 1L
  )
)

locked_basis <- make_locked_basis_info(
  df = 4L,
  internal_knots = 3,
  boundary_knots = cfg$trajectory_window,
  degree = cfg$bs_degree,
  intercept = cfg$bs_intercept
)

# =============================================================================
# 2. Input and locked analysis
# =============================================================================

data_dir <- resolve_data_dir(input_file, project_dir, data_dir)
input_table <- tibble::tibble(
  inoculum_id = inoculum_id,
  genotype_ptype = genotype_ptype,
  file = input_file,
  model_key = inoculum_id
)
input_data <- combine_inocula(input_table, data_dir)

result <- run_locked_single_analysis(
  data = input_data,
  basis_info = locked_basis,
  cfg = cfg,
  output_root = output_root,
  prefix = inoculum_id
)

dirs <- make_output_dirs(output_root)

# =============================================================================
# 3. Compact publication-style outputs
# =============================================================================

application_composite <- plot_locked_application_composite(
  result$data,
  result$prediction,
  result$lobo,
  result$diagnostics
)
save_tiff(
  application_composite,
  file.path(dirs$plots, paste0(inoculum_id, "_fixed_settings_composite.tiff")),
  width = 10.5,
  height = 8.0
)

feature_metrics <- c(
  "peak_logVL",
  "peak_dpi",
  "AUC_0_5",
  "max_positive_slope",
  "max_positive_slope_dpi",
  "most_negative_slope",
  "most_negative_slope_dpi",
  "fitted_0dpi_logVL",
  "max_net_increase",
  "net_peak_dpi",
  "centered_AUC_0_5"
)

feature_table_out <- result$features |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(feature_metrics),
    names_to = "metric",
    values_to = "point_estimate"
  ) |>
  dplyr::left_join(
    result$intervals |>
      dplyr::select(
        model_key,
        metric,
        median,
        lower_2.5,
        upper_97.5,
        valid_draws = n_valid
      ),
    by = c("model_key", "metric")
  ) |>
  dplyr::arrange(metric)

write_csv(
  feature_table_out,
  file.path(dirs$tables, paste0(inoculum_id, "_feature_summary.csv"))
)

# =============================================================================
# 4. Reproducibility records
# =============================================================================

q <- cfg$qc
metadata <- tibble::tibble(
  item = c(
    "analysis_version",
    "analysis_role",
    "inoculum_id",
    "genotype_ptype",
    "input_file",
    "spline_df",
    "internal_knots",
    "boundary_knots",
    "bs_degree",
    "bs_intercept",
    "trajectory_window",
    "auc_window",
    "reference_dpi",
    "feature_grid_points",
    "grid_stability_points",
    "simulation_draws",
    "simulation_seed",
    "qc_min_batches",
    "qc_min_timepoints",
    "qc_require_dpi0",
    "qc_min_followup_dpi",
    "qc_require_all_within_boundary",
    "qc_allow_duplicate_dpi",
    "qc_early_window",
    "qc_middle_window",
    "qc_late_window",
    "qc_min_early_timepoints",
    "qc_min_middle_timepoints",
    "qc_min_late_timepoints",
    "eligible_observations",
    "eligible_batches",
    "R_version",
    "run_timestamp"
  ),
  value = c(
    "10.3-single-locked-specification",
    "Locked-specification application; no model or spline reselection",
    inoculum_id,
    genotype_ptype,
    basename(input_file),
    as.character(locked_basis$spline_df),
    paste(locked_basis$knots, collapse = ";"),
    paste(locked_basis$boundary_knots, collapse = ";"),
    as.character(locked_basis$degree),
    as.character(locked_basis$intercept),
    paste(cfg$trajectory_window, collapse = ";"),
    paste(cfg$auc_window, collapse = ";"),
    as.character(cfg$reference_dpi),
    as.character(cfg$feature_grid_points),
    paste(cfg$grid_points, collapse = ";"),
    as.character(cfg$n_sim),
    as.character(cfg$seed + 101L),
    as.character(q$min_batches),
    as.character(q$min_timepoints),
    as.character(q$require_dpi0),
    as.character(q$min_followup_dpi),
    as.character(q$require_all_within_boundary),
    as.character(q$allow_duplicate_dpi),
    paste(q$early_window, collapse = ";"),
    paste(q$middle_window, collapse = ";"),
    paste(q$late_window, collapse = ";"),
    as.character(q$min_early_timepoints),
    as.character(q$min_middle_timepoints),
    as.character(q$min_late_timepoints),
    as.character(nrow(result$data)),
    as.character(dplyr::n_distinct(result$data$batch)),
    R.version.string,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
)

specification_table <- tibble::tibble(
  component = c(
    "Spline df",
    "Internal knot",
    "Boundary knots",
    "Degree",
    "Intercept",
    "AUC window",
    "Reference dpi"
  ),
  locked_value = c(
    as.character(locked_basis$spline_df),
    paste(locked_basis$knots, collapse = ";"),
    paste(locked_basis$boundary_knots, collapse = ";"),
    as.character(locked_basis$degree),
    as.character(locked_basis$intercept),
    paste(cfg$auc_window, collapse = ";"),
    as.character(cfg$reference_dpi)
  ),
  reselected_from_new_data = FALSE
)

write_csv(metadata, file.path(output_root, "analysis_metadata.csv"))
write_csv(specification_table, file.path(output_root, "locked_specification.csv"))
write_csv(
  file_checksum_table(file.path(data_dir, input_file), labels = input_file),
  file.path(output_root, "input_file_checksum.csv")
)
write_model_summary(
  result$bundle$model,
  file.path(output_root, paste0(inoculum_id, "_model_summary.txt"))
)
write_session_info(file.path(output_root, "sessionInfo.txt"))
write_output_manifest(output_root)

message(
  "Locked-specification single-inoculum analysis completed: ",
  normalizePath(output_root, winslash = "/", mustWork = TRUE)
)
