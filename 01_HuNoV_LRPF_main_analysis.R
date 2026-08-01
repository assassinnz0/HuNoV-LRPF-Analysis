# HuNoV-LRPF publication analysis
#
# This script reproduces the complete manuscript analysis:
#   1. development on the five-inoculum panel;
#   2. locked application to GZJY122;
#   3. fixed-specification six-inoculum update;
#   4. two prespecified GZJY122 sensitivity analyses;
#   5. manuscript and supplementary figures, source tables, model objects,
#      checksums, run metadata, and session information.
#
# Required input files:
#   df_315.csv, df_324.csv, df_17.csv, df_19.csv, df_29.csv, df_122.csv
# Each file must contain the columns batch, dpi, and VL.

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

cfg <- list(
  trajectory_window = c(0, 7),
  auc_window = c(0, 5),
  reference_dpi = 0,
  bs_degree = 3L,
  bs_intercept = FALSE,
  candidate_df = c(3L, 4L),
  bic_tolerance = 2,
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

# Set data_dir to an explicit folder when the CSV files are not beside the script.
data_dir <- NULL
output_root <- file.path(project_dir, paste0("HuNoV_LRPF_v", gsub("\\.", "_", hunov_analysis_version), "_full_results"))

# The five-inoculum development dataset combines mostly previously reported records with additional records incorporated in this study.
development_files <- tibble::tribble(
  ~inoculum_id, ~genotype_ptype, ~file,        ~model_key,
  "GZJY315",   "GII.2[P16]",   "df_315.csv",  "GZJY315",
  "GZJY324",   "GII.17[P17]",  "df_324.csv",  "GZJY324",
  "GZJY17",    "GII.4[P31]",   "df_17.csv", "GZJY17",
  "GZJY19",    "GII.4[P16]",   "df_19.csv",   "GZJY19",
  "GZJY29",    "GII.3[P12]",   "df_29.csv",   "GZJY29"
)

# GZJY122 is analyzed under the fixed df = 4, knot = 3 dpi specification.
application_file <- tibble::tibble(
  inoculum_id = "GZJY122",
  genotype_ptype = "GII.4[P31]",
  file = "df_122.csv",
  model_key = "GZJY122"
)

required_files <- c(development_files$file, application_file$file)
data_dir <- resolve_data_dir(required_files, project_dir, data_dir)
dirs <- make_publication_output_dirs(output_root)

input_paths <- file.path(data_dir, required_files)
write_csv(
  file_checksum_table(input_paths, labels = required_files),
  file.path(dirs$logs, "input_file_checksums.csv")
)

# =============================================================================
# 2. Development dataset and structural quality control
# =============================================================================

development_raw <- combine_inocula(development_files, data_dir)
development_qc <- structural_qc(development_raw, cfg)
development_filtered <- apply_structural_qc(development_raw, development_qc, cfg)
development_data <- development_filtered$data

legacy_outcome_flags <- development_qc |>
  dplyr::transmute(
    model_key,
    inoculum_id,
    genotype_ptype,
    batch,
    observed_0dpi_logVL,
    observed_peak_logVL,
    observed_peak_dpi,
    observed_peak_increase,
    peak_increase_at_least_2_log10 = observed_peak_increase >= 2,
    observed_peak_after_1_dpi = observed_peak_dpi > 1,
    note = "Descriptive only; not used for structural QC"
  )

retained_summary <- development_data |>
  dplyr::group_by(model_key, inoculum_id, genotype_ptype) |>
  dplyr::summarise(
    observations = dplyr::n(),
    valid_batches = dplyr::n_distinct(batch),
    dpi_min = min(dpi),
    dpi_max = max(dpi),
    .groups = "drop"
  )

# =============================================================================
# 3. Development-stage model hierarchy and spline selection
# =============================================================================

model_bundles <- list()
for (spline_df in cfg$candidate_df) {
  basis <- make_basis_info_df(development_data$dpi, spline_df, cfg)
  for (model_type in c("M0", "M1", "M2")) {
    key <- paste0(model_type, "_df", spline_df)
    model_bundles[[key]] <- fit_bundle(
      development_data,
      model_type = model_type,
      basis_info = basis,
      group_col = "inoculum_batch"
    )
  }
}

model_hierarchy <- purrr::map_dfr(model_bundles, model_row) |>
  dplyr::mutate(
    delta_AIC_global = AIC - min(AIC),
    delta_BIC_global = BIC - min(BIC)
  ) |>
  dplyr::arrange(spline_df, model_type)

selected_df <- select_df_by_bic(model_hierarchy, model_type = "M2", tolerance = cfg$bic_tolerance)
primary_bundle <- model_bundles[[paste0("M2_df", selected_df)]]
selected_basis <- primary_bundle$basis

m2_df_selection <- model_hierarchy |>
  dplyr::filter(model_type == "M2") |>
  dplyr::mutate(
    delta_BIC_within_M2 = BIC - min(BIC),
    selection = ifelse(spline_df == selected_df, "Selected within M2", "Not selected")
  ) |>
  dplyr::arrange(BIC)

fixed_df_comparison <- model_hierarchy |>
  dplyr::filter(spline_df == selected_df, model_type %in% c("M1", "M2")) |>
  dplyr::mutate(
    role = dplyr::recode(
      model_type,
      M1 = "Shared-shape comparator",
      M2 = "Feature-extraction model"
    )
  ) |>
  dplyr::arrange(model_type)

lrt_table <- function(reduced, full, label) {
  out <- stats::anova(reduced$model, full$model)
  row <- as.data.frame(out)[2, , drop = FALSE]
  tibble::tibble(
    comparison = label,
    spline_df = full$basis$spline_df,
    chi_square = row$Chisq,
    df = row$`Chi Df`,
    p_value = row$`Pr(>Chisq)`
  )
}

likelihood_ratio_tests <- dplyr::bind_rows(
  lrt_table(model_bundles[[paste0("M0_df", selected_df)]], model_bundles[[paste0("M1_df", selected_df)]], "M0 vs M1"),
  lrt_table(model_bundles[[paste0("M1_df", selected_df)]], model_bundles[[paste0("M2_df", selected_df)]], "M1 vs M2")
)

# =============================================================================
# 4. Development-stage feature extraction and uncertainty
# =============================================================================

development_dpi_grid <- make_dpi_grid(cfg$feature_grid_points, cfg$trajectory_window, c(cfg$auc_window, cfg$reference_dpi))
development_grid <- new_fixed_grid(development_data, development_dpi_grid, selected_basis)
development_prediction <- predict_fixed(primary_bundle, development_grid)
development_features <- feature_table(development_prediction, cfg)
development_draws <- simulate_feature_draws(primary_bundle, development_grid, cfg, cfg$seed)
development_intervals <- summarise_feature_draws(development_draws)
development_pairwise <- pairwise_feature_contrasts(development_draws)
development_derivatives <- derivative_curve(development_prediction)
development_grid_check <- grid_stability(primary_bundle, development_data, cfg)

development_feature_table <- development_features |>
  tidyr::pivot_longer(
    cols = c(
      peak_logVL, peak_dpi, AUC_0_5,
      max_positive_slope, max_positive_slope_dpi,
      most_negative_slope, most_negative_slope_dpi,
      fitted_0dpi_logVL, max_net_increase, net_peak_dpi,
      centered_AUC_0_5
    ),
    names_to = "metric",
    values_to = "point_estimate"
  ) |>
  dplyr::left_join(
    development_intervals |>
      dplyr::select(
        model_key, metric, median, lower_2.5, upper_97.5,
        valid_draws = n_valid
      ),
    by = c("model_key", "metric")
  ) |>
  dplyr::arrange(model_key, metric)

# =============================================================================
# 5. Development-stage LOBO prediction and diagnostics
# =============================================================================

development_lobo_m1 <- run_lobo(
  development_data,
  model_type = "M1",
  basis_info = selected_basis,
  group_col = "inoculum_batch",
  fold_col = "inoculum_batch"
)
development_lobo_m2 <- run_lobo(
  development_data,
  model_type = "M2",
  basis_info = selected_basis,
  group_col = "inoculum_batch",
  fold_col = "inoculum_batch"
)

development_lobo_summary <- dplyr::bind_rows(
  dplyr::mutate(summarise_lobo(development_lobo_m1), model_type = "M1"),
  dplyr::mutate(summarise_lobo(development_lobo_m2), model_type = "M2")
)

development_lobo_by_inoculum <- development_lobo_m2 |>
  dplyr::filter(status == "OK") |>
  dplyr::group_by(model_key, inoculum_id, genotype_ptype) |>
  dplyr::summarise(
    test_observations = dplyr::n(),
    folds = dplyr::n_distinct(fold),
    RMSE = sqrt(mean((logVL - predicted)^2)),
    MAE = mean(abs(logVL - predicted)),
    bias = mean(logVL - predicted),
    correlation = stats::cor(logVL, predicted),
    prediction_interval_coverage = mean(logVL >= prediction_lower & logVL <= prediction_upper),
    .groups = "drop"
  )

development_diagnostics <- residual_diagnostics(primary_bundle)
development_reconstruction <- reconstruction_summary(development_diagnostics)
development_random_intercepts <- random_intercept_blups(primary_bundle)
development_model_status <- model_status(primary_bundle)

summarise_lobo_folds <- function(predictions, suffix) {
  predictions |>
    dplyr::filter(status == "OK") |>
    dplyr::group_by(fold, model_key, inoculum_id, genotype_ptype) |>
    dplyr::summarise(
      n = dplyr::n(),
      RMSE = sqrt(mean((logVL - predicted)^2)),
      MAE = mean(abs(logVL - predicted)),
      PI_coverage = mean(logVL >= prediction_lower & logVL <= prediction_upper),
      singular_training_fit = any(singular %in% TRUE),
      convergence_message = paste(
        unique(convergence_message[!is.na(convergence_message) & nzchar(convergence_message)]),
        collapse = " | "
      ),
      .groups = "drop"
    ) |>
    dplyr::rename_with(
      ~paste0(.x, "_", suffix),
      c(
        "RMSE", "MAE", "PI_coverage",
        "singular_training_fit", "convergence_message"
      )
    )
}

development_lobo_fold_comparison <- summarise_lobo_folds(development_lobo_m1, "M1") |>
  dplyr::left_join(
    summarise_lobo_folds(development_lobo_m2, "M2"),
    by = c("fold", "model_key", "inoculum_id", "genotype_ptype", "n")
  ) |>
  dplyr::mutate(
    delta_RMSE_M2_minus_M1 = RMSE_M2 - RMSE_M1,
    delta_MAE_M2_minus_M1 = MAE_M2 - MAE_M1
  )

development_residual_by_inoculum <- development_diagnostics |>
  dplyr::group_by(model_key, inoculum_id, genotype_ptype) |>
  dplyr::summarise(
    n = dplyr::n(),
    batches = dplyr::n_distinct(batch),
    residual_SD = stats::sd(residual),
    .groups = "drop"
  )

development_reconstruction_publication <- development_residual_by_inoculum |>
  dplyr::left_join(
    development_lobo_by_inoculum |>
      dplyr::select(
        model_key, inoculum_id, genotype_ptype,
        LOBO_RMSE = RMSE,
        LOBO_MAE = MAE,
        LOBO_r = correlation,
        LOBO_PI_coverage = prediction_interval_coverage
      ),
    by = c("model_key", "inoculum_id", "genotype_ptype")
  )

development_lobo_m2_overall <- dplyr::filter(
  development_lobo_summary,
  model_type == "M2"
)

development_reconstruction_publication <- dplyr::bind_rows(
  development_reconstruction_publication,
  tibble::tibble(
    model_key = "Overall",
    inoculum_id = "Overall",
    genotype_ptype = "All tested inocula",
    n = nrow(development_data),
    batches = dplyr::n_distinct(development_data$inoculum_batch),
    residual_SD = development_model_status$residual_sd[1],
    LOBO_RMSE = development_lobo_m2_overall$RMSE[1],
    LOBO_MAE = development_lobo_m2_overall$MAE[1],
    LOBO_r = development_lobo_m2_overall$correlation[1],
    LOBO_PI_coverage = development_lobo_m2_overall$prediction_interval_coverage[1]
  )
)

development_r2 <- tryCatch(
  performance::r2_nakagawa(primary_bundle$model),
  error = function(e) NULL
)
model_explanatory_metrics <- tibble::tibble(
  metric = c(
    "Marginal R2",
    "Conditional R2",
    "Random-intercept variance",
    "Residual variance",
    "Batch ICC"
  ),
  value = c(
    if (is.null(development_r2)) NA_real_ else as.numeric(development_r2$R2_marginal),
    if (is.null(development_r2)) NA_real_ else as.numeric(development_r2$R2_conditional),
    development_model_status$random_intercept_variance,
    development_model_status$residual_sd^2,
    development_model_status$random_intercept_variance /
      (development_model_status$random_intercept_variance + development_model_status$residual_sd^2)
  )
)

# =============================================================================
# 6. Locked application to GZJY122
# =============================================================================

locked_basis <- make_locked_basis_info(
  df = 4L,
  internal_knots = 3,
  boundary_knots = cfg$trajectory_window,
  degree = cfg$bs_degree,
  intercept = cfg$bs_intercept
)

specification_check <- tibble::tibble(
  component = c("Spline df", "Internal knot", "Boundary knots", "Degree", "Intercept"),
  development_value = c(
    as.character(selected_basis$spline_df),
    paste(selected_basis$knots, collapse = ";"),
    paste(selected_basis$boundary_knots, collapse = ";"),
    as.character(selected_basis$degree),
    as.character(selected_basis$intercept)
  ),
  locked_value = c(
    as.character(locked_basis$spline_df),
    paste(locked_basis$knots, collapse = ";"),
    paste(locked_basis$boundary_knots, collapse = ";"),
    as.character(locked_basis$degree),
    as.character(locked_basis$intercept)
  )
) |>
  dplyr::mutate(matches = development_value == locked_value)

if (!all(specification_check$matches)) {
  stop("The development fit does not reproduce the fixed df = 4, knot = 3 dpi specification.", call. = FALSE)
}

application_raw <- combine_inocula(application_file, data_dir)
application_result <- run_locked_single_analysis(
  data = application_raw,
  basis_info = locked_basis,
  cfg = cfg,
  output_root = file.path(output_root, "application_intermediate"),
  prefix = "GZJY122"
)

application_feature_table <- application_result$features |>
  tidyr::pivot_longer(
    cols = c(
      peak_logVL, peak_dpi, AUC_0_5, max_positive_slope,
      most_negative_slope, fitted_0dpi_logVL, max_net_increase,
      centered_AUC_0_5
    ),
    names_to = "metric",
    values_to = "point_estimate"
  ) |>
  dplyr::left_join(
    application_result$intervals |>
      dplyr::select(
        model_key, metric, median, lower_2.5, upper_97.5,
        valid_draws = n_valid
      ),
    by = c("model_key", "metric")
  )


application_lobo_folds <- application_result$lobo |>
  dplyr::filter(status == "OK") |>
  dplyr::group_by(fold) |>
  dplyr::summarise(
    status = dplyr::first(status),
    singular_training_fit = any(singular),
    n = dplyr::n(),
    RMSE = sqrt(mean((logVL - predicted)^2)),
    MAE = mean(abs(logVL - predicted)),
    bias = mean(logVL - predicted),
    PI_coverage = mean(logVL >= prediction_lower & logVL <= prediction_upper),
    .groups = "drop"
  )

application_model_and_features <- dplyr::bind_rows(
  tibble::tibble(
    category = "Model status",
    metric = c("Residual SD", "Random-intercept variance", "Singular fit"),
    point_or_status = c(
      format(application_result$model_status$residual_sd[1], digits = 6),
      format(application_result$model_status$random_intercept_variance[1], digits = 6),
      ifelse(application_result$model_status$singular[1], "Yes", "No")
    ),
    lower_2.5 = NA_real_,
    upper_97.5 = NA_real_,
    valid_draws = NA_integer_
  ),
  application_feature_table |>
    dplyr::transmute(
      category = dplyr::if_else(
        metric %in% c("fitted_0dpi_logVL", "max_net_increase", "centered_AUC_0_5"),
        "Reference-centered feature",
        "Primary feature"
      ),
      metric,
      point_or_status = format(point_estimate, digits = 8),
      lower_2.5,
      upper_97.5,
      valid_draws = valid_draws
    )
)

# =============================================================================
# 6A. GZJY122 sensitivity analyses
# =============================================================================

# These analyses do not reselect the spline specification and do not replace
# the primary mixed-effects model. They assess whether the reported trajectory
# features are materially affected by (i) removing the unresolved random
# intercept and (ii) omitting one batch from the feature-estimation dataset.

application_sensitivity_metrics <- c(
  "peak_logVL",
  "peak_dpi",
  "AUC_0_5",
  "max_positive_slope",
  "most_negative_slope",
  "fitted_0dpi_logVL",
  "max_net_increase",
  "centered_AUC_0_5"
)

application_sensitivity_labels <- c(
  peak_logVL = "Peak log10 RNA load",
  peak_dpi = "Peak time",
  AUC_0_5 = "0–5-dpi log-scale AUC",
  max_positive_slope = "Maximum positive slope",
  most_negative_slope = "Most negative slope",
  fitted_0dpi_logVL = "Fitted 0-dpi reference level",
  max_net_increase = "Maximum increase above 0-dpi reference",
  centered_AUC_0_5 = "0–5-dpi reference-centered AUC"
)

application_sensitivity_units <- c(
  peak_logVL = "log10 copies/µL",
  peak_dpi = "dpi",
  AUC_0_5 = "log10(copies/µL)·day",
  max_positive_slope = "log10(copies/µL)/day",
  most_negative_slope = "log10(copies/µL)/day",
  fitted_0dpi_logVL = "log10 copies/µL",
  max_net_increase = "log10 units",
  centered_AUC_0_5 = "log10(copies/µL)·day"
)

application_sensitivity_dpi_grid <- sort(unique(application_result$prediction$dpi))

# Sensitivity analysis 1: retain the locked spline basis but remove the batch
# random intercept. This ordinary least-squares fit is descriptive only.
application_fixed_only_data <- add_basis(application_result$data, locked_basis)
application_fixed_only_formula <- stats::as.formula(
  paste0(
    "logVL ~ ",
    paste(locked_basis$columns, collapse = " + ")
  )
)
application_fixed_only_model <- stats::lm(
  application_fixed_only_formula,
  data = application_fixed_only_data
)

if (anyNA(stats::coef(application_fixed_only_model))) {
  stop(
    "The fixed-effect-only sensitivity model contains non-estimable coefficients.",
    call. = FALSE
  )
}

application_fixed_only_grid <- new_fixed_grid(
  application_result$data,
  application_sensitivity_dpi_grid,
  locked_basis
)
application_fixed_only_predict <- stats::predict(
  application_fixed_only_model,
  newdata = application_fixed_only_grid,
  se.fit = TRUE
)
application_fixed_only_prediction <- application_fixed_only_grid |>
  dplyr::mutate(
    predicted = as.numeric(application_fixed_only_predict$fit),
    fixed_se = as.numeric(application_fixed_only_predict$se.fit),
    lower = predicted - 1.96 * fixed_se,
    upper = predicted + 1.96 * fixed_se
  )
application_fixed_only_features <- feature_table(
  application_fixed_only_prediction,
  cfg
)

application_fixed_only_status <- tibble::tibble(
  model = "Fixed-effect-only sensitivity model",
  formula = paste(
    deparse(stats::formula(application_fixed_only_model)),
    collapse = " "
  ),
  observations = stats::nobs(application_fixed_only_model),
  residual_sd = stats::sigma(application_fixed_only_model),
  AIC = stats::AIC(application_fixed_only_model),
  BIC = stats::BIC(application_fixed_only_model),
  role = "Sensitivity analysis only; no model or spline reselection"
)

# Sensitivity analysis 2: recreate the same five four-batch training fits used
# by the existing LOBO procedure and extract trajectory features from each fit.
# The current run_lobo() helper returns predictions rather than fitted bundles,
# so these training fits are reconstructed with the identical locked basis,
# formula, optimizer sequence, and held-out batches.
application_sensitivity_folds <- unique(as.character(application_result$lobo$fold))

application_lobo_feature_fit_objects <- stats::setNames(
  purrr::map(application_sensitivity_folds, function(fold) {
    train <- droplevels(
      application_result$data[
        as.character(application_result$data$batch) != fold,
        ,
        drop = FALSE
      ]
    )

    tryCatch(
      {
        bundle <- fit_bundle(
          train,
          model_type = "single",
          basis_info = locked_basis,
          group_col = "batch"
        )
        grid <- new_fixed_grid(
          train,
          application_sensitivity_dpi_grid,
          locked_basis
        )
        prediction <- predict_fixed(bundle, grid)
        features <- feature_table(prediction, cfg)

        list(
          status = "OK",
          error_message = "",
          bundle = bundle,
          prediction = prediction,
          features = features
        )
      },
      error = function(e) {
        list(
          status = "fit_failed",
          error_message = conditionMessage(e),
          bundle = NULL,
          prediction = NULL,
          features = NULL
        )
      }
    )
  }),
  application_sensitivity_folds
)

application_lobo_training_feature_estimates <- purrr::imap_dfr(
  application_lobo_feature_fit_objects,
  function(fit_object, fold) {
    if (!identical(fit_object$status, "OK")) {
      return(
        tibble::tibble(
          held_out_batch = fold,
          status = fit_object$status,
          error_message = fit_object$error_message,
          training_batches = NA_integer_,
          singular_training_fit = NA,
          convergence_message = NA_character_,
          metric = application_sensitivity_metrics,
          estimate = NA_real_
        )
      )
    }

    fit_status <- model_status(fit_object$bundle)
    
    fit_object$features |>
      dplyr::select(dplyr::all_of(application_sensitivity_metrics)) |>
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = "metric",
        values_to = "estimate"
      ) |>
      dplyr::mutate(
        held_out_batch = fold,
        status = "OK",
        error_message = "",
        training_batches = dplyr::n_distinct(
          fit_object$bundle$data$batch
        ),
        singular_training_fit = fit_status$singular[1],
        convergence_message = fit_status$convergence_message[1],
        .before = 1
      )
  }
)

safe_finite_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}

safe_finite_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

application_primary_sensitivity_long <- application_result$features |>
  dplyr::select(dplyr::all_of(application_sensitivity_metrics)) |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "metric",
    values_to = "primary_mixed_model"
  )

application_fixed_only_sensitivity_long <- application_fixed_only_features |>
  dplyr::select(dplyr::all_of(application_sensitivity_metrics)) |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "metric",
    values_to = "fixed_effect_only_model"
  )

application_lobo_feature_range <- application_lobo_training_feature_estimates |>
  dplyr::left_join(
    application_primary_sensitivity_long,
    by = "metric"
  ) |>
  dplyr::group_by(metric) |>
  dplyr::summarise(
    successful_refits = sum(status == "OK" & is.finite(estimate)),
    four_batch_min = safe_finite_min(estimate),
    four_batch_max = safe_finite_max(estimate),
    maximum_absolute_deviation = safe_finite_max(
      abs(estimate - primary_mixed_model)
    ),
    .groups = "drop"
  )

application_sensitivity_summary <- application_primary_sensitivity_long |>
  dplyr::left_join(
    application_fixed_only_sensitivity_long,
    by = "metric"
  ) |>
  dplyr::left_join(
    application_lobo_feature_range,
    by = "metric"
  ) |>
  dplyr::mutate(
    feature = unname(application_sensitivity_labels[metric]),
    unit = unname(application_sensitivity_units[metric]),
    fixed_only_absolute_difference = abs(
      fixed_effect_only_model - primary_mixed_model
    ),
    metric_order = match(metric, application_sensitivity_metrics)
  ) |>
  dplyr::arrange(metric_order)

# =============================================================================
# 7. Fixed-specification six-inoculum update
# =============================================================================

six_data <- dplyr::bind_rows(development_data, application_result$data) |>
  dplyr::mutate(
    model_key = factor(as.character(model_key), levels = unique(as.character(model_key))),
    inoculum_id = factor(as.character(inoculum_id), levels = unique(as.character(inoculum_id))),
    batch = factor(as.character(batch)),
    inoculum_batch = interaction(model_key, batch, drop = TRUE, sep = "__batch__")
  )

six_bundle <- fit_bundle(six_data, "M2", locked_basis, "inoculum_batch")
six_grid <- new_fixed_grid(six_data, development_dpi_grid, locked_basis)
six_prediction <- predict_fixed(six_bundle, six_grid)
six_features <- feature_table(six_prediction, cfg)
six_draws <- simulate_feature_draws(six_bundle, six_grid, cfg, cfg$seed + 6122L)
six_intervals <- summarise_feature_draws(six_draws)
six_model_status <- model_status(six_bundle)

six_feature_long <- six_features |>
  tidyr::pivot_longer(
    cols = c(
      peak_logVL, peak_dpi, AUC_0_5, max_positive_slope,
      most_negative_slope, centered_AUC_0_5
    ),
    names_to = "metric",
    values_to = "point_estimate"
  ) |>
  dplyr::left_join(
    six_intervals |>
      dplyr::select(model_key, metric, lower_2.5, upper_97.5),
    by = c("model_key", "metric")
  ) |>
  dplyr::mutate(
    formatted = sprintf(
      "%.2f (%.2f to %.2f)",
      point_estimate, lower_2.5, upper_97.5
    )
  )

six_publication_features <- six_feature_long |>
  dplyr::select(model_key, inoculum_id, genotype_ptype, metric, formatted) |>
  tidyr::pivot_wider(names_from = metric, values_from = formatted) |>
  dplyr::mutate(
    dataset_role = dplyr::if_else(
      as.character(model_key) == "GZJY122",
      "Locked application",
      "Development"
    )
  ) |>
  dplyr::select(
    inoculum_id, genotype_ptype, dataset_role,
    peak_logVL, peak_dpi, AUC_0_5,
    max_positive_slope, most_negative_slope,
    centered_AUC_0_5
  )

six_contrasts <- pairwise_feature_contrasts(
  dplyr::filter(six_draws, as.character(model_key) %in% c("GZJY122", "GZJY17"))
) |>
  dplyr::filter(
    (inoculum_A == "GZJY122" & inoculum_B == "GZJY17") |
      (inoculum_A == "GZJY17" & inoculum_B == "GZJY122")
  ) |>
  dplyr::mutate(
    median_difference = ifelse(inoculum_A == "GZJY122", median_difference, -median_difference),
    lower_original = lower_2.5,
    upper_original = upper_97.5,
    lower_2.5 = ifelse(inoculum_A == "GZJY122", lower_original, -upper_original),
    upper_97.5 = ifelse(inoculum_A == "GZJY122", upper_original, -lower_original),
    proportion_above_zero = ifelse(inoculum_A == "GZJY122", proportion_above_zero, 1 - proportion_above_zero),
    contrast = "GZJY122 minus GZJY17"
  ) |>
  dplyr::select(contrast, metric, median_difference, lower_2.5, upper_97.5, proportion_above_zero)

six_difference_curve <- fixed_curve_difference(six_bundle, six_grid, "GZJY122", "GZJY17")

# =============================================================================
# 8. Manuscript and supplementary figures
# =============================================================================

main_figures <- list(
  Figure_2_development_data_and_basis = list(
    plot = plot_development_design(development_data, selected_basis, cfg$trajectory_window),
    width = 11.5,
    height = 8.5
  ),
  Figure_3_development_trajectories = list(
    plot = plot_trajectory_facets_publication(development_data, development_prediction),
    width = 10.5,
    height = 7.5
  ),
  Figure_4_development_trajectory_features = list(
    plot = plot_endpoint_composite(development_features),
    width = 12.5,
    height = 10.2
  ),
  Figure_5_GZJY122_locked_application = list(
    plot = plot_locked_application_composite(
      application_result$data,
      application_result$prediction,
      application_result$lobo,
      application_result$diagnostics
    ),
    width = 10.5,
    height = 8.0
  ),
  Figure_6_GZJY122_minus_GZJY17_difference = list(
    plot = plot_curve_difference(six_difference_curve),
    width = 7.5,
    height = 5.2
  )
)

supplementary_figures <- list(
  Figure_S1_development_LOBO = list(
    plot = plot_lobo_combined(development_lobo_m2),
    width = 10.5,
    height = 5.2
  ),
  Figure_S2_development_pooled_diagnostics = list(
    plot = plot_pooled_diagnostics(development_diagnostics),
    width = 10.5,
    height = 8.0
  ),
  Figure_S3_development_inoculum_diagnostics = list(
    plot = plot_inoculum_diagnostics(development_diagnostics),
    width = 12.0,
    height = 11.5
  ),
  Figure_S4_development_derivative_curves = list(
    plot = plot_derivative_curves_publication(development_derivatives),
    width = 9.0,
    height = 5.5
  ),
  Figure_S5_development_centered_trajectories = list(
    plot = plot_centered_trajectories(development_prediction, cfg$reference_dpi),
    width = 10.5,
    height = 7.5
  ),
  Figure_S6_six_inoculum_fixed_specification = list(
    plot = plot_trajectory_facets_publication(six_data, six_prediction),
    width = 10.5,
    height = 8.0
  )
)

save_figure_set(main_figures, dirs$main_figures)
save_figure_set(supplementary_figures, dirs$supplementary_figures)

# =============================================================================
# 9. Manuscript-aligned tables and source data
# =============================================================================

format_fixed <- function(x, digits) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

format_percent <- function(x, digits = 1L) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f%%"), 100 * x))
}

format_interval <- function(point, lower, upper, digits = 2L, separator = "–") {
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f", separator, "%.", digits, "f)"),
    point, lower, upper
  )
}

format_interval_to <- function(point, lower, upper, digits = 2L) {
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f to %.", digits, "f)"),
    point, lower, upper
  )
}

development_order <- as.character(development_files$inoculum_id)
six_order <- c(development_order, "GZJY122")

primary_metrics <- c(
  "peak_logVL", "peak_dpi", "AUC_0_5",
  "max_positive_slope", "most_negative_slope"
)
primary_metric_labels <- c(
  peak_logVL = "Peak_logVL",
  peak_dpi = "Peak_dpi",
  AUC_0_5 = "AUC",
  max_positive_slope = "Max_Rate",
  most_negative_slope = "Max_Decline_Rate"
)
centered_metrics <- c(
  "centered_AUC_0_5", "fitted_0dpi_logVL",
  "max_net_increase", "net_peak_dpi"
)
centered_metric_labels <- c(
  centered_AUC_0_5 = "Baseline_Adjusted_AUC",
  fitted_0dpi_logVL = "Baseline_logVL",
  max_net_increase = "Max_Net_Increase",
  net_peak_dpi = "Net_Peak_dpi"
)

# Main Table 1: dataset roles and retained structure.
table_1_dataset_roles <- dplyr::bind_rows(
  development_files |>
    dplyr::transmute(
      `Analysis role` = "Development",
      `Inoculum(s)` = inoculum_id,
      Genotype = genotype_ptype,
      `No. of inocula` = nrow(development_files),
      `Valid batches` = dplyr::n_distinct(development_data$inoculum_batch),
      Observations = nrow(development_data),
      Purpose = "Model and workflow development"
    ),
  application_file |>
    dplyr::transmute(
      `Analysis role` = "Locked application",
      `Inoculum(s)` = inoculum_id,
      Genotype = genotype_ptype,
      `No. of inocula` = 1L,
      `Valid batches` = dplyr::n_distinct(application_result$data$batch),
      Observations = nrow(application_result$data),
      Purpose = "Locked application"
    )
)

# Main Table 2 is static experimental metadata included so the exported table
# set matches the manuscript's complete Table 1–6 sequence.
table_2_primers_probes <- tibble::tribble(
  ~`Norovirus-GII`, ~Name, ~`Sequence (5′–3′)`, ~Reference,
  "Forward", "QNIF2d", "ATGTTCAGRTGGATGAGRTTCTCWGA", "Stals et al. (2009)",
  "Reverse", "COG2R", "TCGACGCCATCTTCATTCACA", "Kageyama et al. (2003)",
  "Probe", "QNIFS", "FAM-AGCACGTGGGAGGGCGATCG-TAMRA", "Stals et al. (2009)"
)

# Main Table 3: model specification and supporting comparator.
table_3_model_specification <- dplyr::bind_rows(
  m2_df_selection |>
    dplyr::transmute(
      Model = paste0(model_type, ", df = ", spline_df),
      `Internal knot` = dplyr::if_else(
        is.na(internal_knots) | internal_knots %in% c("", "None"),
        "None",
        paste0(internal_knots, " dpi")
      ),
      AIC = format_fixed(AIC, 2),
      BIC = format_fixed(BIC, 2),
      `Residual SD` = format_fixed(residual_sd, 3),
      Role = selection
    ),
  fixed_df_comparison |>
    dplyr::filter(model_type == "M1") |>
    dplyr::transmute(
      Model = paste0(model_type, ", df = ", spline_df),
      `Internal knot` = paste0(paste(selected_basis$knots, collapse = ";"), " dpi"),
      AIC = format_fixed(AIC, 2),
      BIC = format_fixed(BIC, 2),
      `Residual SD` = format_fixed(residual_sd, 3),
      Role = "Compact shared-shape comparator"
    )
) |>
  dplyr::mutate(
    row_order = match(
      Model,
      c("M2, df = 4", "M2, df = 3", "M1, df = 4")
    )
  ) |>
  dplyr::arrange(row_order) |>
  dplyr::select(-row_order)

# Main Table 4: five-row publication table with point estimates and intervals.
table_4_development_features <- development_feature_table |>
  dplyr::filter(metric %in% primary_metrics) |>
  dplyr::mutate(
    formatted = format_interval(point_estimate, lower_2.5, upper_97.5, 2),
    inoculum = publication_label(inoculum_id, genotype_ptype)
  ) |>
  dplyr::select(inoculum_id, inoculum, metric, formatted) |>
  tidyr::pivot_wider(names_from = metric, values_from = formatted) |>
  dplyr::mutate(row_order = match(as.character(inoculum_id), development_order)) |>
  dplyr::arrange(row_order) |>
  dplyr::transmute(
    Inoculum = inoculum,
    `Peak log10 RNA load` = peak_logVL,
    `Peak time (dpi)` = peak_dpi,
    `0–5 dpi AUC` = AUC_0_5,
    `Maximum positive slope` = max_positive_slope,
    `Most negative slope` = most_negative_slope
  )

application_feature_labels <- c(
  peak_logVL = "Peak log10 RNA load",
  peak_dpi = "Fitted peak time",
  AUC_0_5 = "0–5-dpi log-scale AUC",
  max_positive_slope = "Maximum positive fitted slope",
  most_negative_slope = "Most negative fitted slope",
  fitted_0dpi_logVL = "Fitted 0-dpi reference level",
  max_net_increase = "Maximum increase above 0-dpi reference",
  centered_AUC_0_5 = "0–5-dpi reference-centered AUC"
)
application_feature_units <- c(
  peak_logVL = "log10 copies/µL",
  peak_dpi = "dpi",
  AUC_0_5 = "log10(copies/µL)·day",
  max_positive_slope = "log10(copies/µL)/day",
  most_negative_slope = "log10(copies/µL)/day",
  fitted_0dpi_logVL = "log10 copies/µL",
  max_net_increase = "log10 units",
  centered_AUC_0_5 = "log10(copies/µL)·day"
)
application_feature_order <- names(application_feature_labels)

# Main Table 5: human-readable GZJY122 feature table.
table_5_GZJY122_features <- application_feature_table |>
  dplyr::filter(metric %in% application_feature_order) |>
  dplyr::mutate(row_order = match(metric, application_feature_order)) |>
  dplyr::arrange(row_order) |>
  dplyr::transmute(
    Feature = unname(application_feature_labels[metric]),
    `Point estimate` = format_fixed(point_estimate, 2),
    `95% conditional interval` = sprintf(
      "%.2f to %.2f", lower_2.5, upper_97.5
    ),
    Unit = unname(application_feature_units[metric])
  )

# Main Table 6: one-row GZJY122 LOBO summary.
table_6_GZJY122_LOBO <- application_result$lobo_summary |>
  dplyr::transmute(
    Group = "Overall",
    `Test observations` = n_predictions,
    Folds = dplyr::n_distinct(
      application_result$lobo$fold[application_result$lobo$status == "OK"]
    ),
    RMSE = format_fixed(RMSE, 3),
    MAE = format_fixed(MAE, 3),
    Bias = format_fixed(bias, 4),
    r = format_fixed(correlation, 3),
    `95% PI coverage` = format_percent(prediction_interval_coverage, 1)
  )

# Supplementary Table S1.
table_S1_development_QC <- development_qc |>
  dplyr::mutate(
    inoculum_order = match(as.character(inoculum_id), development_order),
    batch_number = suppressWarnings(as.numeric(as.character(batch)))
  ) |>
  dplyr::arrange(inoculum_order, batch_number, batch) |>
  dplyr::transmute(
    Inoculum = inoculum_id,
    `Genotype/P-type` = genotype_ptype,
    Batch = batch,
    `Distinct dpi` = n_timepoints,
    Min = min_dpi,
    Max = max_dpi,
    `Early n` = n_early_timepoints,
    `Middle n` = n_middle_timepoints,
    `Late n` = n_late_timepoints,
    `QC status/reason` = exclusion_reason
  )

# Supplementary Table S2.
table_S2_model_hierarchy <- model_hierarchy |>
  dplyr::mutate(
    model_order = match(model_type, c("M0", "M1", "M2")),
    
    convergence_publication = dplyr::case_when(
      singular ~ "Boundary (singular) fit",
      is.na(convergence_message) |
        convergence_message == "" ~ "",
      TRUE ~ convergence_message
    ),
    
    knot_publication = {
      knot_chr <- as.character(internal_knots)
      
      has_knot <- (
        !is.na(knot_chr) &
          nzchar(knot_chr) &
          knot_chr != "None"
      )
      
      knot_num <- rep(NA_real_, length(knot_chr))
      knot_num[has_knot] <- as.numeric(knot_chr[has_knot])
      
      output <- rep("None", length(knot_chr))
      output[has_knot] <- format_fixed(
        knot_num[has_knot],
        1
      )
      
      output
    }
  ) |>
  dplyr::arrange(
    model_order,
    dplyr::desc(spline_df)
  ) |>
  dplyr::transmute(
    Model = model_type,
    df = spline_df,
    `Internal knot(s)` = knot_publication,
    AIC = format_fixed(AIC, 2),
    BIC = format_fixed(BIC, 2),
    LogLik = format_fixed(logLik, 2),
    Sigma = format_fixed(residual_sd, 3),
    Singular = ifelse(singular, "Yes", "No"),
    Convergence = convergence_publication
  )

# Supplementary Table S3: complete M0/M1/M2 hierarchy at both candidate dfs.
complete_lrt_hierarchy <- purrr::map_dfr(cfg$candidate_df, function(k) {
  bundles <- model_bundles[paste0(c("M0", "M1", "M2"), "_df", k)]
  log_likelihood <- vapply(
    bundles,
    function(bundle) as.numeric(stats::logLik(bundle$model)),
    numeric(1)
  )
  parameters <- vapply(
    bundles,
    function(bundle) as.numeric(attr(stats::logLik(bundle$model), "df")),
    numeric(1)
  )
  chi_square <- c(
    NA_real_,
    2 * (log_likelihood[2] - log_likelihood[1]),
    2 * (log_likelihood[3] - log_likelihood[2])
  )
  test_df <- c(
    NA_real_,
    parameters[2] - parameters[1],
    parameters[3] - parameters[2]
  )
  p_value <- c(
    NA_real_,
    stats::pchisq(chi_square[2], df = test_df[2], lower.tail = FALSE),
    stats::pchisq(chi_square[3], df = test_df[3], lower.tail = FALSE)
  )
  tibble::tibble(
    spline_df = k,
    model = c("m0", "m1", "m2"),
    parameters = parameters,
    AIC = vapply(bundles, function(bundle) stats::AIC(bundle$model), numeric(1)),
    BIC = vapply(bundles, function(bundle) stats::BIC(bundle$model), numeric(1)),
    logLik = log_likelihood,
    chi_square = chi_square,
    test_df = test_df,
    p_value = p_value
  )
})

table_S3_likelihood_ratio_tests <- complete_lrt_hierarchy |>
  dplyr::transmute(
    `Spline df` = spline_df,
    Model = model,
    Parameters = parameters,
    AIC = format_fixed(AIC, 2),
    BIC = format_fixed(BIC, 2),
    LogLik = format_fixed(logLik, 2),
    `χ²` = format_fixed(chi_square, 2),
    `Test df` = ifelse(is.na(test_df), "", as.character(as.integer(test_df))),
    p = ifelse(is.na(p_value), "", format(p_value, scientific = TRUE, digits = 3))
  )

# Supplementary Table S4: five primary trajectory features only.
s4_metric_order <- c(
  "AUC_0_5", "most_negative_slope", "max_positive_slope",
  "peak_dpi", "peak_logVL"
)
table_S4_primary_features <- development_feature_table |>
  dplyr::filter(metric %in% s4_metric_order) |>
  dplyr::mutate(
    metric_order = match(metric, s4_metric_order),
    metric_publication = unname(primary_metric_labels[metric])
  ) |>
  dplyr::arrange(metric_order, dplyr::desc(point_estimate)) |>
  dplyr::transmute(
    Inoculum = inoculum_id,
    `Genotype/P-type` = genotype_ptype,
    Metric = metric_publication,
    `Point estimate` = format_fixed(point_estimate, 3),
    `2.5th` = format_fixed(lower_2.5, 3),
    `97.5th` = format_fixed(upper_97.5, 3),
    `Valid draws` = valid_draws
  )

# Supplementary Table S5: reference-centered summaries, including net peak time.
table_S5_centered_features <- development_feature_table |>
  dplyr::filter(metric %in% centered_metrics) |>
  dplyr::mutate(
    metric_order = match(metric, centered_metrics),
    metric_publication = unname(centered_metric_labels[metric])
  ) |>
  dplyr::arrange(metric_order, dplyr::desc(point_estimate)) |>
  dplyr::transmute(
    Inoculum = inoculum_id,
    `Genotype/P-type` = genotype_ptype,
    Metric = metric_publication,
    `Point estimate` = format_fixed(point_estimate, 3),
    `2.5th` = format_fixed(lower_2.5, 3),
    `97.5th` = format_fixed(upper_97.5, 3),
    `Valid draws` = valid_draws
  )

development_genotype_lookup <- stats::setNames(
  as.character(development_files$genotype_ptype),
  as.character(development_files$model_key)
)

# Supplementary Table S6: the five prespecified primary feature contrasts.
table_S6_pairwise_contrasts <- development_pairwise |>
  dplyr::filter(metric %in% primary_metrics) |>
  dplyr::mutate(
    metric_order = match(metric, primary_metrics),
    Metric = unname(primary_metric_labels[metric]),
    Contrast = paste0(
      unname(development_genotype_lookup[as.character(inoculum_A)]),
      " - ",
      unname(development_genotype_lookup[as.character(inoculum_B)])
    )
  ) |>
  dplyr::arrange(metric_order) |>
  dplyr::transmute(
    Metric,
    Contrast,
    `Median difference` = format_fixed(median_difference, 3),
    `2.5th` = format_fixed(lower_2.5, 3),
    `97.5th` = format_fixed(upper_97.5, 3),
    `Fraction of simulated differences > 0` = format_fixed(
      proportion_above_zero, 3
    )
  )

# Supplementary Table S7: manuscript display; the full fold audit remains in
# source_data and includes M1/M2 MAE, singularity, and convergence fields.
table_S7_LOBO_comparison <- development_lobo_fold_comparison |>
  dplyr::mutate(
    batch_number = sub("^.*__batch__", "", as.character(fold)),
    held_out_publication = paste0(genotype_ptype, "_batch", batch_number),
    inoculum_order = match(as.character(inoculum_id), development_order)
  ) |>
  dplyr::arrange(inoculum_order, suppressWarnings(as.numeric(batch_number)), batch_number) |>
  dplyr::transmute(
    `Held-out batch` = held_out_publication,
    Inoculum = inoculum_id,
    n,
    `RMSE M1` = format_fixed(RMSE_M1, 3),
    `RMSE M2` = format_fixed(RMSE_M2, 3),
    `ΔRMSE M2−M1` = format_fixed(delta_RMSE_M2_minus_M1, 3),
    `ΔMAE M2−M1` = format_fixed(delta_MAE_M2_minus_M1, 3),
    `PI coverage M1` = format_percent(PI_coverage_M1, 1),
    `PI coverage M2` = format_percent(PI_coverage_M2, 1)
  )

# Supplementary Table S8.
table_S8_reconstruction <- development_reconstruction_publication |>
  dplyr::mutate(
    row_order = dplyr::if_else(
      as.character(inoculum_id) == "Overall",
      length(development_order) + 1L,
      match(as.character(inoculum_id), development_order)
    )
  ) |>
  dplyr::arrange(row_order) |>
  dplyr::transmute(
    Inoculum = inoculum_id,
    `Genotype/P-type` = dplyr::if_else(
      as.character(inoculum_id) == "Overall", "", as.character(genotype_ptype)
    ),
    n,
    Batches = batches,
    `Residual SD` = dplyr::if_else(
      as.character(inoculum_id) == "Overall",
      "",
      format_fixed(residual_SD, 3)
    ),
    `LOBO RMSE` = format_fixed(LOBO_RMSE, 3),
    `LOBO MAE` = format_fixed(LOBO_MAE, 3),
    `LOBO r` = format_fixed(LOBO_r, 3),
    `LOBO PI coverage` = format_percent(LOBO_PI_coverage, 1)
  )

# Supplementary Table S9: maximum absolute extraction difference across the
# five development inocula at each requested grid size.
table_S9_grid_stability <- development_grid_check$stability |>
  dplyr::group_by(requested_grid_points, actual_grid_points) |>
  dplyr::summarise(
    reference_grid = max(cfg$grid_points),
    peak_level_diff = max(abs_diff_peak_logVL, na.rm = TRUE),
    peak_time_diff = max(abs_diff_peak_dpi, na.rm = TRUE),
    AUC_diff = max(abs_diff_AUC_0_5, na.rm = TRUE),
    max_rate_diff = max(abs_diff_max_positive_slope, na.rm = TRUE),
    decline_rate_diff = max(abs_diff_most_negative_slope, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::transmute(
    `Requested grid` = requested_grid_points,
    `Actual grid` = actual_grid_points,
    Reference = reference_grid,
    `Peak level diff` = format_fixed(peak_level_diff, 6),
    `Peak-time diff` = format_fixed(peak_time_diff, 6),
    `AUC diff` = format_fixed(AUC_diff, 6),
    `Max-rate diff` = format_fixed(max_rate_diff, 6),
    `Decline-rate diff` = format_fixed(decline_rate_diff, 6)
  )

# Supplementary Table S10.
table_S10_GZJY122_QC <- application_result$qc |>
  dplyr::mutate(batch_number = suppressWarnings(as.numeric(as.character(batch)))) |>
  dplyr::arrange(batch_number, batch) |>
  dplyr::transmute(
    Batch = batch,
    `Distinct dpi` = n_timepoints,
    Min = min_dpi,
    Max = max_dpi,
    `Early n` = n_early_timepoints,
    `Middle n` = n_middle_timepoints,
    `Late n` = n_late_timepoints,
    `0-dpi logVL` = format_fixed(observed_0dpi_logVL, 3),
    `Observed peak` = format_fixed(observed_peak_logVL, 3),
    `Peak dpi` = observed_peak_dpi,
    QC = exclusion_reason
  )

s11_metric_labels <- c(
  peak_logVL = "Peak_logVL",
  peak_dpi = "Peak_dpi",
  AUC_0_5 = "AUC",
  max_positive_slope = "Max_Rate",
  most_negative_slope = "Max_Decline_Rate",
  fitted_0dpi_logVL = "Baseline_logVL",
  max_net_increase = "Max_Net_Increase",
  centered_AUC_0_5 = "Baseline_Adjusted_AUC"
)

# Supplementary Table S11.
table_S11_GZJY122_model_features <- dplyr::bind_rows(
  tibble::tibble(
    Category = "Model status",
    Metric = c("Residual SD", "Random-intercept variance", "Singular fit"),
    `Point/status` = c(
      format_fixed(application_result$model_status$residual_sd[1], 3),
      format_fixed(application_result$model_status$random_intercept_variance[1], 4),
      ifelse(application_result$model_status$singular[1], "Yes", "No")
    ),
    `2.5th` = "",
    `97.5th` = "",
    `Valid draws` = ""
  ),
  application_feature_table |>
    dplyr::filter(metric %in% names(s11_metric_labels)) |>
    dplyr::mutate(
      row_order = match(metric, names(s11_metric_labels)),
      Category = dplyr::if_else(
        metric %in% c(
          "fitted_0dpi_logVL", "max_net_increase", "centered_AUC_0_5"
        ),
        "Reference-centered feature",
        "Primary feature"
      )
    ) |>
    dplyr::arrange(row_order) |>
    dplyr::transmute(
      Category,
      Metric = unname(s11_metric_labels[metric]),
      `Point/status` = format_fixed(point_estimate, 3),
      `2.5th` = format_fixed(lower_2.5, 3),
      `97.5th` = format_fixed(upper_97.5, 3),
      `Valid draws` = as.character(valid_draws)
    )
)

# Supplementary Table S12.
table_S12_GZJY122_grid <- application_result$grid_stability$stability |>
  dplyr::transmute(
    `Requested grid` = requested_grid_points,
    `Actual grid` = actual_grid_points,
    `Peak-level diff` = format_fixed(abs_diff_peak_logVL, 6),
    `Peak-time diff` = format_fixed(abs_diff_peak_dpi, 6),
    `AUC diff` = format_fixed(abs_diff_AUC_0_5, 6),
    `Max-rate diff` = format_fixed(abs_diff_max_positive_slope, 6),
    `Decline-rate diff` = format_fixed(abs_diff_most_negative_slope, 6)
  )

# Supplementary Table S13.
table_S13_GZJY122_LOBO <- application_lobo_folds |>
  dplyr::mutate(batch_number = suppressWarnings(as.numeric(as.character(fold)))) |>
  dplyr::arrange(batch_number, fold) |>
  dplyr::transmute(
    `Held-out batch` = fold,
    Status = status,
    `Singular training fit` = ifelse(singular_training_fit, "Yes", "No"),
    n,
    RMSE = format_fixed(RMSE, 3),
    MAE = format_fixed(MAE, 3),
    Bias = format_fixed(bias, 3),
    `95% PI coverage` = format_percent(PI_coverage, 1)
  )

# Supplementary Table S14.
model_metric_labels <- c(
  "Marginal R2" = "Marginal R² (fixed effects)",
  "Conditional R2" = "Conditional R² (fixed + random effects)",
  "Random-intercept variance" = "Random-intercept variance",
  "Residual variance" = "Residual variance",
  "Batch ICC" = "Batch ICC"
)
table_S14_model_metrics <- model_explanatory_metrics |>
  dplyr::transmute(
    Metric = unname(model_metric_labels[metric]),
    Value = format_fixed(value, 3)
  )

# Supplementary Table S15.
table_S15_six_inoculum_features <- six_publication_features |>
  dplyr::mutate(
    row_order = match(as.character(inoculum_id), six_order),
    dataset_role = dplyr::if_else(
      as.character(inoculum_id) == "GZJY122",
      "Locked application",
      "Development"
    )
  ) |>
  dplyr::arrange(row_order) |>
  dplyr::transmute(
    Inoculum = inoculum_id,
    `Genotype/P-type` = genotype_ptype,
    `Dataset role` = dataset_role,
    `Peak logVL` = peak_logVL,
    `Peak dpi` = peak_dpi,
    AUC = AUC_0_5,
    `Max rate` = max_positive_slope,
    `Most negative rate` = most_negative_slope,
    `Centered AUC` = centered_AUC_0_5
  )

s16_metric_labels <- c(
  peak_logVL = "Peak_logVL",
  peak_dpi = "Peak_dpi",
  AUC_0_5 = "AUC",
  max_positive_slope = "Max_Rate",
  most_negative_slope = "Max_Decline_Rate",
  fitted_0dpi_logVL = "Baseline_logVL",
  max_net_increase = "Max_Net_Increase",
  net_peak_dpi = "Net_Peak_dpi",
  centered_AUC_0_5 = "Baseline_Adjusted_AUC"
)

# Supplementary Table S16.
table_S16_GZJY122_GZJY17_contrasts <- six_contrasts |>
  dplyr::filter(metric %in% names(s16_metric_labels)) |>
  dplyr::mutate(row_order = match(metric, names(s16_metric_labels))) |>
  dplyr::arrange(row_order) |>
  dplyr::transmute(
    Metric = unname(s16_metric_labels[metric]),
    `Median difference` = format_fixed(median_difference, 3),
    `2.5th` = format_fixed(lower_2.5, 3),
    `97.5th` = format_fixed(upper_97.5, 3),
    `Fraction of simulated differences > 0` = format_fixed(
      proportion_above_zero, 3
    )
  )

# Supplementary Table S17: GZJY122 sensitivity analyses. The fixed-effect-
# only model retains the locked basis; the four-batch ranges come from the same
# fold definitions used in the primary LOBO assessment.
table_S17_GZJY122_sensitivity <- application_sensitivity_summary |>
  dplyr::transmute(
    Feature = feature,
    `Primary mixed model` = format_fixed(primary_mixed_model, 3),
    `Fixed-effect-only model` = format_fixed(fixed_effect_only_model, 3),
    `Absolute difference` = format_fixed(
      fixed_only_absolute_difference,
      6
    ),
    `Four-batch refit range` = dplyr::if_else(
      successful_refits > 0L,
      sprintf("%.3f to %.3f", four_batch_min, four_batch_max),
      ""
    ),
    `Maximum absolute deviation` = dplyr::if_else(
      successful_refits > 0L,
      format_fixed(maximum_absolute_deviation, 6),
      ""
    ),
    `Successful four-batch refits` = successful_refits,
    Unit = unit
  )

# Figure 4 machine-readable source tables corresponding to the v9.3 panels.
figure4_bubble_source <- development_features |>
  dplyr::transmute(
    Inoculum = inoculum_id,
    `Genotype/P-type` = genotype_ptype,
    `Peak time (dpi)` = peak_dpi,
    `Peak log10 RNA load` = peak_logVL,
    `0–5 dpi AUC` = AUC_0_5
  )

figure4_rate_source <- development_features |>
  dplyr::transmute(
    Inoculum = inoculum_id,
    `Genotype/P-type` = genotype_ptype,
    `Maximum positive slope` = max_positive_slope,
    `Most negative slope` = most_negative_slope
  )

figure4_heatmap_source <- development_features |>
  dplyr::transmute(
    inoculum_id,
    genotype_ptype,
    Peak = peak_logVL,
    AUC = AUC_0_5,
    Positive_slope = max_positive_slope,
    Decline_strength = abs(most_negative_slope)
  ) |>
  tidyr::pivot_longer(
    cols = c(Peak, AUC, Positive_slope, Decline_strength),
    names_to = "Endpoint",
    values_to = "Raw_value"
  ) |>
  dplyr::group_by(Endpoint) |>
  dplyr::mutate(
    Endpoint_mean = mean(Raw_value, na.rm = TRUE),
    Endpoint_SD = stats::sd(Raw_value, na.rm = TRUE),
    Z_score = dplyr::if_else(
      is.finite(Endpoint_SD) & Endpoint_SD > 0,
      (Raw_value - Endpoint_mean) / Endpoint_SD,
      0
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::transmute(
    Inoculum = inoculum_id,
    `Genotype/P-type` = genotype_ptype,
    Endpoint,
    Raw_value,
    Z_score
  )

main_tables <- list(
  Table_1_dataset_roles = table_1_dataset_roles,
  Table_2_RT_qPCR_primers_probes = table_2_primers_probes,
  Table_3_model_specification = table_3_model_specification,
  Table_4_development_features = table_4_development_features,
  Table_5_GZJY122_features = table_5_GZJY122_features,
  Table_6_GZJY122_LOBO = table_6_GZJY122_LOBO
)

supplementary_tables <- list(
  Table_S1_development_structural_QC = table_S1_development_QC,
  Table_S2_complete_model_hierarchy = table_S2_model_hierarchy,
  Table_S3_likelihood_ratio_tests = table_S3_likelihood_ratio_tests,
  Table_S4_development_primary_features = table_S4_primary_features,
  Table_S5_development_centered_features = table_S5_centered_features,
  Table_S6_development_pairwise_contrasts = table_S6_pairwise_contrasts,
  Table_S7_fixed_df_M1_M2_LOBO_by_fold = table_S7_LOBO_comparison,
  Table_S8_development_reconstruction = table_S8_reconstruction,
  Table_S9_development_grid_stability = table_S9_grid_stability,
  Table_S10_GZJY122_structural_QC = table_S10_GZJY122_QC,
  Table_S11_GZJY122_model_and_features = table_S11_GZJY122_model_features,
  Table_S12_GZJY122_grid_stability = table_S12_GZJY122_grid,
  Table_S13_GZJY122_LOBO_folds = table_S13_GZJY122_LOBO,
  Table_S14_model_explanatory_metrics = table_S14_model_metrics,
  Table_S15_six_inoculum_features = table_S15_six_inoculum_features,
  Table_S16_GZJY122_minus_GZJY17_contrasts =
    table_S16_GZJY122_GZJY17_contrasts,
  Table_S17_GZJY122_sensitivity = table_S17_GZJY122_sensitivity
)

# Fail early if a later code change breaks the manuscript table contract.
validate_table_contract <- function(tables, expected_columns, expected_rows) {
  missing_tables <- setdiff(names(expected_columns), names(tables))
  if (length(missing_tables) > 0L) {
    stop(
      "Missing manuscript table object(s): ",
      paste(missing_tables, collapse = ", "),
      call. = FALSE
    )
  }

  purrr::iwalk(expected_columns, function(columns, table_name) {
    actual_columns <- names(tables[[table_name]])
    if (!identical(actual_columns, columns)) {
      stop(
        "Column mismatch in ", table_name, ".\nExpected: ",
        paste(columns, collapse = " | "), "\nActual: ",
        paste(actual_columns, collapse = " | "),
        call. = FALSE
      )
    }
    if (!is.null(expected_rows[[table_name]]) &&
        nrow(tables[[table_name]]) != expected_rows[[table_name]]) {
      stop(
        "Row-count mismatch in ", table_name, ": expected ",
        expected_rows[[table_name]], ", found ",
        nrow(tables[[table_name]]),
        call. = FALSE
      )
    }
  })
  invisible(TRUE)
}

expected_main_columns <- list(
  Table_1_dataset_roles = c(
    "Analysis role", "Inoculum(s)", "Genotype", "No. of inocula",
    "Valid batches", "Observations", "Purpose"
  ),
  Table_2_RT_qPCR_primers_probes = c(
    "Norovirus-GII", "Name", "Sequence (5′–3′)", "Reference"
  ),
  Table_3_model_specification = c(
    "Model", "Internal knot", "AIC", "BIC", "Residual SD", "Role"
  ),
  Table_4_development_features = c(
    "Inoculum", "Peak log10 RNA load", "Peak time (dpi)", "0–5 dpi AUC",
    "Maximum positive slope", "Most negative slope"
  ),
  Table_5_GZJY122_features = c(
    "Feature", "Point estimate", "95% conditional interval", "Unit"
  ),
  Table_6_GZJY122_LOBO = c(
    "Group", "Test observations", "Folds", "RMSE", "MAE", "Bias", "r",
    "95% PI coverage"
  )
)

expected_supplementary_columns <- list(
  Table_S1_development_structural_QC = c(
    "Inoculum", "Genotype/P-type", "Batch", "Distinct dpi", "Min", "Max",
    "Early n", "Middle n", "Late n", "QC status/reason"
  ),
  Table_S2_complete_model_hierarchy = c(
    "Model", "df", "Internal knot(s)", "AIC", "BIC", "LogLik", "Sigma",
    "Singular", "Convergence"
  ),
  Table_S3_likelihood_ratio_tests = c(
    "Spline df", "Model", "Parameters", "AIC", "BIC", "LogLik", "χ²",
    "Test df", "p"
  ),
  Table_S4_development_primary_features = c(
    "Inoculum", "Genotype/P-type", "Metric", "Point estimate", "2.5th",
    "97.5th", "Valid draws"
  ),
  Table_S5_development_centered_features = c(
    "Inoculum", "Genotype/P-type", "Metric", "Point estimate", "2.5th",
    "97.5th", "Valid draws"
  ),
  Table_S6_development_pairwise_contrasts = c(
    "Metric", "Contrast", "Median difference", "2.5th", "97.5th",
    "Fraction of simulated differences > 0"
  ),
  Table_S7_fixed_df_M1_M2_LOBO_by_fold = c(
    "Held-out batch", "Inoculum", "n", "RMSE M1", "RMSE M2",
    "ΔRMSE M2−M1", "ΔMAE M2−M1", "PI coverage M1", "PI coverage M2"
  ),
  Table_S8_development_reconstruction = c(
    "Inoculum", "Genotype/P-type", "n", "Batches", "Residual SD",
    "LOBO RMSE", "LOBO MAE", "LOBO r", "LOBO PI coverage"
  ),
  Table_S9_development_grid_stability = c(
    "Requested grid", "Actual grid", "Reference", "Peak level diff",
    "Peak-time diff", "AUC diff", "Max-rate diff", "Decline-rate diff"
  ),
  Table_S10_GZJY122_structural_QC = c(
    "Batch", "Distinct dpi", "Min", "Max", "Early n", "Middle n", "Late n",
    "0-dpi logVL", "Observed peak", "Peak dpi", "QC"
  ),
  Table_S11_GZJY122_model_and_features = c(
    "Category", "Metric", "Point/status", "2.5th", "97.5th", "Valid draws"
  ),
  Table_S12_GZJY122_grid_stability = c(
    "Requested grid", "Actual grid", "Peak-level diff", "Peak-time diff",
    "AUC diff", "Max-rate diff", "Decline-rate diff"
  ),
  Table_S13_GZJY122_LOBO_folds = c(
    "Held-out batch", "Status", "Singular training fit", "n", "RMSE", "MAE",
    "Bias", "95% PI coverage"
  ),
  Table_S14_model_explanatory_metrics = c("Metric", "Value"),
  Table_S15_six_inoculum_features = c(
    "Inoculum", "Genotype/P-type", "Dataset role", "Peak logVL", "Peak dpi",
    "AUC", "Max rate", "Most negative rate", "Centered AUC"
  ),
  Table_S16_GZJY122_minus_GZJY17_contrasts = c(
    "Metric", "Median difference", "2.5th", "97.5th",
    "Fraction of simulated differences > 0"
  ),
  Table_S17_GZJY122_sensitivity = c(
    "Feature", "Primary mixed model", "Fixed-effect-only model",
    "Absolute difference", "Four-batch refit range",
    "Maximum absolute deviation", "Successful four-batch refits", "Unit"
  )
)

expected_main_rows <- list(
  Table_1_dataset_roles = 6L,
  Table_2_RT_qPCR_primers_probes = 3L,
  Table_3_model_specification = 3L,
  Table_4_development_features = 5L,
  Table_5_GZJY122_features = 8L,
  Table_6_GZJY122_LOBO = 1L
)

expected_supplementary_rows <- list(
  Table_S1_development_structural_QC = nrow(development_qc),
  Table_S2_complete_model_hierarchy = 6L,
  Table_S3_likelihood_ratio_tests = 6L,
  Table_S4_development_primary_features = 25L,
  Table_S5_development_centered_features = 20L,
  Table_S6_development_pairwise_contrasts = 50L,
  Table_S7_fixed_df_M1_M2_LOBO_by_fold =
    dplyr::n_distinct(development_data$inoculum_batch),
  Table_S8_development_reconstruction = 6L,
  Table_S9_development_grid_stability = length(cfg$grid_points),
  Table_S10_GZJY122_structural_QC =
    dplyr::n_distinct(application_result$data$batch),
  Table_S11_GZJY122_model_and_features = 11L,
  Table_S12_GZJY122_grid_stability = length(cfg$grid_points),
  Table_S13_GZJY122_LOBO_folds =
    dplyr::n_distinct(application_result$lobo$fold),
  Table_S14_model_explanatory_metrics = 5L,
  Table_S15_six_inoculum_features = 6L,
  Table_S16_GZJY122_minus_GZJY17_contrasts = 9L,
  Table_S17_GZJY122_sensitivity = length(application_sensitivity_metrics)
)

validate_table_contract(main_tables, expected_main_columns, expected_main_rows)
validate_table_contract(
  supplementary_tables,
  expected_supplementary_columns,
  expected_supplementary_rows
)

source_tables <- list(
  Figure_4A_peak_time_AUC_bubble_source = figure4_bubble_source,
  Figure_4B_change_rate_source = figure4_rate_source,
  Figure_4C_standardized_endpoint_source = figure4_heatmap_source,
  development_all_feature_intervals = development_feature_table,
  development_LOBO_M1_M2_fold_audit = development_lobo_fold_comparison,
  development_raw_data = development_raw,
  development_model_data = development_data,
  development_legacy_outcome_flags = legacy_outcome_flags,
  development_fixed_prediction_curves = development_prediction,
  development_feature_simulation_draws = development_draws,
  development_derivative_curves = development_derivatives,
  development_LOBO_M1_predictions = development_lobo_m1,
  development_LOBO_M2_predictions = development_lobo_m2,
  development_residual_diagnostics = development_diagnostics,
  development_random_intercept_BLUPs = development_random_intercepts,
  application_raw_data = application_raw,
  application_model_data = application_result$data,
  application_fixed_prediction_curve = application_result$prediction,
  application_feature_simulation_draws = application_result$draws,
  application_residual_diagnostics = application_result$diagnostics,
  application_derivative_curve = application_result$derivatives,
  application_fixed_effect_only_model_status = application_fixed_only_status,
  application_fixed_effect_only_prediction_curve =
    application_fixed_only_prediction,
  application_fixed_effect_only_features = application_fixed_only_features,
  application_LOBO_training_feature_estimates =
    application_lobo_training_feature_estimates,
  application_sensitivity_feature_summary = application_sensitivity_summary,
  six_inoculum_model_data = six_data,
  six_inoculum_fixed_prediction_curves = six_prediction,
  six_inoculum_feature_simulation_draws = six_draws,
  GZJY122_minus_GZJY17_difference_curve = six_difference_curve
)

write_named_tables(main_tables, dirs$main_tables)
write_named_tables(supplementary_tables, dirs$supplementary_tables)
write_named_tables(source_tables, dirs$source_data)

saveRDS(primary_bundle, file.path(dirs$model_objects, "development_primary_M2_df4_bundle.rds"))
saveRDS(application_result$bundle, file.path(dirs$model_objects, "GZJY122_locked_single_bundle.rds"))
saveRDS(
  application_fixed_only_model,
  file.path(dirs$model_objects, "GZJY122_fixed_effect_only_sensitivity_model.rds")
)
saveRDS(
  application_lobo_feature_fit_objects,
  file.path(dirs$model_objects, "GZJY122_LOBO_training_feature_fits.rds")
)
saveRDS(six_bundle, file.path(dirs$model_objects, "six_inoculum_fixed_M2_bundle.rds"))

write_model_summary(primary_bundle$model, file.path(dirs$logs, "development_primary_model_summary.txt"))
write_model_summary(application_result$bundle$model, file.path(dirs$logs, "GZJY122_locked_model_summary.txt"))
write_model_summary(
  application_fixed_only_model,
  file.path(dirs$logs, "GZJY122_fixed_effect_only_sensitivity_summary.txt")
)
write_model_summary(six_bundle$model, file.path(dirs$logs, "six_inoculum_model_summary.txt"))
write_session_info(file.path(dirs$logs, "sessionInfo.txt"))
write_csv(specification_check, file.path(dirs$logs, "development_vs_locked_specification.csv"))
publication_table_notes <- tibble::tribble(
  ~table, ~note,
  "Table S6",
  paste(
    "The fraction of simulated differences above zero is descriptive;",
    "it is not a posterior probability or p-value. No multiplicity adjustment was applied."
  ),
  "Table S7",
  paste(
    "The manuscript display table omits full MAE and singularity fields;",
    "these are retained in source_data/development_LOBO_M1_M2_fold_audit.csv."
  ),
  "Table S8",
  "The overall row reports aggregate M2 LOBO metrics; inoculum rows report inoculum-specific reconstruction summaries.",
  "Table S14",
  "R2 values are Nakagawa-type mixed-model summaries; ICC is based on the random-intercept and residual variances.",
  "Table S16",
  paste(
    "The fraction of simulated differences above zero is descriptive;",
    "it is not a posterior probability or p-value. No multiplicity adjustment was applied."
  ),
  "Table S17",
  paste(
    "Sensitivity analyses retained the locked spline basis and did not trigger",
    "model or spline reselection. The fixed-effect-only model is descriptive,",
    "and the four-batch ranges summarize feature estimates from the five LOBO",
    "training fits rather than held-out prediction performance."
  )
)
write_csv(
  publication_table_notes,
  file.path(dirs$logs, "publication_table_notes.csv")
)
write_run_metadata(
  cfg,
  development_data,
  application_result$data,
  selected_basis,
  file.path(dirs$logs, "analysis_metadata.csv"),
  version = hunov_analysis_version
)
write_output_manifest(output_root)

message(
  "HuNoV-LRPF full publication analysis completed: ",
  normalizePath(output_root, winslash = "/", mustWork = TRUE)
)
