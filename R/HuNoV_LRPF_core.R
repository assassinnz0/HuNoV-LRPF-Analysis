# HuNoV-LRPF core functions
# Shared by the manuscript analysis and the standalone locked-application script.

require_hunov_packages <- function() {
  packages <- c(
    "dplyr", "tidyr", "purrr", "tibble", "readr", "ggplot2",
    "splines", "lme4", "MASS", "DescTools", "Matrix", "performance"
  )
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing required packages: ", paste(missing, collapse = ", "),
      ". Install them before running the analysis.",
      call. = FALSE
    )
  }
}

resolve_data_dir <- function(required_files, project_dir, configured_dir = NULL, interactive_picker = TRUE) {
  candidates <- unique(c(
    configured_dir,
    Sys.getenv("HUNOV_DATA_DIR", unset = ""),
    file.path(project_dir, "data"),
    file.path(project_dir, "input"),
    project_dir,
    getwd()
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates) & dir.exists(candidates)]
  is_complete <- function(path) all(file.exists(file.path(path, required_files)))
  complete <- candidates[vapply(candidates, is_complete, logical(1))]
  if (length(complete) > 0L) {
    return(normalizePath(complete[1], winslash = "/", mustWork = TRUE))
  }

  if (interactive_picker && interactive()) {
    selected <- if (.Platform$OS.type == "windows") {
      utils::choose.dir(default = project_dir, caption = "Select the folder containing the HuNoV CSV files")
    } else {
      chosen_file <- tryCatch(file.choose(), error = function(e) NA_character_)
      if (is.na(chosen_file)) NA_character_ else dirname(chosen_file)
    }
    if (!is.na(selected) && nzchar(selected) && dir.exists(selected) && is_complete(selected)) {
      return(normalizePath(selected, winslash = "/", mustWork = TRUE))
    }
  }

  stop(
    "Could not locate all required input files: ", paste(required_files, collapse = ", "),
    ". Put them in the script directory, a data/ or input/ subdirectory, ",
    "set data_dir in the configuration, or define HUNOV_DATA_DIR.",
    call. = FALSE
  )
}

make_output_dirs <- function(root) {
  dirs <- list(
    root = root,
    tables = file.path(root, "tables"),
    plots = file.path(root, "plots"),
    objects = file.path(root, "objects")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}

write_csv <- function(x, path) {
  readr::write_csv(x, path, na = "")
  invisible(path)
}

save_tiff <- function(plot, path, width = 8, height = 6, dpi = 600) {
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
  invisible(path)
}

theme_hunov <- function(base_size = 11, base_family = "Arial") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank()
    )
}

read_inoculum_csv <- function(file, inoculum_id, genotype_ptype, model_key = inoculum_id) {
  dat <- readr::read_csv(file, show_col_types = FALSE, progress = FALSE)
  required <- c("batch", "dpi", "VL")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop(basename(file), " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  batch_chr <- trimws(as.character(dat$batch))
  if (any(is.na(batch_chr) | !nzchar(batch_chr))) {
    stop(basename(file), " contains missing or empty batch identifiers.", call. = FALSE)
  }
  dat <- dat |>
    dplyr::transmute(
      inoculum_id = as.character(inoculum_id),
      genotype_ptype = as.character(genotype_ptype),
      model_key = as.character(model_key),
      source_file = basename(file),
      batch_raw = batch_chr,
      batch = factor(batch_chr),
      dpi = as.numeric(dpi),
      VL = as.numeric(VL)
    )
  if (any(!is.finite(dat$dpi))) stop(basename(file), " contains non-numeric dpi values.", call. = FALSE)
  if (any(!is.finite(dat$VL) | dat$VL <= 0)) {
    stop(basename(file), " contains missing, zero, negative, or non-finite VL values.", call. = FALSE)
  }
  dat |>
    dplyr::mutate(logVL = log10(VL)) |>
    dplyr::arrange(batch, dpi)
}

combine_inocula <- function(input_table, data_dir) {
  purrr::pmap_dfr(
    input_table,
    function(inoculum_id, genotype_ptype, file, model_key = inoculum_id) {
      read_inoculum_csv(
        file = file.path(data_dir, file),
        inoculum_id = inoculum_id,
        genotype_ptype = genotype_ptype,
        model_key = model_key
      )
    }
  ) |>
    dplyr::mutate(
      model_key = factor(model_key, levels = unique(model_key)),
      inoculum_id = factor(inoculum_id, levels = unique(inoculum_id)),
      batch = factor(batch),
      inoculum_batch = interaction(model_key, batch, drop = TRUE, sep = "__batch__")
    ) |>
    dplyr::arrange(model_key, batch, dpi)
}

structural_qc <- function(data, cfg) {
  q <- cfg$qc
  data |>
    dplyr::group_by(model_key, inoculum_id, genotype_ptype, batch) |>
    dplyr::summarise(
      n_obs = dplyr::n(),
      n_timepoints = dplyr::n_distinct(dpi),
      n_duplicate_rows = n_obs - n_timepoints,
      min_dpi = min(dpi),
      max_dpi = max(dpi),
      all_within_boundary = all(dpi >= cfg$trajectory_window[1] & dpi <= cfg$trajectory_window[2]),
      has_dpi0 = any(dpi == 0),
      n_early_timepoints = dplyr::n_distinct(dpi[dpi >= q$early_window[1] & dpi <= q$early_window[2]]),
      n_middle_timepoints = dplyr::n_distinct(dpi[dpi >= q$middle_window[1] & dpi <= q$middle_window[2]]),
      n_late_timepoints = dplyr::n_distinct(dpi[dpi >= q$late_window[1] & dpi <= q$late_window[2]]),
      observed_0dpi_logVL = if (any(dpi == 0)) mean(logVL[dpi == 0]) else NA_real_,
      observed_peak_logVL = max(logVL),
      observed_peak_dpi = dpi[which.max(logVL)][1],
      observed_peak_increase = observed_peak_logVL - observed_0dpi_logVL,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      eligible =
        n_timepoints >= q$min_timepoints &
        (!q$require_dpi0 | has_dpi0) &
        max_dpi >= q$min_followup_dpi &
        (!q$require_all_within_boundary | all_within_boundary) &
        (q$allow_duplicate_dpi | n_duplicate_rows == 0L) &
        n_early_timepoints >= q$min_early_timepoints &
        n_middle_timepoints >= q$min_middle_timepoints &
        n_late_timepoints >= q$min_late_timepoints,
      exclusion_reason = dplyr::case_when(
        eligible ~ "Included",
        n_timepoints < q$min_timepoints ~ "Too few distinct time points",
        q$require_dpi0 & !has_dpi0 ~ "Missing 0 dpi",
        max_dpi < q$min_followup_dpi ~ "Insufficient late follow-up",
        q$require_all_within_boundary & !all_within_boundary ~ "Observation outside analysis window",
        !q$allow_duplicate_dpi & n_duplicate_rows > 0L ~ "Duplicate batch-dpi rows",
        n_early_timepoints < q$min_early_timepoints ~ "Insufficient early-window coverage",
        n_middle_timepoints < q$min_middle_timepoints ~ "Insufficient middle-window coverage",
        n_late_timepoints < q$min_late_timepoints ~ "Insufficient late-window coverage",
        TRUE ~ "Other structural QC failure"
      )
    )
}

apply_structural_qc <- function(data, qc_table, cfg) {
  eligible <- qc_table |>
    dplyr::filter(eligible) |>
    dplyr::transmute(model_key = as.character(model_key), batch = as.character(batch))
  filtered <- data |>
    dplyr::mutate(model_key_chr = as.character(model_key), batch_chr = as.character(batch)) |>
    dplyr::inner_join(eligible, by = c("model_key_chr" = "model_key", "batch_chr" = "batch")) |>
    dplyr::select(-model_key_chr, -batch_chr) |>
    droplevels()

  batch_counts <- filtered |>
    dplyr::group_by(model_key, inoculum_id, genotype_ptype) |>
    dplyr::summarise(n_batches = dplyr::n_distinct(batch), .groups = "drop")
  if (any(batch_counts$n_batches < cfg$qc$min_batches)) {
    stop("At least one inoculum has fewer than ", cfg$qc$min_batches, " eligible batches.", call. = FALSE)
  }
  list(data = filtered, counts = batch_counts)
}

make_basis_info_df <- function(dpi, spline_df, cfg) {
  mat <- splines::bs(
    dpi,
    df = spline_df,
    degree = cfg$bs_degree,
    Boundary.knots = cfg$trajectory_window,
    intercept = cfg$bs_intercept
  )
  list(
    spline_df = as.integer(spline_df),
    degree = cfg$bs_degree,
    knots = as.numeric(attr(mat, "knots")),
    boundary_knots = as.numeric(attr(mat, "Boundary.knots")),
    intercept = cfg$bs_intercept,
    columns = paste0("B", seq_len(ncol(mat)))
  )
}

make_locked_basis_info <- function(df, internal_knots, boundary_knots, degree = 3L, intercept = FALSE) {
  probe <- splines::bs(
    seq(boundary_knots[1], boundary_knots[2], length.out = 20),
    knots = internal_knots,
    degree = degree,
    Boundary.knots = boundary_knots,
    intercept = intercept
  )
  if (ncol(probe) != as.integer(df)) {
    stop(
      "Declared spline df does not match the number of generated basis columns.",
      call. = FALSE
    )
  }
  list(
    spline_df = as.integer(df),
    degree = as.integer(degree),
    knots = as.numeric(internal_knots),
    boundary_knots = as.numeric(boundary_knots),
    intercept = isTRUE(intercept),
    columns = paste0("B", seq_len(ncol(probe)))
  )
}

basis_matrix <- function(dpi, basis_info) {
  mat <- splines::bs(
    dpi,
    knots = basis_info$knots,
    degree = basis_info$degree,
    Boundary.knots = basis_info$boundary_knots,
    intercept = basis_info$intercept
  )
  out <- as.data.frame(mat)
  names(out) <- basis_info$columns
  out
}

add_basis <- function(data, basis_info) {
  old <- grep("^B[0-9]+$", names(data), value = TRUE)
  clean <- dplyr::select(data, -dplyr::any_of(old))
  dplyr::bind_cols(clean, basis_matrix(clean$dpi, basis_info))
}

make_formula <- function(model_type, basis_info, group_col) {
  basis <- paste(basis_info$columns, collapse = " + ")
  fixed <- switch(
    model_type,
    single = basis,
    M0 = basis,
    M1 = paste("model_key +", basis),
    M2 = paste0("model_key * (", basis, ")"),
    stop("Unknown model type: ", model_type, call. = FALSE)
  )
  stats::as.formula(paste0("logVL ~ ", fixed, " + (1 | ", group_col, ")"))
}

convergence_message <- function(model) {
  msg <- model@optinfo$conv$lme4$messages
  if (is.null(msg)) "" else paste(msg, collapse = " | ")
}

fit_ml_model <- function(formula, data) {
  optimizers <- c("bobyqa", "Nelder_Mead")
  last_fit <- NULL
  last_error <- NULL
  for (opt in optimizers) {
    fit <- tryCatch(
      lme4::lmer(
        formula,
        data = data,
        REML = FALSE,
        control = lme4::lmerControl(
          optimizer = opt,
          optCtrl = list(maxfun = 1e7),
          calc.derivs = TRUE
        )
      ),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (!is.null(fit)) {
      last_fit <- fit
      if (identical(convergence_message(fit), "")) return(fit)
    }
  }
  if (!is.null(last_fit)) return(last_fit)
  stop("Model fitting failed: ", last_error$message, call. = FALSE)
}

fit_bundle <- function(data, model_type, basis_info, group_col) {
  data_basis <- add_basis(data, basis_info)
  formula <- make_formula(model_type, basis_info, group_col)
  model <- fit_ml_model(formula, data_basis)
  list(
    model = model,
    data = data_basis,
    basis = basis_info,
    model_type = model_type,
    group_col = group_col,
    formula = formula
  )
}

model_row <- function(bundle) {
  tibble::tibble(
    model_type = bundle$model_type,
    spline_df = bundle$basis$spline_df,
    internal_knots = if (length(bundle$basis$knots) == 0L) "None" else paste(bundle$basis$knots, collapse = ";"),
    formula = paste(deparse(stats::formula(bundle$model)), collapse = " "),
    n = stats::nobs(bundle$model),
    logLik = as.numeric(stats::logLik(bundle$model)),
    AIC = stats::AIC(bundle$model),
    BIC = stats::BIC(bundle$model),
    residual_sd = stats::sigma(bundle$model),
    singular = lme4::isSingular(bundle$model, tol = 1e-4),
    convergence_message = convergence_message(bundle$model)
  )
}

select_df_by_bic <- function(model_table, model_type = "M2", tolerance = 2) {
  tab <- dplyr::filter(model_table, .data$model_type == model_type)
  best <- min(tab$BIC)
  tab |>
    dplyr::filter(BIC <= best + tolerance) |>
    dplyr::arrange(spline_df) |>
    dplyr::slice(1) |>
    dplyr::pull(spline_df)
}

make_dpi_grid <- function(n_points, window, required_points = NULL) {
  required_points <- as.numeric(c(window, required_points))
  required_points <- required_points[
    is.finite(required_points) &
      required_points >= window[1] &
      required_points <= window[2]
  ]
  sort(unique(c(
    seq(window[1], window[2], length.out = n_points),
    required_points
  )))
}

new_fixed_grid <- function(data, dpi_grid, basis_info) {
  levels_key <- levels(data$model_key)
  grid <- tidyr::expand_grid(model_key = levels_key, dpi = dpi_grid) |>
    dplyr::mutate(model_key = factor(model_key, levels = levels_key))
  lookup <- data |>
    dplyr::distinct(model_key, inoculum_id, genotype_ptype)
  grid |>
    dplyr::left_join(lookup, by = "model_key") |>
    dplyr::mutate(
      batch = factor(levels(data$batch)[1], levels = levels(data$batch)),
      inoculum_batch = factor(levels(data$inoculum_batch)[1], levels = levels(data$inoculum_batch))
    ) |>
    add_basis(basis_info)
}

fixed_model_matrix <- function(model, newdata) {
  fixed_formula <- lme4::nobars(stats::formula(model))
  fixed_terms <- stats::delete.response(stats::terms(fixed_formula))
  stats::model.matrix(fixed_terms, newdata)
}

positive_definite_vcov <- function(model) {
  vc <- as.matrix(stats::vcov(model))
  eig <- tryCatch(eigen(vc, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
  if (any(!is.finite(eig)) || any(eig <= 1e-10)) vc <- as.matrix(Matrix::nearPD(vc)$mat)
  vc
}

predict_fixed <- function(bundle, newdata) {
  X <- fixed_model_matrix(bundle$model, newdata)
  beta <- lme4::fixef(bundle$model)
  vc <- positive_definite_vcov(bundle$model)
  fit <- as.numeric(X %*% beta)
  se <- sqrt(pmax(diag(X %*% vc %*% t(X)), 0))
  newdata |>
    dplyr::mutate(predicted = fit, fixed_se = se, lower = fit - 1.96 * se, upper = fit + 1.96 * se)
}

safe_auc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2L) return(NA_real_)
  DescTools::AUC(x[ok], y[ok], method = "trapezoid")
}

extract_features <- function(curve, auc_window) {
  curve <- curve |>
    dplyr::filter(is.finite(dpi), is.finite(predicted)) |>
    dplyr::arrange(dpi)
  if (nrow(curve) < 3L) return(tibble::tibble())
  peak_i <- which.max(curve$predicted)
  rates <- diff(curve$predicted) / diff(curve$dpi)
  mid <- (curve$dpi[-1] + curve$dpi[-nrow(curve)]) / 2
  auc_data <- dplyr::filter(curve, dpi >= auc_window[1], dpi <= auc_window[2])
  max_i <- which.max(rates)
  min_i <- which.min(rates)
  max_rate <- rates[max_i]
  min_rate <- rates[min_i]
  tibble::tibble(
    peak_logVL = curve$predicted[peak_i],
    peak_dpi = curve$dpi[peak_i],
    AUC_0_5 = safe_auc(auc_data$dpi, auc_data$predicted),
    max_positive_slope = max(max_rate, 0),
    max_positive_slope_dpi = if (max_rate > 0) mid[max_i] else NA_real_,
    most_negative_slope = min(min_rate, 0),
    most_negative_slope_dpi = if (min_rate < 0) mid[min_i] else NA_real_
  )
}

extract_centered_features <- function(curve, reference_dpi, auc_window) {
  curve <- curve |>
    dplyr::filter(is.finite(dpi), is.finite(predicted)) |>
    dplyr::arrange(dpi)
  ref_i <- which.min(abs(curve$dpi - reference_dpi))
  baseline <- curve$predicted[ref_i]
  delta <- curve$predicted - baseline
  peak_i <- which.max(delta)
  auc_i <- curve$dpi >= auc_window[1] & curve$dpi <= auc_window[2]
  tibble::tibble(
    fitted_0dpi_logVL = baseline,
    max_net_increase = delta[peak_i],
    net_peak_dpi = curve$dpi[peak_i],
    centered_AUC_0_5 = safe_auc(curve$dpi[auc_i], delta[auc_i])
  )
}

feature_table <- function(prediction, cfg) {
  primary <- prediction |>
    dplyr::group_by(model_key, inoculum_id, genotype_ptype) |>
    dplyr::group_modify(~extract_features(.x, cfg$auc_window)) |>
    dplyr::ungroup()
  centered <- prediction |>
    dplyr::group_by(model_key, inoculum_id, genotype_ptype) |>
    dplyr::group_modify(~extract_centered_features(.x, cfg$reference_dpi, cfg$auc_window)) |>
    dplyr::ungroup()
  dplyr::left_join(primary, centered, by = c("model_key", "inoculum_id", "genotype_ptype"))
}

simulate_feature_draws <- function(bundle, prediction_grid, cfg, seed = cfg$seed) {
  set.seed(seed)
  X <- fixed_model_matrix(bundle$model, prediction_grid)
  beta_draws <- MASS::mvrnorm(cfg$n_sim, mu = lme4::fixef(bundle$model), Sigma = positive_definite_vcov(bundle$model))
  values <- X %*% t(beta_draws)
  split_index <- split(seq_len(nrow(prediction_grid)), prediction_grid$model_key)

  purrr::imap_dfr(split_index, function(rows, key) {
    meta <- prediction_grid[rows[1], c("model_key", "inoculum_id", "genotype_ptype")]
    purrr::map_dfr(seq_len(cfg$n_sim), function(i) {
      curve <- tibble::tibble(dpi = prediction_grid$dpi[rows], predicted = values[rows, i])
      dplyr::bind_cols(
        tibble::tibble(simulation = i),
        meta,
        extract_features(curve, cfg$auc_window),
        extract_centered_features(curve, cfg$reference_dpi, cfg$auc_window)
      )
    })
  })
}

summarise_feature_draws <- function(draws) {
  id_cols <- c("model_key", "inoculum_id", "genotype_ptype", "simulation")
  draws |>
    tidyr::pivot_longer(-dplyr::all_of(id_cols), names_to = "metric", values_to = "value") |>
    dplyr::group_by(model_key, inoculum_id, genotype_ptype, metric) |>
    dplyr::summarise(
      n_valid = sum(is.finite(value)),
      median = stats::median(value, na.rm = TRUE),
      lower_2.5 = stats::quantile(value, 0.025, na.rm = TRUE, names = FALSE),
      upper_97.5 = stats::quantile(value, 0.975, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
}

random_intercept_variance <- function(model, group_col) {
  vc <- as.data.frame(lme4::VarCorr(model))
  out <- vc |>
    dplyr::filter(grp == group_col, var1 == "(Intercept)") |>
    dplyr::pull(vcov)
  if (length(out) == 0L || !is.finite(out[1])) 0 else as.numeric(out[1])
}

predict_new_batch <- function(bundle, newdata) {
  new_basis <- add_basis(newdata, bundle$basis)
  X <- fixed_model_matrix(bundle$model, new_basis)
  beta <- lme4::fixef(bundle$model)
  vc <- positive_definite_vcov(bundle$model)
  pred <- as.numeric(X %*% beta)
  fixed_se <- sqrt(pmax(diag(X %*% vc %*% t(X)), 0))
  pred_se <- sqrt(fixed_se^2 + random_intercept_variance(bundle$model, bundle$group_col) + stats::sigma(bundle$model)^2)
  newdata |>
    dplyr::mutate(
      predicted = pred,
      residual = logVL - pred,
      squared_error = residual^2,
      absolute_error = abs(residual),
      prediction_lower = pred - 1.96 * pred_se,
      prediction_upper = pred + 1.96 * pred_se,
      covered = logVL >= prediction_lower & logVL <= prediction_upper
    )
}

run_lobo <- function(data, model_type, basis_info, group_col, fold_col) {
  folds <- unique(as.character(data[[fold_col]]))
  purrr::map_dfr(folds, function(fold) {
    test <- data[as.character(data[[fold_col]]) == fold, , drop = FALSE]
    train <- droplevels(data[as.character(data[[fold_col]]) != fold, , drop = FALSE])
    bundle <- tryCatch(fit_bundle(train, model_type, basis_info, group_col), error = function(e) NULL)
    if (is.null(bundle)) {
      return(test |>
        dplyr::mutate(fold = fold, status = "fit_failed", predicted = NA_real_, residual = NA_real_,
                      squared_error = NA_real_, absolute_error = NA_real_, prediction_lower = NA_real_,
                      prediction_upper = NA_real_, covered = NA))
    }
    predict_new_batch(bundle, test) |>
      dplyr::mutate(
        fold = fold,
        status = "OK",
        singular = lme4::isSingular(bundle$model, tol = 1e-4),
        convergence_message = convergence_message(bundle$model)
      )
  })
}

summarise_lobo <- function(predictions) {
  valid <- dplyr::filter(predictions, status == "OK", is.finite(predicted))
  tibble::tibble(
    n_predictions = nrow(valid),
    RMSE = sqrt(mean(valid$squared_error)),
    MAE = mean(valid$absolute_error),
    bias = mean(valid$residual),
    correlation = stats::cor(valid$logVL, valid$predicted),
    prediction_interval_coverage = mean(valid$covered)
  )
}

residual_diagnostics <- function(bundle) {
  bundle$data |>
    dplyr::mutate(
      fitted = stats::fitted(bundle$model),
      residual = stats::residuals(bundle$model),
      standardized_residual = residual / stats::sigma(bundle$model)
    )
}

model_status <- function(bundle) {
  vc <- as.data.frame(lme4::VarCorr(bundle$model))
  tibble::tibble(
    n = stats::nobs(bundle$model),
    residual_sd = stats::sigma(bundle$model),
    random_intercept_variance = random_intercept_variance(bundle$model, bundle$group_col),
    singular = lme4::isSingular(bundle$model, tol = 1e-4),
    convergence_message = convergence_message(bundle$model),
    formula = paste(deparse(stats::formula(bundle$model)), collapse = " ")
  )
}

grid_stability <- function(bundle, data, cfg) {
  results <- purrr::map_dfr(cfg$grid_points, function(n) {
    dpi <- make_dpi_grid(n, cfg$trajectory_window, c(cfg$auc_window, cfg$reference_dpi))
    grid <- new_fixed_grid(data, dpi, bundle$basis)
    pred <- predict_fixed(bundle, grid)
    feature_table(pred, cfg) |>
      dplyr::mutate(requested_grid_points = n, actual_grid_points = length(dpi))
  })
  ref_n <- max(cfg$grid_points)
  metrics <- c("peak_logVL", "peak_dpi", "AUC_0_5", "max_positive_slope", "most_negative_slope")
  ref <- results |>
    dplyr::filter(requested_grid_points == ref_n) |>
    dplyr::select(model_key, dplyr::all_of(metrics)) |>
    dplyr::rename_with(~paste0(.x, "_reference"), dplyr::all_of(metrics))
  stability <- results |>
    dplyr::left_join(ref, by = "model_key")
  for (metric in metrics) {
    stability[[paste0("abs_diff_", metric)]] <- abs(stability[[metric]] - stability[[paste0(metric, "_reference")]])
  }
  list(metrics = results, stability = stability)
}

plot_trajectory <- function(data, prediction, facet = TRUE) {
  p <- ggplot2::ggplot(data, ggplot2::aes(dpi, logVL, group = batch)) +
    ggplot2::geom_line(colour = "grey65", alpha = 0.45, linewidth = 0.35) +
    ggplot2::geom_point(colour = "#1A1A1A", size = 1.5) +
    ggplot2::geom_ribbon(
      data = prediction,
      ggplot2::aes(dpi, ymin = lower, ymax = upper),
      inherit.aes = FALSE,
      fill = "#56B4E9",
      alpha = 0.22
    ) +
    ggplot2::geom_line(
      data = prediction,
      ggplot2::aes(dpi, predicted),
      inherit.aes = FALSE,
      colour = "#0072B2",
      linewidth = 0.95
    ) +
    ggplot2::labs(x = "Days post injection (dpi)", y = expression(log[10] * " RNA load")) +
    theme_hunov()
  if (facet && dplyr::n_distinct(data$model_key) > 1L) {
    p <- p + ggplot2::facet_wrap(~inoculum_id, scales = "free_y")
  }
  p
}

plot_residuals <- function(diagnostics) {
  ggplot2::ggplot(diagnostics, ggplot2::aes(fitted, standardized_residual)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey35") +
    ggplot2::geom_point(size = 1.6, colour = "#0072B2", alpha = 0.82) +
    ggplot2::labs(x = "Fitted value", y = "Standardized residual") +
    theme_hunov()
}

plot_qq <- function(diagnostics) {
  ggplot2::ggplot(diagnostics, ggplot2::aes(sample = standardized_residual)) +
    ggplot2::stat_qq(size = 1.4, colour = "#0072B2", alpha = 0.82) +
    ggplot2::stat_qq_line(linewidth = 0.55, colour = "grey30") +
    ggplot2::labs(x = "Theoretical quantile", y = "Standardized residual") +
    theme_hunov()
}

plot_lobo <- function(predictions) {
  ggplot2::ggplot(dplyr::filter(predictions, status == "OK"), ggplot2::aes(logVL, predicted)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.5, colour = "grey35") +
    ggplot2::geom_point(size = 1.7, colour = "#0072B2", alpha = 0.82) +
    ggplot2::labs(x = "Observed log10 RNA load", y = "Predicted log10 RNA load") +
    theme_hunov()
}


derivative_curve <- function(prediction) {
  prediction |>
    dplyr::arrange(model_key, dpi) |>
    dplyr::group_by(model_key, inoculum_id, genotype_ptype) |>
    dplyr::group_modify(~tibble::tibble(
      dpi_mid = (.x$dpi[-1] + .x$dpi[-nrow(.x)]) / 2,
      fitted_slope = diff(.x$predicted) / diff(.x$dpi)
    )) |>
    dplyr::ungroup()
}

random_intercept_blups <- function(bundle) {
  group <- bundle$group_col
  out <- as.data.frame(lme4::ranef(bundle$model)[[group]])
  tibble::rownames_to_column(out, group) |>
    dplyr::rename(random_intercept_BLUP = `(Intercept)`)
}

reconstruction_summary <- function(diagnostics) {
  diagnostics |>
    dplyr::group_by(model_key, inoculum_id, genotype_ptype) |>
    dplyr::summarise(
      n = dplyr::n(),
      RMSE = sqrt(mean(residual^2)),
      MAE = mean(abs(residual)),
      bias = mean(residual),
      correlation = stats::cor(logVL, fitted),
      .groups = "drop"
    )
}

pairwise_feature_contrasts <- function(draws) {
  metrics <- setdiff(names(draws), c("model_key", "inoculum_id", "genotype_ptype", "simulation"))
  keys <- unique(as.character(draws$model_key))
  pairs <- utils::combn(keys, 2, simplify = FALSE)
  purrr::map_dfr(pairs, function(pair) {
    a <- draws |>
      dplyr::filter(as.character(model_key) == pair[1]) |>
      dplyr::select(simulation, dplyr::all_of(metrics)) |>
      dplyr::rename_with(~paste0(.x, "_A"), -simulation)
    b <- draws |>
      dplyr::filter(as.character(model_key) == pair[2]) |>
      dplyr::select(simulation, dplyr::all_of(metrics)) |>
      dplyr::rename_with(~paste0(.x, "_B"), -simulation)
    joined <- dplyr::inner_join(a, b, by = "simulation")
    purrr::map_dfr(metrics, function(metric) {
      difference <- joined[[paste0(metric, "_A")]] - joined[[paste0(metric, "_B")]]
      tibble::tibble(
        inoculum_A = pair[1],
        inoculum_B = pair[2],
        metric = metric,
        median_difference = stats::median(difference, na.rm = TRUE),
        lower_2.5 = stats::quantile(difference, 0.025, na.rm = TRUE, names = FALSE),
        upper_97.5 = stats::quantile(difference, 0.975, na.rm = TRUE, names = FALSE),
        proportion_above_zero = mean(difference > 0, na.rm = TRUE)
      )
    })
  })
}

fixed_curve_difference <- function(bundle, prediction_grid, key_a, key_b) {
  a <- dplyr::filter(prediction_grid, as.character(model_key) == key_a) |>
    dplyr::arrange(dpi)
  b <- dplyr::filter(prediction_grid, as.character(model_key) == key_b) |>
    dplyr::arrange(dpi)
  if (nrow(a) != nrow(b) || !isTRUE(all.equal(a$dpi, b$dpi))) {
    stop("The two prediction grids do not share the same dpi values.", call. = FALSE)
  }
  Xa <- fixed_model_matrix(bundle$model, a)
  Xb <- fixed_model_matrix(bundle$model, b)
  Xd <- Xa - Xb
  beta <- lme4::fixef(bundle$model)
  vc <- positive_definite_vcov(bundle$model)
  estimate <- as.numeric(Xd %*% beta)
  se <- sqrt(pmax(diag(Xd %*% vc %*% t(Xd)), 0))
  tibble::tibble(
    dpi = a$dpi,
    inoculum_A = key_a,
    inoculum_B = key_b,
    difference = estimate,
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se
  )
}

plot_curve_difference <- function(difference) {
  ggplot2::ggplot(difference, ggplot2::aes(dpi, difference)) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey35"
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = "#CC79A7",
      alpha = 0.22
    ) +
    ggplot2::geom_line(linewidth = 0.9, colour = "#7A3E9D") +
    ggplot2::labs(
      x = "Days post injection (dpi)",
      y = "Fitted log10 RNA-load difference"
    ) +
    theme_hunov()
}

run_locked_single_analysis <- function(data, basis_info, cfg, output_root, prefix = "application") {
  dirs <- make_output_dirs(output_root)
  qc <- structural_qc(data, cfg)
  filtered <- apply_structural_qc(data, qc, cfg)
  dat <- filtered$data |>
    dplyr::mutate(batch = factor(batch), inoculum_batch = interaction(model_key, batch, drop = TRUE))

  bundle <- fit_bundle(dat, "single", basis_info, "batch")
  dpi_grid <- make_dpi_grid(cfg$feature_grid_points, cfg$trajectory_window, c(cfg$auc_window, cfg$reference_dpi))
  grid <- new_fixed_grid(dat, dpi_grid, basis_info)
  prediction <- predict_fixed(bundle, grid)
  features <- feature_table(prediction, cfg)
  draws <- simulate_feature_draws(bundle, grid, cfg, seed = cfg$seed + 101L)
  intervals <- summarise_feature_draws(draws)
  lobo <- run_lobo(dat, "single", basis_info, "batch", "batch")
  lobo_summary <- summarise_lobo(lobo)
  diagnostics <- residual_diagnostics(bundle)
  reconstruction <- reconstruction_summary(diagnostics)
  random_intercepts <- random_intercept_blups(bundle)
  derivatives <- derivative_curve(prediction)
  stability <- grid_stability(bundle, dat, cfg)

  write_csv(qc, file.path(dirs$tables, paste0(prefix, "_01_structural_qc.csv")))
  write_csv(filtered$counts, file.path(dirs$tables, paste0(prefix, "_02_batch_counts.csv")))
  status <- model_status(bundle)
  write_csv(status, file.path(dirs$tables, paste0(prefix, "_03_model_status.csv")))
  write_csv(prediction, file.path(dirs$tables, paste0(prefix, "_04_fitted_curve.csv")))
  write_csv(features, file.path(dirs$tables, paste0(prefix, "_05_feature_estimates.csv")))
  write_csv(intervals, file.path(dirs$tables, paste0(prefix, "_06_feature_intervals.csv")))
  write_csv(draws, file.path(dirs$tables, paste0(prefix, "_06b_feature_simulation_draws.csv")))
  write_csv(lobo, file.path(dirs$tables, paste0(prefix, "_07_lobo_predictions.csv")))
  write_csv(lobo_summary, file.path(dirs$tables, paste0(prefix, "_08_lobo_summary.csv")))
  write_csv(diagnostics, file.path(dirs$tables, paste0(prefix, "_09_residual_diagnostics.csv")))
  write_csv(reconstruction, file.path(dirs$tables, paste0(prefix, "_10_reconstruction_summary.csv")))
  write_csv(random_intercepts, file.path(dirs$tables, paste0(prefix, "_11_random_intercept_BLUPs.csv")))
  write_csv(derivatives, file.path(dirs$tables, paste0(prefix, "_12_derivative_curve.csv")))
  write_csv(stability$metrics, file.path(dirs$tables, paste0(prefix, "_13_grid_metrics.csv")))
  write_csv(stability$stability, file.path(dirs$tables, paste0(prefix, "_14_grid_stability.csv")))
  saveRDS(bundle, file.path(dirs$objects, paste0(prefix, "_model_bundle.rds")))

  save_tiff(
    plot_trajectory(dat, prediction, facet = FALSE),
    file.path(dirs$plots, paste0(prefix, "_trajectory.tiff")),
    7.5, 5.5
  )
  save_tiff(plot_residuals(diagnostics), file.path(dirs$plots, paste0(prefix, "_residuals.tiff")), 6, 4.8)
  save_tiff(plot_qq(diagnostics), file.path(dirs$plots, paste0(prefix, "_qq.tiff")), 6, 4.8)
  save_tiff(plot_lobo(lobo), file.path(dirs$plots, paste0(prefix, "_lobo.tiff")), 6, 5)

  list(
    data = dat,
    qc = qc,
    bundle = bundle,
    model_status = status,
    prediction = prediction,
    features = features,
    draws = draws,
    intervals = intervals,
    lobo = lobo,
    lobo_summary = lobo_summary,
    diagnostics = diagnostics,
    reconstruction = reconstruction,
    random_intercepts = random_intercepts,
    derivatives = derivatives,
    grid_stability = stability
  )
}
