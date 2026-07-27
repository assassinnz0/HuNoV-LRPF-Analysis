# HuNoV-LRPF output and reproducibility helpers

make_publication_output_dirs <- function(root) {
  dirs <- list(
    root = root,
    main_tables = file.path(root, "main_tables"),
    supplementary_tables = file.path(root, "supplementary_tables"),
    main_figures = file.path(root, "main_figures"),
    supplementary_figures = file.path(root, "supplementary_figures"),
    source_data = file.path(root, "source_data"),
    model_objects = file.path(root, "model_objects"),
    logs = file.path(root, "logs")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}

write_text <- function(lines, path) {
  writeLines(lines, con = path, useBytes = TRUE)
  invisible(path)
}

write_model_summary <- function(model, path) {
  write_text(capture.output(summary(model)), path)
}

write_session_info <- function(path) {
  write_text(capture.output(utils::sessionInfo()), path)
}

file_checksum_table <- function(paths, labels = basename(paths)) {
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  tibble::tibble(
    label = labels,
    file = basename(paths),
    md5 = unname(tools::md5sum(paths)),
    bytes = as.numeric(file.info(paths)$size),
    modified_time = format(file.info(paths)$mtime, "%Y-%m-%d %H:%M:%S %Z")
  )
}

save_figure_set <- function(figures, directory) {
  purrr::iwalk(figures, function(spec, name) {
    if (!is.list(spec) || is.null(spec$plot)) {
      stop("Each figure specification must contain a plot element: ", name, call. = FALSE)
    }
    width <- if (is.null(spec$width)) 8 else spec$width
    height <- if (is.null(spec$height)) 6 else spec$height
    save_tiff(
      spec$plot,
      file.path(directory, paste0(name, ".tiff")),
      width = width,
      height = height
    )
  })
  invisible(names(figures))
}

write_named_tables <- function(tables, directory) {
  purrr::iwalk(tables, function(x, name) {
    write_csv(x, file.path(directory, paste0(name, ".csv")))
  })
  invisible(names(tables))
}

write_output_manifest <- function(root, path = file.path(root, "OUTPUT_MANIFEST.csv")) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]

  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  file_norm <- normalizePath(files, winslash = "/", mustWork = TRUE)
  keep <- file_norm != path_norm
  files <- files[keep]
  file_norm <- file_norm[keep]

  root_norm <- normalizePath(root, winslash = "/", mustWork = TRUE)
  rel <- sub(paste0("^", root_norm, "/?"), "", file_norm)
  manifest <- tibble::tibble(
    relative_path = rel,
    bytes = as.numeric(file.info(files)$size),
    md5 = unname(tools::md5sum(files))
  ) |>
    dplyr::arrange(relative_path)
  write_csv(manifest, path)
  manifest
}

write_run_metadata <- function(cfg, development_data, application_data, selected_basis, output_path, version = "10.3") {
  metadata <- tibble::tibble(
    item = c(
      "analysis_version",
      "primary_model_family",
      "selected_spline_df",
      "internal_knots",
      "boundary_knots",
      "trajectory_window",
      "auc_window",
      "feature_grid_points",
      "grid_stability_points",
      "simulation_draws",
      "development_observations",
      "development_batches",
      "application_observations",
      "application_batches",
      "R_version",
      "run_timestamp"
    ),
    value = c(
      version,
      "M2",
      as.character(selected_basis$spline_df),
      paste(selected_basis$knots, collapse = ";"),
      paste(selected_basis$boundary_knots, collapse = ";"),
      paste(cfg$trajectory_window, collapse = ";"),
      paste(cfg$auc_window, collapse = ";"),
      as.character(cfg$feature_grid_points),
      paste(cfg$grid_points, collapse = ";"),
      as.character(cfg$n_sim),
      as.character(nrow(development_data)),
      as.character(dplyr::n_distinct(development_data$inoculum_batch)),
      as.character(nrow(application_data)),
      as.character(dplyr::n_distinct(application_data$batch)),
      R.version.string,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
  )
  write_csv(metadata, output_path)
}
