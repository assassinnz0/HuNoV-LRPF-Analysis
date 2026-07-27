# HuNoV-LRPF publication figures
# Plot constructors are separated from statistical estimation so that the main
# analysis remains readable and every figure can be regenerated from source data.

require_hunov_plot_packages <- function() {
  packages <- c("ggplot2", "dplyr", "tidyr", "forcats", "patchwork", "scales")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing plotting packages: ", paste(missing, collapse = ", "),
      ". Install them before running the analysis.",
      call. = FALSE
    )
  }
}


hunov_inoculum_palette <- c(
  "GZJY315" = "#0072B2",
  "GZJY324" = "#E69F00",
  "GZJY17"  = "#009E73",
  "GZJY19"  = "#D55E00",
  "GZJY29"  = "#CC79A7",
  "GZJY122" = "#56B4E9"
)

hunov_basis_palette <- c(
  "B1" = "#D55E00",
  "B2" = "#009E73",
  "B3" = "#0072B2",
  "B4" = "#CC79A7"
)

hunov_rate_palette <- c(
  "Maximum positive slope" = "#0072B2",
  "Most negative slope" = "#D55E00"
)


# Figure 4 retains the endpoint-summary design and palette used in legacy v9.3.
hunov_endpoint_palette_v93 <- c(
  "GZJY315" = "#2C7FB8",
  "GZJY324" = "#7FCDBB",
  "GZJY17"  = "#41B6C4",
  "GZJY19"  = "#F28E2B",
  "GZJY29"  = "#6A51A3"
)

publication_multiline_label <- function(inoculum_id, genotype_ptype) {
  paste0(as.character(inoculum_id), "\n", as.character(genotype_ptype))
}

publication_label <- function(inoculum_id, genotype_ptype) {
  paste0(as.character(inoculum_id), " | ", as.character(genotype_ptype))
}

add_publication_label <- function(data) {
  data |>
    dplyr::mutate(
      display_label = publication_label(inoculum_id, genotype_ptype),
      display_label = factor(display_label, levels = unique(display_label))
    )
}


plot_data_coverage <- function(data) {
  dat <- add_publication_label(data) |>
    dplyr::mutate(
      batch = factor(as.character(batch), levels = rev(sort(unique(as.character(batch))))),
      tile_text_colour = dplyr::if_else(logVL < 3, "white", "grey10")
    )

  ggplot2::ggplot(dat, ggplot2::aes(x = dpi, y = batch, fill = logVL)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f", logVL), colour = tile_text_colour),
      size = 2.45
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_viridis_c(option = "viridis", end = 0.96) +
    ggplot2::facet_wrap(~display_label, ncol = 2, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = seq(0, 7, by = 1), expand = c(0, 0)) +
    ggplot2::labs(
      x = "Days post injection (dpi)",
      y = "Batch",
      fill = expression(log[10] * " RNA load")
    ) +
    theme_hunov(base_size = 10) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7.5),
      strip.text = ggplot2::element_text(size = 9, face = "bold"),
      panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.4)
    )
}


plot_bspline_basis <- function(basis_info, trajectory_window = c(0, 7), n = 401L) {
  dpi <- seq(trajectory_window[1], trajectory_window[2], length.out = n)
  basis <- basis_matrix(dpi, basis_info)
  dat <- dplyr::bind_cols(tibble::tibble(dpi = dpi), basis) |>
    tidyr::pivot_longer(-dpi, names_to = "basis", values_to = "value")

  ggplot2::ggplot(dat, ggplot2::aes(dpi, value, colour = basis, linetype = basis)) +
    ggplot2::geom_vline(
      xintercept = basis_info$knots,
      linetype = 3,
      linewidth = 0.5,
      colour = "grey35"
    ) +
    ggplot2::geom_line(linewidth = 0.95) +
    ggplot2::scale_colour_manual(values = hunov_basis_palette, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = seq(trajectory_window[1], trajectory_window[2], 1)) +
    ggplot2::labs(
      x = "Days post injection (dpi)",
      y = "B-spline basis value",
      colour = "Basis",
      linetype = "Basis"
    ) +
    theme_hunov() +
    ggplot2::theme(legend.position = "bottom")
}

plot_development_design <- function(data, basis_info, trajectory_window = c(0, 7)) {
  # Figure 2 uses an explicit two-row design:
  #   A. full-width development-data coverage panel
  #   B. centered spline-basis panel with sufficient width for axis and legend text
  # Manual panel tags are used so that the spacer columns do not affect labeling.

  p_coverage <- plot_data_coverage(data) +
    ggplot2::labs(tag = "A") +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 12),
      plot.tag.position = c(0, 1),
      plot.margin = ggplot2::margin(t = 8, r = 8, b = 4, l = 8)
    )

  p_basis <- plot_bspline_basis(basis_info, trajectory_window) +
    ggplot2::labs(tag = "B") +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 12),
      plot.tag.position = c(0, 1),
      plot.margin = ggplot2::margin(t = 4, r = 8, b = 8, l = 8)
    )

  basis_row <- patchwork::plot_spacer() + p_basis + patchwork::plot_spacer() +
    patchwork::plot_layout(widths = c(0.40, 1.20, 0.40))

  p_coverage / basis_row +
    patchwork::plot_layout(heights = c(2.60, 1.15)) &
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}


plot_trajectory_overlay <- function(data, prediction) {
  dat <- add_publication_label(data)
  pred <- add_publication_label(prediction)

  ggplot2::ggplot(dat, ggplot2::aes(dpi, logVL, group = interaction(display_label, batch))) +
    ggplot2::geom_point(ggplot2::aes(colour = inoculum_id), alpha = 0.55, size = 1.3) +
    ggplot2::geom_ribbon(
      data = pred,
      ggplot2::aes(x = dpi, ymin = lower, ymax = upper, group = display_label, fill = inoculum_id),
      inherit.aes = FALSE,
      alpha = 0.12,
      colour = NA
    ) +
    ggplot2::geom_line(
      data = pred,
      ggplot2::aes(x = dpi, y = predicted, colour = inoculum_id, linetype = display_label),
      inherit.aes = FALSE,
      linewidth = 0.95
    ) +
    ggplot2::scale_colour_manual(values = hunov_inoculum_palette) +
    ggplot2::scale_fill_manual(values = hunov_inoculum_palette) +
    ggplot2::scale_x_continuous(breaks = 0:7) +
    ggplot2::labs(
      x = "Days post injection (dpi)",
      y = expression(log[10] * " RNA load (copies/" * mu * "L)"),
      colour = "Inoculum",
      linetype = "Inoculum"
    ) +
    theme_hunov() +
    ggplot2::theme(legend.position = "bottom")
}


plot_trajectory_facets_publication <- function(data, prediction) {
  dat <- add_publication_label(data)
  pred <- add_publication_label(prediction)

  ggplot2::ggplot(dat, ggplot2::aes(dpi, logVL, group = batch)) +
    ggplot2::geom_line(colour = "grey65", alpha = 0.42, linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(colour = inoculum_id), size = 1.35, alpha = 0.82) +
    ggplot2::geom_ribbon(
      data = pred,
      ggplot2::aes(x = dpi, ymin = lower, ymax = upper, fill = inoculum_id),
      inherit.aes = FALSE,
      alpha = 0.20,
      colour = NA
    ) +
    ggplot2::geom_line(
      data = pred,
      ggplot2::aes(x = dpi, y = predicted, colour = inoculum_id),
      inherit.aes = FALSE,
      linewidth = 0.95
    ) +
    ggplot2::scale_colour_manual(values = hunov_inoculum_palette, guide = "none") +
    ggplot2::scale_fill_manual(values = hunov_inoculum_palette, guide = "none") +
    ggplot2::facet_wrap(~display_label, ncol = 2) +
    ggplot2::scale_x_continuous(breaks = 0:7) +
    ggplot2::labs(
      x = "Days post injection (dpi)",
      y = expression(log[10] * " RNA load (copies/" * mu * "L)")
    ) +
    theme_hunov(base_size = 10.5)
}


plot_derivative_curves_publication <- function(derivatives) {
  dat <- add_publication_label(derivatives)
  display_levels <- levels(dat$display_label)
  inoculum_by_display <- dat |>
    dplyr::distinct(display_label, inoculum_id) |>
    dplyr::mutate(display_label = as.character(display_label))
  display_palette <- stats::setNames(
    hunov_inoculum_palette[as.character(inoculum_by_display$inoculum_id)],
    inoculum_by_display$display_label
  )

  ggplot2::ggplot(dat, ggplot2::aes(
    dpi_mid, fitted_slope,
    colour = display_label,
    linetype = display_label
  )) +
    ggplot2::geom_hline(yintercept = 0, linetype = 3, linewidth = 0.4, colour = "grey35") +
    ggplot2::geom_line(linewidth = 0.95) +
    ggplot2::scale_colour_manual(
      name = "Inoculum",
      values = display_palette,
      breaks = display_levels,
      drop = FALSE
    ) +
    ggplot2::scale_linetype_discrete(
      name = "Inoculum",
      breaks = display_levels,
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(breaks = 0:7) +
    ggplot2::labs(
      x = "Days post injection (dpi)",
      y = expression("Fitted slope (" * log[10] * " copies/" * mu * "L/day)")
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    theme_hunov() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    )
}

metric_label_map <- c(
  peak_logVL = "Peak log10 RNA load",
  peak_dpi = "Fitted peak time",
  AUC_0_5 = "0-5 dpi log-scale AUC",
  max_positive_slope = "Maximum positive fitted slope",
  most_negative_slope = "Most negative fitted slope",
  fitted_0dpi_logVL = "Fitted 0-dpi level",
  max_net_increase = "Maximum increase above 0 dpi",
  centered_AUC_0_5 = "0-5 dpi centered AUC"
)

plot_feature_forest <- function(intervals, metrics = names(metric_label_map)[1:5]) {
  dat <- intervals |>
    dplyr::filter(metric %in% metrics) |>
    add_publication_label() |>
    dplyr::mutate(
      metric_label = unname(metric_label_map[metric]),
      metric_label = factor(metric_label, levels = unname(metric_label_map[metrics]))
    )

  ggplot2::ggplot(dat, ggplot2::aes(median, display_label)) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = lower_2.5, xmax = upper_97.5),
      height = 0.18,
      linewidth = 0.45
    ) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~metric_label, scales = "free_x", ncol = 2) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_hunov(base_size = 10.5) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
}


plot_centered_trajectories <- function(prediction, reference_dpi = 0) {
  dat <- prediction |>
    dplyr::group_by(model_key, inoculum_id, genotype_ptype) |>
    dplyr::mutate(
      reference = predicted[which.min(abs(dpi - reference_dpi))],
      centered = predicted - reference,
      lower_centered = lower - reference,
      upper_centered = upper - reference
    ) |>
    dplyr::ungroup() |>
    add_publication_label()

  ggplot2::ggplot(dat, ggplot2::aes(dpi, centered)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 3, linewidth = 0.4, colour = "grey35") +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower_centered, ymax = upper_centered, fill = inoculum_id),
      alpha = 0.18,
      colour = NA
    ) +
    ggplot2::geom_line(ggplot2::aes(colour = inoculum_id), linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = hunov_inoculum_palette, guide = "none") +
    ggplot2::scale_fill_manual(values = hunov_inoculum_palette, guide = "none") +
    ggplot2::facet_wrap(~display_label, ncol = 2) +
    ggplot2::scale_x_continuous(breaks = 0:7) +
    ggplot2::labs(
      x = "Days post injection (dpi)",
      y = expression("Change from fitted 0 dpi (" * log[10] * " units)")
    ) +
    theme_hunov(base_size = 10.5)
}


plot_lobo_observed_predicted <- function(predictions, facet = FALSE) {
  dat <- predictions |>
    dplyr::filter(status == "OK") |>
    add_publication_label()

  p <- ggplot2::ggplot(dat, ggplot2::aes(predicted, logVL)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.5, colour = "grey35") +
    ggplot2::geom_point(ggplot2::aes(colour = inoculum_id), size = 1.5, alpha = 0.78) +
    ggplot2::scale_colour_manual(values = hunov_inoculum_palette, guide = "none") +
    ggplot2::labs(
      x = expression("LOBO-predicted " * log[10] * " RNA load"),
      y = expression("Observed " * log[10] * " RNA load")
    ) +
    theme_hunov()
  if (facet && dplyr::n_distinct(dat$model_key) > 1L) {
    p <- p + ggplot2::facet_wrap(~display_label, ncol = 2)
  }
  p
}


plot_lobo_residual_time <- function(predictions, facet = FALSE) {
  dat <- predictions |>
    dplyr::filter(status == "OK") |>
    dplyr::mutate(cv_residual = logVL - predicted) |>
    add_publication_label()

  p <- ggplot2::ggplot(dat, ggplot2::aes(dpi, cv_residual)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey35") +
    ggplot2::geom_point(ggplot2::aes(colour = inoculum_id), size = 1.4, alpha = 0.76) +
    ggplot2::geom_smooth(
      method = "loess", formula = y ~ x, se = FALSE,
      linewidth = 0.75, colour = "#0072B2"
    ) +
    ggplot2::scale_colour_manual(values = hunov_inoculum_palette, guide = "none") +
    ggplot2::scale_x_continuous(breaks = 0:7) +
    ggplot2::labs(x = "Days post injection (dpi)", y = "LOBO residual") +
    theme_hunov()
  if (facet && dplyr::n_distinct(dat$model_key) > 1L) {
    p <- p + ggplot2::facet_wrap(~display_label, ncol = 2)
  }
  p
}

plot_lobo_combined <- function(predictions, facet = FALSE) {
  plot_lobo_observed_predicted(predictions, facet = facet) +
    plot_lobo_residual_time(predictions, facet = facet) +
    patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(tag_levels = "A")
}


plot_pooled_diagnostics <- function(diagnostics) {
  dat <- add_publication_label(diagnostics)

  p1 <- ggplot2::ggplot(dat, ggplot2::aes(fitted, standardized_residual)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey35") +
    ggplot2::geom_point(ggplot2::aes(colour = inoculum_id), size = 1.35, alpha = 0.72) +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.75, colour = "#0072B2") +
    ggplot2::scale_colour_manual(values = hunov_inoculum_palette, guide = "none") +
    ggplot2::labs(x = "Fitted value", y = "Standardized residual") +
    theme_hunov()

  p2 <- ggplot2::ggplot(dat, ggplot2::aes(sample = standardized_residual)) +
    ggplot2::stat_qq(ggplot2::aes(colour = inoculum_id), size = 1.25, alpha = 0.75) +
    ggplot2::stat_qq_line(linewidth = 0.5, colour = "grey30") +
    ggplot2::scale_colour_manual(values = hunov_inoculum_palette, guide = "none") +
    ggplot2::labs(x = "Theoretical quantile", y = "Standardized residual") +
    theme_hunov()

  p3 <- ggplot2::ggplot(dat, ggplot2::aes(fitted, logVL)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.5, colour = "grey35") +
    ggplot2::geom_point(ggplot2::aes(colour = inoculum_id), size = 1.35, alpha = 0.72) +
    ggplot2::scale_colour_manual(values = hunov_inoculum_palette, guide = "none") +
    ggplot2::labs(x = "Fitted value", y = "Observed log10 RNA load") +
    theme_hunov()

  p4 <- ggplot2::ggplot(dat, ggplot2::aes(dpi, standardized_residual)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey35") +
    ggplot2::geom_point(ggplot2::aes(colour = inoculum_id), size = 1.35, alpha = 0.72) +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.75, colour = "#0072B2") +
    ggplot2::scale_colour_manual(values = hunov_inoculum_palette, guide = "none") +
    ggplot2::scale_x_continuous(breaks = 0:7) +
    ggplot2::labs(x = "Days post injection (dpi)", y = "Standardized residual") +
    theme_hunov()

  (p1 + p2) / (p3 + p4) + patchwork::plot_annotation(tag_levels = "A")
}


plot_inoculum_diagnostics <- function(diagnostics) {
  # Supplementary Figure S3 is intentionally monochrome. The inoculum identity
  # is carried by the facet strips, so colour is unnecessary and can distract
  # from the residual patterns being assessed.
  dat <- add_publication_label(diagnostics)

  p1 <- ggplot2::ggplot(dat, ggplot2::aes(fitted, standardized_residual)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = "grey40") +
    ggplot2::geom_point(size = 0.9, alpha = 0.72, colour = "black") +
    ggplot2::geom_smooth(
      method = "loess", formula = y ~ x, se = FALSE,
      linewidth = 0.55, colour = "grey20"
    ) +
    ggplot2::facet_wrap(~display_label, nrow = 1) +
    ggplot2::labs(x = "Fitted value", y = "Standardized residual") +
    theme_hunov(base_size = 8.5)

  p2 <- ggplot2::ggplot(dat, ggplot2::aes(sample = standardized_residual)) +
    ggplot2::stat_qq(size = 0.9, alpha = 0.72, colour = "black") +
    ggplot2::stat_qq_line(linewidth = 0.4, colour = "grey30") +
    ggplot2::facet_wrap(~display_label, nrow = 1) +
    ggplot2::labs(x = "Theoretical quantile", y = "Standardized residual") +
    theme_hunov(base_size = 8.5)

  p3 <- ggplot2::ggplot(dat, ggplot2::aes(fitted, logVL)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4, colour = "grey40") +
    ggplot2::geom_point(size = 0.9, alpha = 0.72, colour = "black") +
    ggplot2::facet_wrap(~display_label, nrow = 1) +
    ggplot2::labs(x = "Fitted value", y = "Observed log10 RNA load") +
    theme_hunov(base_size = 8.5)

  p4 <- ggplot2::ggplot(dat, ggplot2::aes(dpi, standardized_residual)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = "grey40") +
    ggplot2::geom_point(size = 0.9, alpha = 0.72, colour = "black") +
    ggplot2::geom_smooth(
      method = "loess", formula = y ~ x, se = FALSE,
      linewidth = 0.55, colour = "grey20"
    ) +
    ggplot2::scale_x_continuous(breaks = 0:7) +
    ggplot2::facet_wrap(~display_label, nrow = 1) +
    ggplot2::labs(x = "Days post injection (dpi)", y = "Standardized residual") +
    theme_hunov(base_size = 8.5)

  p1 / p2 / p3 / p4 +
    patchwork::plot_layout(heights = c(1, 1, 1, 1)) +
    patchwork::plot_annotation(tag_levels = "A")
}

plot_variance_components <- function(model_status_table) {
  dat <- tibble::tibble(
    component = c("Batch random intercept", "Residual"),
    variance = c(
      model_status_table$random_intercept_variance[1],
      model_status_table$residual_sd[1]^2
    )
  )
  ggplot2::ggplot(dat, ggplot2::aes(component, variance)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::labs(x = NULL, y = "Estimated variance") +
    theme_hunov() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 15, hjust = 1))
}


plot_endpoint_bubble <- function(features) {
  # Figure 4A follows the legacy v9.3 peak-time/AUC bubble implementation.
  dat <- features |>
    dplyr::mutate(
      inoculum_id = factor(as.character(inoculum_id), levels = unique(as.character(inoculum_id))),
      label = as.character(inoculum_id)
    )

  bubble_x_breaks <- seq(3.1, 3.7, by = 0.1)
  bubble_x_limits <- range(bubble_x_breaks)
  bubble_y_limits <- range(dat$peak_logVL, na.rm = TRUE) + c(-0.60, 0.45)
  palette <- hunov_endpoint_palette_v93[levels(dat$inoculum_id)]

  ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = peak_dpi,
      y = peak_logVL,
      size = AUC_0_5,
      fill = inoculum_id
    )
  ) +
    ggplot2::geom_point(shape = 21, colour = "white", stroke = 0.4, alpha = 0.83) +
    ggplot2::geom_text(
      ggplot2::aes(y = peak_logVL - 0.22, label = label),
      size = 3.8,
      vjust = 1,
      hjust = 0.5,
      colour = "grey20",
      show.legend = FALSE
    ) +
    ggplot2::annotate(
      "text",
      x = bubble_x_limits[1] + 0.02,
      y = bubble_y_limits[1] + 0.18,
      label = "Bubble size = 0–5-dpi log-scale AUC",
      hjust = 0,
      colour = "#6B7A90",
      size = 3.6
    ) +
    ggplot2::scale_size_continuous(range = c(8, 18), guide = "none") +
    ggplot2::scale_fill_manual(values = palette, guide = "none") +
    ggplot2::scale_x_continuous(breaks = bubble_x_breaks, minor_breaks = NULL) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 6), minor_breaks = NULL) +
    ggplot2::coord_cartesian(xlim = bubble_x_limits, ylim = bubble_y_limits, clip = "off") +
    ggplot2::labs(
      x = "Fitted peak time (dpi)",
      y = expression("Peak " * log[10] * " RNA load")
    ) +
    theme_hunov(base_size = 13) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(colour = "grey88", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 28, 10, 10)
    )
}


plot_change_rate_lollipop <- function(features) {
  # Figure 4B follows the legacy v9.3 connected increase/decline display.
  endpoint <- features |>
    dplyr::mutate(
      inoculum_id = factor(as.character(inoculum_id), levels = unique(as.character(inoculum_id))),
      display_multiline = publication_multiline_label(inoculum_id, genotype_ptype)
    ) |>
    dplyr::arrange(inoculum_id)

  segment_data <- endpoint |>
    dplyr::transmute(
      inoculum_id,
      x_decline = most_negative_slope,
      x_increase = max_positive_slope
    )

  point_data <- endpoint |>
    dplyr::select(inoculum_id, max_positive_slope, most_negative_slope) |>
    tidyr::pivot_longer(
      cols = c(max_positive_slope, most_negative_slope),
      names_to = "rate_type",
      values_to = "rate"
    ) |>
    dplyr::mutate(
      rate_type = dplyr::recode(
        rate_type,
        max_positive_slope = "Maximum positive slope",
        most_negative_slope = "Most negative slope"
      ),
      rate_type = factor(
        rate_type,
        levels = c("Maximum positive slope", "Most negative slope")
      )
    )

  y_labels <- stats::setNames(endpoint$display_multiline, as.character(endpoint$inoculum_id))

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = segment_data,
      ggplot2::aes(
        x = x_decline,
        xend = x_increase,
        y = inoculum_id,
        yend = inoculum_id
      ),
      colour = "grey82",
      linewidth = 1.4,
      alpha = 0.95
    ) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.8) +
    ggplot2::geom_point(
      data = point_data,
      ggplot2::aes(x = rate, y = inoculum_id, colour = rate_type),
      size = 4.0,
      alpha = 0.96
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Maximum positive slope" = "#2C7FB8",
        "Most negative slope" = "#F28E2B"
      ),
      breaks = c("Maximum positive slope", "Most negative slope")
    ) +
    ggplot2::scale_y_discrete(
      limits = rev(levels(endpoint$inoculum_id)),
      labels = y_labels
    ) +
    ggplot2::scale_x_continuous(breaks = seq(-2, 3, by = 1), minor_breaks = NULL) +
    ggplot2::labs(
      x = expression(log[10] * " RNA-load change per day"),
      y = NULL,
      colour = NULL
    ) +
    theme_hunov(base_size = 13) +
    ggplot2::theme(
      legend.position = "right",
      legend.justification = "center",
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.box.margin = ggplot2::margin(0, 0, 0, 10),
      plot.margin = ggplot2::margin(10, 28, 10, 10)
    )
}


plot_standardized_endpoint_heatmap <- function(features) {
  # Figure 4C follows the legacy v9.3 four-endpoint standardized heatmap.
  endpoint <- features |>
    dplyr::mutate(
      inoculum_id = factor(as.character(inoculum_id), levels = unique(as.character(inoculum_id))),
      display_multiline = publication_multiline_label(inoculum_id, genotype_ptype)
    )

  dat <- endpoint |>
    dplyr::transmute(
      inoculum_id,
      Peak = peak_logVL,
      AUC = AUC_0_5,
      `Positive\nslope` = max_positive_slope,
      `Decline\nstrength` = abs(most_negative_slope)
    ) |>
    tidyr::pivot_longer(
      cols = -inoculum_id,
      names_to = "endpoint",
      values_to = "raw_value"
    ) |>
    dplyr::group_by(endpoint) |>
    dplyr::mutate(
      endpoint_mean = mean(raw_value, na.rm = TRUE),
      endpoint_sd = stats::sd(raw_value, na.rm = TRUE),
      z_score = dplyr::if_else(
        is.finite(endpoint_sd) & endpoint_sd > 0,
        (raw_value - endpoint_mean) / endpoint_sd,
        0
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      endpoint = factor(
        endpoint,
        levels = c("Peak", "AUC", "Positive\nslope", "Decline\nstrength")
      ),
      z_label = sprintf("%+.1f", z_score)
    )

  y_labels <- stats::setNames(endpoint$display_multiline, as.character(endpoint$inoculum_id))
  heatmap_limit <- max(abs(dat$z_score), na.rm = TRUE)
  heatmap_limit <- max(1.5, ceiling(heatmap_limit * 10) / 10)

  ggplot2::ggplot(dat, ggplot2::aes(x = endpoint, y = inoculum_id, fill = z_score)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = z_label), size = 3.7, colour = "black") +
    ggplot2::scale_y_discrete(
      limits = rev(levels(endpoint$inoculum_id)),
      labels = y_labels
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-heatmap_limit, heatmap_limit),
      name = "Z-score"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_hunov(base_size = 13) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = "right"
    )
}


plot_endpoint_composite <- function(features) {
  # Legacy v9.3 Figure 4 arrangement:
  #   A. peak-time/AUC bubble plot
  #   B. connected positive/negative rate plot
  #   C. standardized four-endpoint heatmap
  endpoint_bubble_panel <- plot_endpoint_bubble(features) +
    ggplot2::labs(tag = "A") +
    ggplot2::theme(plot.tag.position = c(0.01, 0.99))

  endpoint_rate_panel <- plot_change_rate_lollipop(features) +
    ggplot2::labs(tag = "B") +
    ggplot2::theme(
      plot.tag.position = c(0.01, 0.99),
      legend.position = "right"
    )

  endpoint_heatmap_panel <- plot_standardized_endpoint_heatmap(features) +
    ggplot2::labs(tag = "C") +
    ggplot2::theme(plot.tag.position = c(0.01, 0.99))

  endpoint_bubble_panel /
    (endpoint_rate_panel | endpoint_heatmap_panel) +
    patchwork::plot_layout(heights = c(1.65, 1), widths = c(1, 1.15)) +
    patchwork::plot_annotation()
}

plot_locked_application_composite <- function(data, prediction, lobo, diagnostics) {
  # Main Figure 5 is intentionally monochrome. Line weight, point shape, grey
  # bands, and reference lines carry the visual hierarchy.
  lobo_data <- dplyr::filter(lobo, status == "OK")

  p1 <- ggplot2::ggplot(data, ggplot2::aes(dpi, logVL, group = batch)) +
    ggplot2::geom_line(colour = "grey65", alpha = 0.50, linewidth = 0.38) +
    ggplot2::geom_point(colour = "black", size = 1.55, alpha = 0.82) +
    ggplot2::geom_ribbon(
      data = prediction,
      ggplot2::aes(dpi, ymin = lower, ymax = upper),
      inherit.aes = FALSE,
      fill = "grey75",
      alpha = 0.45,
      colour = NA
    ) +
    ggplot2::geom_line(
      data = prediction,
      ggplot2::aes(dpi, predicted),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 1.0
    ) +
    ggplot2::scale_x_continuous(breaks = 0:7) +
    ggplot2::labs(
      x = "Days post injection (dpi)",
      y = expression(log[10] * " RNA load")
    ) +
    theme_hunov()

  p2 <- ggplot2::ggplot(lobo_data, ggplot2::aes(predicted, logVL)) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0,
      linetype = 2, linewidth = 0.5, colour = "grey35"
    ) +
    ggplot2::geom_point(
      shape = 21, fill = "white", colour = "black",
      stroke = 0.55, size = 1.8, alpha = 0.92
    ) +
    ggplot2::labs(
      x = expression("LOBO-predicted " * log[10] * " RNA load"),
      y = expression("Observed " * log[10] * " RNA load")
    ) +
    theme_hunov()

  p3 <- ggplot2::ggplot(diagnostics, ggplot2::aes(fitted, standardized_residual)) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = 2,
      linewidth = 0.4, colour = "grey35"
    ) +
    ggplot2::geom_point(
      shape = 21, fill = "white", colour = "black",
      stroke = 0.5, size = 1.7, alpha = 0.90
    ) +
    ggplot2::labs(x = "Fitted value", y = "Standardized residual") +
    theme_hunov()

  p4 <- ggplot2::ggplot(diagnostics, ggplot2::aes(sample = standardized_residual)) +
    ggplot2::stat_qq(
      shape = 21, fill = "white", colour = "black",
      stroke = 0.5, size = 1.55, alpha = 0.90
    ) +
    ggplot2::stat_qq_line(linewidth = 0.55, colour = "grey30") +
    ggplot2::labs(x = "Theoretical quantile", y = "Standardized residual") +
    theme_hunov()

  (p1 + p2) / (p3 + p4) +
    patchwork::plot_annotation(tag_levels = "A")
}

