#' Plot standardized CPUE index and CV
#'
#' Creates a two-panel plot showing the standardized CPUE index and its
#' coefficient of variation (CV) over time.
#'
#' @param index_df A data.frame returned by [get_index()] containing at least
#'   `time`, `index`, and `cv`. If it also contains `region`, each region is
#'   plotted separately.
#' @param time_values Optional vector of labels for the x-axis. Must have the
#'   same length as the number of unique time values in `index_df`.
#' @param time_positions Optional numeric vector of x-axis positions. Must have
#'   the same length as the number of unique time values in `index_df`.
#' @param x_text_angle Optional numeric angle for x-axis text. If `NULL`,
#'   `plot_index()` uses `90` for year-month-like labels and `45` otherwise.
#' @param x_text_hjust Optional horizontal justification for x-axis text. If
#'   `NULL`, a value matching `x_text_angle` is chosen automatically.
#' @param x_text_vjust Optional vertical justification for x-axis text. If
#'   `NULL`, a value matching `x_text_angle` is chosen automatically.
#'
#' @return A ggplot object.
#' @export
plot_index <- function(
    index_df,
    time_values = NULL,
    time_positions = NULL,
    x_text_angle = NULL,
    x_text_hjust = NULL,
    x_text_vjust = NULL
) {
  req <- c("time", "index", "cv")
  miss <- setdiff(req, names(index_df))
  if (length(miss)) {
    stop(
      "`index_df` is missing required columns: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  time_levels <- sort(unique(index_df$time))
  n_time <- length(time_levels)
  if (!is.null(time_values) && length(time_values) != n_time) {
    stop(
      "`time_values` must have the same length as the number of unique `time` values.",
      call. = FALSE
    )
  }

  if (!is.null(time_positions)) {
    if (!is.numeric(time_positions) || length(time_positions) != n_time) {
      stop(
        "`time_positions` must be a numeric vector with the same length as the number of unique `time` values.",
        call. = FALSE
      )
    }
  }

  time_map <- stats::setNames(
    if (is.null(time_positions)) seq_len(n_time) else time_positions,
    as.character(time_levels)
  )
  time_seq <- unname(time_map[as.character(index_df$time)])
  region <- if ("region" %in% names(index_df)) as.character(index_df$region) else "total"

  plot_dat <- rbind(
    data.frame(
      time = time_seq,
      region = region,
      value = index_df$index,
      panel = "Standardized CPUE"
    ),
    data.frame(
      time = time_seq,
      region = region,
      value = index_df$cv,
      panel = "CV"
    )
  )

  plot_dat$panel <- factor(
    plot_dat$panel,
    levels = c("Standardized CPUE", "CV")
  )

  plot_dat$region <- factor(plot_dat$region, levels = unique(region))

  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = time, y = value, group = region)) +
    ggplot2::geom_line(
      data = subset(plot_dat, panel == "Standardized CPUE"),
      colour = "#0F766E",
      linewidth = 1.1,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = subset(plot_dat, panel == "Standardized CPUE"),
      colour = "#0F766E",
      fill = "white",
      shape = 21,
      stroke = 0.7,
      size = 2.2
    ) +
    ggplot2::geom_line(
      data = subset(plot_dat, panel == "CV"),
      colour = "#B45309",
      linewidth = 1,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = subset(plot_dat, panel == "CV"),
      colour = "#B45309",
      fill = "white",
      shape = 21,
      stroke = 0.7,
      size = 1.9
    ) +
    ggplot2::labs(x = "Time", y = NULL) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey96", colour = "grey80"),
      strip.text.x = ggplot2::element_text(face = "bold"),
      strip.placement = "outside",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey92"),
      panel.border = ggplot2::element_rect(colour = "grey40"),
      axis.title.x = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(colour = "grey20"),
      axis.text.y = ggplot2::element_text(colour = "grey20")
    )

  if ("region" %in% names(index_df)) {
    p <- p +
      ggplot2::facet_grid(panel ~ region, scales = "free_y", switch = "y") +
      ggplot2::theme(
        strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold"),
        strip.background.y = ggplot2::element_rect(fill = "grey98", colour = "grey80")
      )
  } else {
    p <- p +
      ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y", strip.position = "top")
  }

  if (!is.null(time_values) || !is.null(time_positions)) {
    x_breaks <- unname(time_map)
    x_labels <- if (is.null(time_values)) time_levels else time_values
    x_labels_chr <- as.character(x_labels)
    has_month_labels <- any(grepl("-", x_labels_chr, fixed = TRUE))

    if (is.null(x_text_angle)) {
      x_text_angle <- if (has_month_labels || any(nchar(x_labels_chr) > 4L)) 90 else 45
    }

    angle_norm <- ((x_text_angle %% 360) + 360) %% 360
    if (is.null(x_text_hjust)) {
      x_text_hjust <- 1
    }
    if (is.null(x_text_vjust)) {
      x_text_vjust <- if (angle_norm %in% c(90, 270)) 0.5 else 1
    }

    p <- p +
      ggplot2::scale_x_continuous(
        breaks = x_breaks,
        labels = x_labels
      ) +
      ggplot2::labs(x = if (has_month_labels) "Year-Month" else "Year") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = x_text_angle,
          hjust = x_text_hjust,
          vjust = x_text_vjust,
          colour = "grey20"
        )
      )
  }

  p
}
