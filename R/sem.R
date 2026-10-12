#' Assemble the data for the contribution SEM.
#'
#' One row per pollinator and site, with its trait-space position, network
#' role, the uniqueness of the plants it visits, and its contributions to
#' robustness and plant FD. Pollinators without traits are dropped.
#' Contributions are multiplied by 10 to bring them to the scale of the
#' other variables, which helps the optimiser.
#'
#' @param trait_network data from `prepare_trait_network_data()`.
#' @param robustness_contributions from `compute_robustness_contributions()`.
#' @param fd_contributions from `compute_fd_contributions()`.
#' @param plant_support from `compute_plant_support()`.
#'
#' @return data frame, as `piecewiseSEM` requires.
#' @export
prepare_sem_data <- function(trait_network, robustness_contributions,
                             fd_contributions, plant_support) {
  by <- c("site", "species")
  trait_network |>
    dplyr::left_join(
      dplyr::select(robustness_contributions, site, species, delta_r_first),
      by = by
    ) |>
    dplyr::left_join(
      dplyr::select(fd_contributions, site, species, delta_fd_first),
      by = by
    ) |>
    dplyr::left_join(
      dplyr::select(plant_support, site, species, partner_uniqueness),
      by = by
    ) |>
    dplyr::filter(!is.na(originality_z)) |>
    dplyr::mutate(delta_r = 10 * delta_r_first, delta_fd = 10 * delta_fd_first) |>
    as.data.frame()
}

#' Fit the piecewise SEM of pollinator contributions.
#'
#' Two paths from pollinator originality to its contribution to plant FD:
#' - originality -> normalised degree -> \eqn{\Delta R} -> \eqn{\Delta FD};
#' - originality -> uniqueness of the plants visited -> \eqn{\Delta FD}.
#'
#' Given \eqn{\Delta R}, the plant uniqueness path acts on the part of
#' \eqn{\Delta FD} that \eqn{\Delta R} does not explain. Each equation is a
#' `glmmTMB` mixed model: binomial for normalised degree (weighted by the
#' number of plants), Gaussian otherwise. Random effects of site and species
#' are kept unless their variance is zero: plant uniqueness has no species
#' effect, and \eqn{\Delta FD} no site effect.
#'
#' @param data from `prepare_sem_data()`.
#'
#' @return a `piecewiseSEM::psem` object. `glmmTMB` must be attached
#'   (`piecewiseSEM` reads the model calls).
#' @export
fit_contribution_sem <- function(data) {
  piecewiseSEM::psem(
    glmmTMB(
      normalised_degree ~ originality_z + (1 | site) + (1 | species),
      weights = n_plants, family = "binomial", data = data
    ),
    glmmTMB(
      partner_uniqueness ~ originality_z + (1 | site),
      data = data
    ),
    glmmTMB(
      delta_r ~ normalised_degree + (1 | site) + (1 | species),
      data = data
    ),
    glmmTMB(
      delta_fd ~ delta_r + partner_uniqueness + (1 | species),
      data = data
    ),
    data = data
  )
}

#' Plot the contribution SEM as a path diagram.
#'
#' Arrows are labelled with standardised coefficients, their width scales
#' with the coefficient's absolute value and their colour with its sign.
#' Dashed arrows have p > 0.05.
#'
#' @param sem a `piecewiseSEM::psem` object, from `fit_contribution_sem()`.
#' @param positive,negative colours of positive and negative paths.
#'
#' @return a ggplot.
#' @export
plot_sem_paths <- function(sem, positive = "#2166ac", negative = "#b2182b") {
  nodes <- tibble::tribble(
    ~name, ~label, ~x, ~y, ~half_w, ~half_h,
    "originality_z", "Originality", 0, 0.8, 0.27, 0.09,
    "normalised_degree", "Normalised\ndegree", 1.2, 1.6, 0.29, 0.13,
    "delta_r", "\u0394R", 2.4, 1.6, 0.12, 0.09,
    "partner_uniqueness", "Uniqueness of\nplants visited", 1.8, 0, 0.35, 0.13,
    "delta_fd", "\u0394FD", 3.4, 0.8, 0.15, 0.09
  )
  gap <- 0.04
  # Distance from a node's centre to its box edge, along direction (dx, dy).
  to_edge <- function(half_w, half_h, dx, dy) {
    pmin(half_w / abs(dx), half_h / abs(dy)) * sqrt(dx^2 + dy^2) + gap
  }
  coefs <- piecewiseSEM::coefs(sem)
  paths <- tibble::tibble(
    from = coefs$Predictor, to = coefs$Response,
    std = coefs$Std.Estimate, p = coefs$P.Value
  ) |>
    dplyr::left_join(
      dplyr::select(
        nodes, from = name, x0 = x, y0 = y, w0 = half_w, h0 = half_h
      ),
      by = "from"
    ) |>
    dplyr::left_join(
      dplyr::select(
        nodes, to = name, x1 = x, y1 = y, w1 = half_w, h1 = half_h
      ),
      by = "to"
    ) |>
    dplyr::mutate(
      dx = x1 - x0,
      dy = y1 - y0,
      len = sqrt(dx^2 + dy^2),
      xs = x0 + to_edge(w0, h0, dx, dy) * dx / len,
      ys = y0 + to_edge(w0, h0, dx, dy) * dy / len,
      xe = x1 - to_edge(w1, h1, dx, dy) * dx / len,
      ye = y1 - to_edge(w1, h1, dx, dy) * dy / len,
      sign = ifelse(std > 0, "positive", "negative"),
      significant = p <= 0.05,
      label = sprintf("%.2f", std)
    )
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = paths,
      ggplot2::aes(
        x = xs, y = ys, xend = xe, yend = ye, color = sign,
        linewidth = abs(std), linetype = significant
      ),
      arrow = grid::arrow(length = grid::unit(2.5, "mm"), type = "closed")
    ) +
    ggplot2::geom_label(
      data = paths,
      ggplot2::aes(x = (xs + xe) / 2, y = (ys + ye) / 2, label = label,
                   color = sign),
      size = 3.2, linewidth = 0, fill = "white", show.legend = FALSE
    ) +
    ggplot2::geom_label(
      data = nodes, ggplot2::aes(x = x, y = y, label = label),
      size = 3.5, label.padding = grid::unit(2.5, "mm"), fill = "grey95",
      linewidth = 0.3, lineheight = 0.9
    ) +
    ggplot2::scale_color_manual(
      values = c(positive = positive, negative = negative), name = NULL
    ) +
    ggplot2::scale_linewidth(range = c(0.4, 2), guide = "none") +
    ggplot2::scale_linetype_manual(
      values = c(`TRUE` = "solid", `FALSE` = "dashed"),
      labels = c(`TRUE` = "p ≤ 0.05", `FALSE` = "p > 0.05"), name = NULL
    ) +
    ggplot2::coord_fixed(xlim = c(-0.4, 3.8), ylim = c(-0.3, 1.9)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom")
}
