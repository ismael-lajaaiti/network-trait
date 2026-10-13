#' Assemble the data for the contribution SEM.
#'
#' One row per pollinator and site, organised along two axes, trait space
#' vs network and pollinator vs plant:
#' - pollinator traits: breadth and position (originality);
#' - plant traits: mean breadth and originality of the plants visited;
#' - pollinator network role: normalised degree;
#' - plant network role: \eqn{\Delta R}, the plants that depend on the
#'   pollinator;
#' - outcome: \eqn{\Delta FD}.
#'
#' Pollinators without traits or without a breadth value are dropped.
#' Contributions are multiplied by 10 to bring them to the scale of the
#' other variables, which helps the optimiser.
#'
#' @param trait_network data from `prepare_trait_network_data()`.
#' @param robustness_contributions from `compute_robustness_contributions()`.
#' @param fd_contributions from `compute_fd_contributions()`.
#' @param plant_support from `compute_plant_support()`.
#' @param breadth pollinator breadth, one row per pollinator: `species` and
#'   `breadth`, e.g. from `get_social_breadth()`.
#' @param plant_breadth plant breadth, one row per plant: `species` and
#'   `breadth`, e.g. from `get_floral_breadth()`.
#' @param web_list named list of interaction matrices, plants as rows and
#'   pollinators as columns.
#'
#' @return data frame, as `piecewiseSEM` requires.
#' @export
prepare_sem_data <- function(trait_network, robustness_contributions,
                             fd_contributions, plant_support, breadth,
                             plant_breadth, web_list) {
  by <- c("site", "species")
  partner_breadth <- lapply(names(web_list), function(s) {
    links <- web_list[[s]] > 0
    value <- plant_breadth$breadth[match(rownames(links), plant_breadth$species)]
    tibble::tibble(
      site = s, species = colnames(links),
      partner_breadth = colSums(links * value) / colSums(links)
    )
  }) |>
    dplyr::bind_rows()
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
      dplyr::select(
        plant_support, site, species, partner_originality, partner_uniqueness
      ),
      by = by
    ) |>
    dplyr::left_join(partner_breadth, by = by) |>
    dplyr::left_join(
      dplyr::select(breadth, species, breadth),
      by = "species"
    ) |>
    dplyr::filter(!is.na(originality_z), !is.na(breadth)) |>
    dplyr::mutate(
      breadth_z = as.numeric(scale(breadth)),
      delta_r = 10 * delta_r_first,
      delta_fd = 10 * delta_fd_first
    ) |>
    as.data.frame()
}

#' Pollinator niche breadth from sociality.
#'
#' Social pollinators maintain colonies over a long season, so they forage
#' over many flowering periods and must collect from whatever flowers, a
#' classic reason why social bees are generalists. Breadth is 1 for social
#' species and 0 for solitary ones.
#'
#' @param trait_table pollinator traits, from `prepare_pollinator_traits()`.
#'
#' @return data tibble with one row per pollinator: `species` and `breadth`.
#' @export
get_social_breadth <- function(trait_table) {
  tibble::tibble(
    species = rownames(trait_table),
    breadth = as.numeric(trait_table$soc_sol == "social")
  )
}

#' Plant niche breadth from inflorescence symmetry.
#'
#' Radially symmetric (actinomorphic) inflorescences are open to a wide
#' range of pollinators, whereas bilaterally symmetric (zygomorphic) ones
#' restrict access to pollinators of matching shape and behaviour, a classic
#' marker of specialised pollination. We use the symmetry of the
#' inflorescence, the unit pollinators land on, rather than of the single
#' flower: the heads of *Taraxacum*, *Cirsium* or *Trifolium* are made of
#' zygomorphic florets but are open, radially symmetric resources. Breadth
#' is 1 for actinomorphic inflorescences and 0 otherwise.
#'
#' @param trait_table plant traits, from `prepare_plant_traits()`.
#'
#' @return data tibble with one row per plant: `species` and `breadth`.
#' @export
get_floral_breadth <- function(trait_table) {
  tibble::tibble(
    species = rownames(trait_table),
    breadth = as.numeric(trait_table$inflorescence_symmetry == "actinomorphic")
  )
}

#' Fit the grid SEM of pollinator contributions.
#'
#' Variables sit on two axes, trait space vs network and pollinator vs
#' plant, and paths only run along an axis:
#' - down, pollinator: breadth and position -> normalised degree;
#' - across, trait space: pollinator position -> position of the plants
#'   visited (trait matching), pollinator breadth -> breadth of the plants
#'   visited (specialists can only reach broad plants);
#' - across, network: normalised degree -> \eqn{\Delta R} (dependence);
#' - down, plant: breadth and position of the plants visited ->
#'   \eqn{\Delta R} (narrow, peripheral plants have few pollinators);
#' - outcome: \eqn{\Delta R} and position of the plants visited ->
#'   \eqn{\Delta FD} (how many plants depend on the pollinator, and which).
#'
#' Breadth and position are correlated within each guild (correlated errors
#' for the plants visited). Each equation is a `glmmTMB` mixed model:
#' binomial for normalised degree (weighted by the number of plants),
#' Gaussian otherwise. \eqn{\Delta FD} sums to zero within a site, so its
#' equation has no site effect.
#'
#' @param data from `prepare_sem_data()`, or simulated pollinators.
#' @param species_effect add a species random effect (`FALSE` when species
#'   occur in one site only, as in simulations).
#'
#' @return a `piecewiseSEM::psem` object. `glmmTMB` must be attached
#'   (`piecewiseSEM` reads the model calls).
#' @export
fit_grid_sem <- function(data, species_effect = TRUE) {
  re <- if (species_effect) "+ (1 | site) + (1 | species)" else "+ (1 | site)"
  f <- function(formula) stats::as.formula(paste(formula, re))
  degree_f <- f("normalised_degree ~ breadth_z + originality_z")
  originality_f <- f("partner_originality ~ originality_z")
  breadth_f <- f("partner_breadth ~ breadth_z")
  delta_r_f <- f("delta_r ~ normalised_degree + partner_breadth + partner_originality")
  delta_fd <- if (species_effect) {
    glmmTMB(delta_fd ~ delta_r + partner_originality + (1 | species), data = data)
  } else {
    stats::lm(delta_fd ~ delta_r + partner_originality, data = data)
  }
  piecewiseSEM::psem(
    glmmTMB(degree_f, weights = n_plants, family = "binomial", data = data),
    glmmTMB(originality_f, data = data),
    glmmTMB(breadth_f, data = data),
    glmmTMB(delta_r_f, data = data),
    delta_fd,
    piecewiseSEM::`%~~%`(partner_breadth, partner_originality),
    data = data
  )
}

#' Fit the grid SEM to the Coux data.
#'
#' @param data from `prepare_sem_data()`.
#'
#' @return a `piecewiseSEM::psem` object, see `fit_grid_sem()`.
#' @export
fit_contribution_sem <- function(data) {
  fit_grid_sem(data, species_effect = TRUE)
}

#' Summarise the contribution SEM.
#'
#' `piecewiseSEM` refits the models for its d-separation tests and then
#' looks up the binomial weights (`n_plants`) in the global environment, not
#' in the data. We set them there for the duration of the call and restore
#' the previous state afterwards.
#'
#' @param sem from `fit_contribution_sem()`.
#' @param data from `prepare_sem_data()`, the data `sem` was fitted to.
#'
#' @return a `summary.psem` object.
#' @export
summarise_contribution_sem <- function(sem, data) {
  env <- globalenv()
  had <- exists("n_plants", envir = env, inherits = FALSE)
  if (had) old <- get("n_plants", envir = env)
  assign("n_plants", data$n_plants, envir = env)
  on.exit(if (had) assign("n_plants", old, envir = env) else rm("n_plants", envir = env))
  summary(sem, .progressBar = FALSE)
}

#' Plot the grid SEM as a path diagram.
#'
#' Variables are laid out on the two axes of the SEM: trait space (top) vs
#' network (bottom), pollinator (left) vs plant (right), with
#' \eqn{\Delta FD} as the outcome on the right. Arrows are labelled with
#' standardised coefficients, their width scales with the coefficient's
#' absolute value and their colour with its sign; dashed arrows have
#' p > 0.05. Grey double-headed segments show the correlation between
#' breadth and position within each guild.
#'
#' @param sem a `piecewiseSEM::psem` object, from `fit_grid_sem()`.
#' @param data the data `sem` was fitted to.
#' @param breadth_label,plant_breadth_label labels of the breadth nodes.
#' @param positive,negative colours of positive and negative paths.
#'
#' @return a ggplot.
#' @export
plot_sem_paths <- function(sem, data,
                           breadth_label = "Pollinator\nbreadth",
                           plant_breadth_label = "Breadth of\nplants visited",
                           positive = "#2166ac", negative = "#b2182b") {
  nodes <- tibble::tribble(
    ~name, ~label, ~x, ~y, ~half_w, ~half_h,
    "breadth_z", breadth_label, -0.35, 1.95, 0.36, 0.14,
    "originality_z", "Pollinator\noriginality", 0.55, 1.4, 0.34, 0.14,
    "partner_breadth", plant_breadth_label, 2.45, 1.95, 0.42, 0.14,
    "partner_originality", "Originality of\nplants visited", 3.35, 1.4, 0.42, 0.14,
    "normalised_degree", "Normalised\ndegree", 0.1, 0, 0.34, 0.14,
    "delta_r", "ΔR", 2.9, 0, 0.13, 0.1,
    "delta_fd", "ΔFD", 4.85, 0.7, 0.17, 0.1
  )
  gap <- 0.04
  # Distance from a node's centre to its box edge, along direction (dx, dy).
  to_edge <- function(half_w, half_h, dx, dy) {
    pmin(half_w / abs(dx), half_h / abs(dy)) * sqrt(dx^2 + dy^2) + gap
  }
  geometry <- function(edges) {
    edges |>
      dplyr::left_join(
        dplyr::select(nodes, from = name, x0 = x, y0 = y, w0 = half_w, h0 = half_h),
        by = "from"
      ) |>
      dplyr::left_join(
        dplyr::select(nodes, to = name, x1 = x, y1 = y, w1 = half_w, h1 = half_h),
        by = "to"
      ) |>
      dplyr::mutate(
        dx = x1 - x0, dy = y1 - y0, len = sqrt(dx^2 + dy^2),
        xs = x0 + to_edge(w0, h0, dx, dy) * dx / len,
        ys = y0 + to_edge(w0, h0, dx, dy) * dy / len,
        xe = x1 - to_edge(w1, h1, dx, dy) * dx / len,
        ye = y1 - to_edge(w1, h1, dx, dy) * dy / len
      )
  }
  coefs <- piecewiseSEM::coefs(sem)
  directed <- !startsWith(coefs$Response, "~~")
  paths <- tibble::tibble(
    from = coefs$Predictor[directed], to = coefs$Response[directed],
    std = coefs$Std.Estimate[directed], p = coefs$P.Value[directed]
  ) |>
    geometry() |>
    dplyr::mutate(
      sign = ifelse(std > 0, "positive", "negative"),
      significant = p <= 0.05,
      label = sprintf("%.2f", std)
    )
  correlations <- tibble::tibble(
    from = c("breadth_z", "partner_breadth"),
    to = c("originality_z", "partner_originality"),
    r = c(
      stats::cor(data$breadth_z, data$originality_z),
      coefs$Std.Estimate[!directed]
    )
  ) |>
    geometry() |>
    dplyr::mutate(label = sprintf("r = %.2f", r))
  cells <- tibble::tibble(
    xmin = c(-0.85, 1.9, -0.85, 1.9), xmax = c(1.05, 3.9, 1.05, 3.9),
    ymin = c(1.05, 1.05, -0.45, -0.45), ymax = c(2.45, 2.45, 0.45, 0.45)
  )
  ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = cells,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "grey97", color = "grey85"
    ) +
    ggplot2::annotate(
      "text", x = c(0.1, 2.9), y = 2.65, label = c("Pollinator", "Plant"),
      fontface = "bold", size = 4
    ) +
    ggplot2::annotate(
      "text", x = -1.05, y = c(1.75, 0), label = c("Trait space", "Network"),
      fontface = "bold", size = 4, angle = 90
    ) +
    ggplot2::geom_segment(
      data = correlations,
      ggplot2::aes(x = xs, y = ys, xend = xe, yend = ye),
      color = "grey60", linetype = "dotted",
      arrow = grid::arrow(length = grid::unit(1.8, "mm"), ends = "both", type = "closed")
    ) +
    ggplot2::geom_text(
      data = correlations,
      ggplot2::aes(x = (xs + xe) / 2 + 0.28, y = (ys + ye) / 2 + 0.05, label = label),
      color = "grey40", size = 2.8
    ) +
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
      ggplot2::aes(x = (xs + xe) / 2, y = (ys + ye) / 2, label = label, color = sign),
      size = 3, linewidth = 0, fill = "white", show.legend = FALSE
    ) +
    ggplot2::geom_label(
      data = nodes, ggplot2::aes(x = x, y = y, label = label),
      size = 3.3, label.padding = grid::unit(2.2, "mm"), fill = "white",
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
    ggplot2::coord_fixed(xlim = c(-1.2, 5.2), ylim = c(-0.55, 2.8)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom")
}
