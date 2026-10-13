#' Simulate a plant-pollinator network with niche overlap.
#'
#' Plants and pollinators have positions in a shared trait space (`d`
#' standard normal traits) and a niche breadth \eqn{\sigma}, log-normal
#' within each guild. A plant and a pollinator interact in proportion to the
#' overlap of their Gaussian niches, which decreases with their trait
#' distance \eqn{d}: \eqn{\exp(-d^2 / 2(\sigma_p^2 + \sigma_a^2))}. A share
#' `lambda` of the preference is uniform instead, so `lambda = 0` is pure
#' trait matching and `lambda = 1` a random network with the same expected
#' number of links. Within each guild, log breadth is correlated with the
#' distance to the guild's centroid with coefficient `rho`: peripheral
#' species are generalists if `rho > 0`, specialists if `rho < 0`. Links are
#' Bernoulli draws whose probabilities are scaled so that plants have
#' `k_plant` pollinators on average. A species left without a partner gets
#' its most likely one.
#'
#' @param n_plants,n_pollinators number of plants and pollinators.
#' @param lambda share of interactions independent of traits, in [0, 1].
#' @param rho correlation between log breadth and distance to the centroid,
#'   in [-1, 1].
#' @param k_plant expected number of pollinators per plant.
#' @param breadth_median,breadth_sd median and standard deviation (log
#'   scale) of niche breadth.
#' @param d number of traits.
#'
#' @return list with `links` (logical matrix, plants as rows and pollinators
#'   as columns), `plant_traits`, `pollinator_traits` (matrices, one row per
#'   species), `plant_breadth` and `pollinator_breadth`.
#' @export
simulate_network <- function(n_plants, n_pollinators, lambda, rho = 0,
                             k_plant = 2, breadth_median = 0.45,
                             breadth_sd = 0.8, d = 2) {
  plant_traits <- matrix(stats::rnorm(n_plants * d), n_plants)
  pollinator_traits <- matrix(stats::rnorm(n_pollinators * d), n_pollinators)
  rownames(plant_traits) <- paste0("plant", seq_len(n_plants))
  rownames(pollinator_traits) <- paste0("pollinator", seq_len(n_pollinators))
  colnames(plant_traits) <- colnames(pollinator_traits) <- paste0("trait", seq_len(d))
  draw_breadth <- function(traits) {
    periphery <- as.numeric(scale(sqrt(rowSums(sweep(traits, 2, colMeans(traits))^2))))
    noise <- stats::rnorm(nrow(traits))
    exp(log(breadth_median) + breadth_sd * (rho * periphery + sqrt(1 - rho^2) * noise))
  }
  plant_breadth <- draw_breadth(plant_traits)
  pollinator_breadth <- draw_breadth(pollinator_traits)
  dist2 <- outer(rowSums(plant_traits^2), rowSums(pollinator_traits^2), "+") -
    2 * plant_traits %*% t(pollinator_traits)
  width2 <- outer(plant_breadth^2, pollinator_breadth^2, "+")
  preference <- exp(-dist2 / (2 * width2))
  preference <- (1 - lambda) * preference / mean(preference) + lambda
  prob <- k_plant / n_pollinators * preference / mean(preference)
  prob[prob > 1] <- 1
  links <- matrix(stats::runif(length(prob)) < prob, n_plants)
  for (i in which(rowSums(links) == 0)) links[i, which.max(prob[i, ])] <- TRUE
  for (j in which(colSums(links) == 0)) links[which.max(prob[, j]), j] <- TRUE
  dimnames(links) <- list(rownames(plant_traits), rownames(pollinator_traits))
  list(
    links = links, plant_traits = plant_traits,
    pollinator_traits = pollinator_traits,
    plant_breadth = plant_breadth, pollinator_breadth = pollinator_breadth
  )
}

#' Robustness and FD areas for one pollinator removal order, fast version.
#'
#' Same quantities as `compute_atn_area()` and `compute_fd_area()`, computed
#' from the step at which each plant dies (the position, in the removal
#' order, of its last pollinator), so that large networks stay tractable.
#'
#' @param links logical interaction matrix, plants as rows and pollinators
#'   as columns.
#' @param membership edge-by-plant matrix of the FD dendrogram, from
#'   `get_edge_membership()`, with columns in the order of `links` rows.
#' @param edge_length branch lengths of the dendrogram.
#' @param sequence column indices of `links`, in removal order.
#'
#' @return named vector: `r` (robustness) and `r_fd` (FD area).
#' @export
compute_removal_curves <- function(links, membership, edge_length, sequence) {
  n <- ncol(links)
  n_plants <- nrow(links)
  position <- match(seq_len(n), sequence)
  death <- apply(links, 1, function(row) max(position[row]))
  # FD after each successive plant death
  membership <- membership * 1
  alive_below <- rowSums(membership) -
    t(apply(membership[, order(death), drop = FALSE], 1, cumsum))
  alive_below <- cbind(rowSums(membership), alive_below)
  n_alive <- matrix(
    n_plants - 0:n_plants, nrow(alive_below), n_plants + 1, byrow = TRUE
  )
  fd <- colSums(edge_length * (alive_below > 0 & alive_below < n_alive))
  n_dead <- vapply(0:n, function(k) sum(death <= k), numeric(1))
  c(
    r = compute_curve_area((n_plants - n_dead) / n_plants),
    r_fd = compute_curve_area(fd[n_dead + 1] / fd[1])
  )
}

#' Simulate one network and describe its pollinators.
#'
#' Pollinator variables mirror those used for the Coux SEM: niche breadth,
#' originality (distance to the centroid of the pollinators' traits), the
#' mean breadth and originality of the plants visited, normalised degree,
#' and the
#' contributions \eqn{\Delta R_i} and \eqn{\Delta FD_i}, computed as in
#' `compute_removal_contributions()` with `nrep` random orders. Uniqueness
#' (distance to the nearest species of the same guild) is also returned, to
#' compare position metrics.
#'
#' @inheritParams simulate_network
#' @param network identifier of the network.
#' @param pollinator_ratio number of pollinators per plant (0.65 in Coux).
#' @param nrep number of random removal orders.
#'
#' @return data tibble with one row per pollinator.
#' @export
simulate_pollinator_data <- function(network, n_plants, lambda, rho = 0,
                                     pollinator_ratio = 0.65, k_plant = 2,
                                     nrep = 50) {
  n_pollinators <- max(2, round(pollinator_ratio * n_plants))
  sim <- simulate_network(n_plants, n_pollinators, lambda, rho, k_plant)
  links <- sim$links
  tree <- ape::as.phylo(stats::hclust(stats::dist(sim$plant_traits), "average"))
  membership <- get_edge_membership(tree)[, rownames(links), drop = FALSE]
  orders <- replicate(nrep, sample(n_pollinators), simplify = FALSE)
  first <- t(vapply(seq_len(n_pollinators), function(i) {
    rowMeans(vapply(orders, function(o) {
      compute_removal_curves(links, membership, tree$edge.length, c(i, setdiff(o, i)))
    }, numeric(2)))
  }, numeric(2)))
  originality <- function(traits) {
    sqrt(rowSums(sweep(traits, 2, colMeans(traits))^2))
  }
  uniqueness <- function(traits) {
    dist <- as.matrix(stats::dist(traits))
    diag(dist) <- Inf
    apply(dist, 1, min)
  }
  plant_originality <- originality(sim$plant_traits)
  plant_uniqueness <- uniqueness(sim$plant_traits)
  partner_mean <- function(value) colSums(links * value) / colSums(links)
  tibble::tibble(
    site = network,
    species = paste(network, colnames(links), sep = "_"),
    n_plants = n_plants,
    lambda = lambda,
    rho = rho,
    breadth = log(sim$pollinator_breadth),
    originality = originality(sim$pollinator_traits),
    uniqueness = uniqueness(sim$pollinator_traits),
    degree = colSums(links),
    normalised_degree = colSums(links) / n_plants,
    partner_breadth = partner_mean(log(sim$plant_breadth)),
    partner_originality = partner_mean(plant_originality),
    partner_uniqueness = partner_mean(plant_uniqueness),
    delta_r = mean(first[, "r"]) - first[, "r"],
    delta_fd = mean(first[, "r_fd"]) - first[, "r_fd"]
  )
}

#' Simulate replicate networks for one parameter setting.
#'
#' @param setting one-row data frame with `n_plants`, `lambda`, `rho` and
#'   `n_rep`.
#' @param cores number of cores for `parallel::mclapply()`.
#' @param seed random seed; each setting gets its own stream.
#'
#' @return data tibble with one row per pollinator and network.
#' @export
simulate_setting <- function(setting, cores = 6, seed = 1) {
  # Reproducible parallel streams; restore the previous generator afterwards.
  old_kind <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old_kind[1]))
  set.seed(
    seed + 1000 * setting$n_plants + round(100 * setting$lambda) +
      round(10 * (setting$rho + 1))
  )
  parallel::mclapply(seq_len(setting$n_rep), function(rep) {
    network <- sprintf(
      "p%d_l%.2f_rho%+.1f_r%02d", setting$n_plants, setting$lambda,
      setting$rho, rep
    )
    simulate_pollinator_data(
      network, setting$n_plants, setting$lambda, setting$rho
    ) |>
      dplyr::mutate(rep = rep)
  }, mc.cores = cores, mc.set.seed = TRUE) |>
    dplyr::bind_rows()
}

#' Fit the grid SEM to one batch of simulated networks.
#'
#' Same structure as for Coux (`fit_grid_sem()`), with breadths on the log
#' scale. Each simulated species occurs in one network only, so there is a
#' site random effect but no species effect.
#'
#' @param data pollinators of the networks in the batch, from
#'   `simulate_setting()`.
#'
#' @return list with the `piecewiseSEM::psem` object `sem` and the data it
#'   was fitted to, `data`.
#' @export
fit_simulated_sem <- function(data) {
  data <- data |>
    dplyr::mutate(
      breadth_z = as.numeric(scale(breadth)),
      originality_z = as.numeric(scale(originality)),
      delta_r = delta_r / stats::sd(delta_r),
      delta_fd = delta_fd / stats::sd(delta_fd)
    ) |>
    as.data.frame()
  list(sem = fit_grid_sem(data, species_effect = FALSE), data = data)
}

#' Fit the grid SEM to simulated networks, by batch.
#'
#' Fits `fit_simulated_sem()` separately to batches of `batch_size`
#' networks per parameter setting, like the Coux sites.
#'
#' @param sim_data from `simulate_setting()`.
#' @param batch_size number of networks per SEM.
#'
#' @return data tibble with one row per SEM and path: setting, `batch`,
#'   `path`, standardised estimate `std`, `p`, Fisher's C p-value
#'   `fisher_p`, `missing`, the omitted paths flagged by the d-separation
#'   tests (p < 0.05) separated by "; ", `breadth_position`, the
#'   correlation between pollinator breadth and originality, and `error`,
#'   the error message of a batch that failed (other columns then `NA`).
#' @export
fit_simulated_sems <- function(sim_data, batch_size = 15) {
  sim_data |>
    dplyr::mutate(batch = (rep - 1) %/% batch_size + 1) |>
    dplyr::group_by(n_plants, lambda, rho, batch) |>
    dplyr::group_modify(.keep = TRUE, function(data, key) {
      fit <- fit_simulated_sem(data)
      summary <- tryCatch(
        summarise_contribution_sem(fit$sem, fit$data),
        error = function(e) e
      )
      # A batch whose models fail to converge is recorded, not fitted.
      if (inherits(summary, "error")) {
        return(tibble::tibble(error = conditionMessage(summary)))
      }
      coefs <- summary$coefficients
      tibble::tibble(
        path = paste(coefs$Predictor, "->", coefs$Response),
        std = coefs$Std.Estimate,
        p = coefs$P.Value,
        fisher_p = summary$Cstat$P.Value,
        missing = paste(
          sub(" \\+ \\.\\.\\.$", "", summary$dTable$Independ.Claim[summary$dTable$P.Value < 0.05]),
          collapse = "; "
        ),
        breadth_position = stats::cor(fit$data$breadth_z, fit$data$originality_z)
      )
    }) |>
    dplyr::ungroup()
}

#' Compare position metrics in simulated networks.
#'
#' For each network: how well a pollinator's position predicts the same
#' metric of the plants it visits (matching), and how well that plant
#' metric explains \eqn{\Delta FD_i} beyond \eqn{\Delta R_i} (Spearman
#' correlation with the residual), for originality and uniqueness.
#'
#' @param sim_data from `simulate_setting()`.
#'
#' @return data tibble with one row per network and metric.
#' @export
compare_position_metrics <- function(sim_data) {
  sim_data |>
    dplyr::group_by(n_plants, lambda, rho, site) |>
    dplyr::filter(stats::sd(delta_r) > 0) |>
    dplyr::summarise(
      residual = list(stats::resid(stats::lm(delta_fd ~ delta_r))),
      matching_originality = stats::cor(originality, partner_originality, method = "spearman"),
      matching_uniqueness = stats::cor(uniqueness, partner_uniqueness, method = "spearman"),
      fd_originality = stats::cor(residual[[1]], partner_originality, method = "spearman"),
      fd_uniqueness = stats::cor(residual[[1]], partner_uniqueness, method = "spearman"),
      .groups = "drop"
    ) |>
    dplyr::select(-residual) |>
    tidyr::pivot_longer(
      c(dplyr::starts_with("matching"), dplyr::starts_with("fd")),
      names_to = c(".value", "metric"), names_sep = "_"
    )
}

#' Plot a simulated network in trait space.
#'
#' Plants and pollinators in the shared trait space (first two traits),
#' links as segments, plant size proportional to plant breadth. Three
#' pollinators, with the narrowest, the median and the widest niche, are
#' highlighted with their links and a circle of radius \eqn{\sigma}, their
#' niche breadth.
#'
#' @param sim from `simulate_network()`.
#' @param plant_color,pollinator_color colours of the two guilds.
#'
#' @return a ggplot.
#' @export
plot_simulated_network <- function(sim, plant_color = "#eb6834",
                                   pollinator_color = "#2a78d6") {
  plants <- tibble::as_tibble(sim$plant_traits[, 1:2], rownames = "plant") |>
    dplyr::mutate(breadth = sim$plant_breadth)
  pollinators <- tibble::as_tibble(sim$pollinator_traits[, 1:2], rownames = "pollinator") |>
    dplyr::mutate(breadth = sim$pollinator_breadth)
  ranks <- order(sim$pollinator_breadth)
  focal <- tibble::tibble(
    pollinator = rownames(sim$pollinator_traits)[ranks[c(1, ceiling(length(ranks) / 2), length(ranks))]],
    niche = factor(c("narrow", "median", "wide"), c("narrow", "median", "wide"))
  )
  links <- which(sim$links, arr.ind = TRUE)
  segments <- tibble::tibble(
    pollinator = rownames(sim$pollinator_traits)[links[, 2]],
    x = sim$plant_traits[links[, 1], 1], y = sim$plant_traits[links[, 1], 2],
    xend = sim$pollinator_traits[links[, 2], 1],
    yend = sim$pollinator_traits[links[, 2], 2]
  ) |>
    dplyr::left_join(focal, by = "pollinator")
  angle <- seq(0, 2 * pi, length.out = 90)
  circles <- pollinators |>
    dplyr::inner_join(focal, by = "pollinator") |>
    dplyr::reframe(
      niche = niche,
      x = trait1 + breadth * cos(angle), y = trait2 + breadth * sin(angle),
      .by = pollinator
    )
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = dplyr::filter(segments, is.na(niche)),
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = "grey85", linewidth = 0.3
    ) +
    ggplot2::geom_path(
      data = circles, ggplot2::aes(x, y, group = pollinator, linetype = niche),
      color = pollinator_color
    ) +
    ggplot2::geom_segment(
      data = dplyr::filter(segments, !is.na(niche)),
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linetype = niche),
      color = pollinator_color, linewidth = 0.6
    ) +
    ggplot2::geom_point(
      data = plants, ggplot2::aes(trait1, trait2, size = breadth),
      color = plant_color, shape = 15
    ) +
    ggplot2::scale_size_area(max_size = 5, name = "Plant breadth") +
    ggplot2::geom_point(
      data = pollinators, ggplot2::aes(trait1, trait2), color = pollinator_color,
      size = 2.2, alpha = 0.5
    ) +
    ggplot2::geom_point(
      data = dplyr::inner_join(pollinators, focal, by = "pollinator"),
      ggplot2::aes(trait1, trait2), color = pollinator_color, size = 3
    ) +
    ggplot2::scale_linetype_manual(
      values = c(narrow = "dotted", median = "dashed", wide = "solid"),
      name = "Pollinator niche"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Trait 1", y = "Trait 2")
}

#' Plot simulated SEM paths against trait matching.
#'
#' One panel per path of the grid SEM, grouped by axis: standardised
#' coefficient against the share of interactions independent of traits
#' (`lambda`), coloured by number of plants (at `rho = 0`) or by the
#' breadth-position correlation `rho` (averaged over richness), with the
#' Coux estimate as a dashed line.
#'
#' @param sim_sems from `fit_simulated_sems()`.
#' @param coux_coefs coefficient table of the Coux SEM, `$coefficients` of
#'   `summarise_contribution_sem()`.
#' @param colour_by `"n_plants"` or `"rho"`.
#'
#' @return a ggplot.
#' @export
plot_simulated_sems <- function(sim_sems, coux_coefs,
                                colour_by = c("n_plants", "rho")) {
  colour_by <- match.arg(colour_by)
  labels <- c(
    "breadth_z -> normalised_degree" = "Down: breadth → degree",
    "originality_z -> normalised_degree" = "Down: position → degree",
    "partner_breadth -> delta_r" = "Down: plant breadth → ΔR",
    "partner_originality -> delta_r" = "Down: plant position → ΔR",
    "breadth_z -> partner_breadth" = "Across: breadth → plant breadth",
    "originality_z -> partner_originality" = "Across: position → plant position",
    "normalised_degree -> delta_r" = "Across: degree → ΔR",
    "delta_r -> delta_fd" = "Outcome: ΔR → ΔFD",
    "partner_originality -> delta_fd" = "Outcome: plant position → ΔFD"
  )
  coux <- tibble::tibble(
    path = paste(coux_coefs$Predictor, "->", coux_coefs$Response),
    std = coux_coefs$Std.Estimate
  ) |>
    dplyr::filter(path %in% names(labels)) |>
    dplyr::mutate(path = factor(labels[path], labels))
  sims <- sim_sems |>
    dplyr::filter(!is.na(path), path %in% names(labels)) |>
    dplyr::mutate(path = factor(labels[path], labels))
  if (colour_by == "n_plants") {
    sims <- sims |>
      dplyr::filter(rho == 0) |>
      dplyr::mutate(group = factor(n_plants))
    colour_name <- "Plants"
  } else {
    sims <- sims |>
      dplyr::group_by(path, lambda, rho, batch) |>
      dplyr::summarise(std = mean(std), .groups = "drop") |>
      dplyr::mutate(group = factor(rho))
    colour_name <- "ρ (breadth–position)"
  }
  ggplot2::ggplot(sims, ggplot2::aes(lambda, std, color = group)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey80") +
    ggplot2::geom_hline(
      data = coux, ggplot2::aes(yintercept = std),
      linetype = "dashed", color = "grey30"
    ) +
    ggplot2::geom_point(alpha = 0.5, size = 1.2) +
    ggplot2::stat_summary(fun = mean, geom = "line") +
    ggplot2::facet_wrap(~path, scales = "free_y", ncol = 4) +
    ggplot2::scale_color_viridis_d(end = 0.85, name = colour_name) +
    ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1)) +
    ggplot2::labs(
      x = "Share of interactions independent of traits (λ)",
      y = "Standardised coefficient"
    )
}
