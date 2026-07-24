#' Read and format ata from Coux et al. 2016.
#'
#' @param dir path where the data is stored.
#'
#' @return list of data frames.
#' @export
read_coux_data <- function(dir) {
  data <- list() # Init.
  data$interaction <- read.csv(here::here(dir, "interactions.csv")) |>
    dplyr::select(-X) |>
    dplyr::rename(
      site = Site, pollinator = Pol_sp, plant = Plant_sp, links = Links
    )
  data$plant_trait <- read.csv(here::here(dir, "plant_traits.csv")) |>
    dplyr::rename(plant = X)
  data$plant_abundance <- read.csv(
    here::here(dir, "plant_abundances_bin.csv")
  ) |>
    dplyr::rename(site = X)
  data$pollinator_trait <- read.csv(
    here::here(dir, "pollinator_traits.csv")
  ) |>
    dplyr::rename(pollinator = X)
  data$pollinator_abundance <- read.csv(
    here::here(dir, "pollinator_abundances.csv")
  ) |> dplyr::rename(site = Site)
  data
}

#' Clean species names in Coux et al. 2016 data.
#'
#' Lowercases plant and pollinator names throughout (pollinator names are
#' capitalised in the raw data, plant names are inconsistently so), and
#' reconciles `alnus_glutinosa`/`alnus_serrulate` into a single `alnus_sp`.
#' `alnus_glutinosa` appears only at site `blahamlak`, with no matching
#' trait/abundance row; `alnus_serrulate` appears at `crobirlee` and
#' `crostelee` and joins fine there. No site uses both names, so we treat
#' all three as the same plant recorded under two names rather than
#' verifying species identity from the raw data.
#'
#' @param data list of data frames, as returned by `read_coux_data()`.
#'
#' @return the same list, with cleaned species names.
#' @export
clean_coux_data <- function(data) {
  clean_name <- function(x) {
    x <- tolower(x)
    dplyr::recode_values(
      x,
      from = c("alnus_glutinosa", "alnus_serrulate"),
      to = c("alnus_sp", "alnus_sp"),
      default = x
    )
  }

  data$interaction <- data$interaction |>
    dplyr::mutate(
      plant = clean_name(plant), pollinator = clean_name(pollinator)
    )
  data$plant_trait <- data$plant_trait |>
    dplyr::mutate(plant = clean_name(plant))
  data$pollinator_trait <- data$pollinator_trait |>
    dplyr::mutate(pollinator = clean_name(pollinator))
  data$plant_abundance <- data$plant_abundance |>
    dplyr::rename_with(clean_name, .cols = -site)
  data$pollinator_abundance <- data$pollinator_abundance |>
    dplyr::rename_with(clean_name, .cols = -site)

  data
}

#' Convert interaction table to interaction matrices.
#'
#' Pollinators as columns. Plants as rows. Grouped by site.
#'
#' @param interaction_table.
#'
#' @return named list of interaction matrices.
#' @export
get_interaction_matrix <- function(interaction_table) {
  web_list <- interaction_table |>
    dplyr::group_by(site) |>
    dplyr::group_map(~ {
      web <- .x |>
        tidyr::pivot_wider(names_from = pollinator, values_from = links)
      mat <- as.matrix(web |> select(-plant))
      rownames(mat) <- web$plant
      mat
    })
  names(web_list) <- interaction_table |>
    dplyr::pull(site) |>
    unique()
  web_list
}

#' Compute per-species network metrics
#'
#' Kept to four metrics, following Coux et al. (2016)'s supplementary
#' analysis: `degree` (number of partners) and `normalised_degree` (degree
#' scaled by the number of possible partners at that site, partially
#' correcting for the site-richness confound noted elsewhere in this
#' report) for structural embeddedness; Blüthgen's `d_prime`
#' (null-model-corrected specialisation, deviation of a species'
#' interaction distribution from random given partner interaction
#' frequency) and `hs` (raw Shannon entropy of a species' interaction
#' distribution, uncorrected for partner frequency -- `bipartite`'s
#' "partner diversity") for specialisation/generality. Other
#' `bipartite::specieslevel()` indices are dropped: centrality measures
#' (betweenness, closeness) are unreliable on many of these networks ("too
#' few nodes" on several sites), and the other specialisation indices
#' (PDI, species specificity, resource range) are, unlike `d'`, not
#' corrected for partner frequency.
#'
#' @param web_list named list of interaction matrices, as returned by
#'   `get_interaction_matrix()`.
#'
#' @return data tibble of `degree`, `normalised_degree`, `d_prime` and
#'   `hs` per species and site.
#' @export
compute_network_metrics <- function(web_list) {
  lapply(web_list, function(web) {
    metrics <- web |> bipartite::specieslevel(
      index = c("degree", "normalised degree", "d", "partner diversity")
    )
    rbind(
      metrics[[1]] |>
        tibble::rownames_to_column("species") |>
        dplyr::mutate(guild = "pollinator"),
      metrics[[2]] |>
        tibble::rownames_to_column("species") |>
        dplyr::mutate(guild = "plant")
    )
  }) |>
    tibble::enframe() |>
    dplyr::rename(site = name) |>
    tidyr::unnest(value) |>
    dplyr::rename(
      normalised_degree = normalised.degree, d_prime = d, hs = partner.diversity
    )
}

#' Type plant trait columns for Gower dissimilarity.
#'
#' Ordinal scores (flower count, fragrance, nectar amount) become ordered
#' factors, single/multi-category traits become factors, and the four
#' flowering-season columns are left as 0/1 integers to be passed as
#' asymmetric binary traits to `gowdis()` (co-absence from a season is not
#' evidence of similarity).
#'
#' @param plant_trait data frame as returned by `clean_coux_data()`.
#'
#' @return data frame of typed traits with `plant` as row names.
#' @export
prepare_plant_traits <- function(plant_trait) {
  plant_trait |>
    dplyr::mutate(
      growth_form = factor(growth_form),
      annual_perennial = factor(annual_perennial),
      single_inflo = factor(single_inflo),
      type_inflorescence = factor(type_inflorescence),
      pollnec_access = factor(pollnec_access),
      flowers_per_inflorescence = factor(
        flowers_per_inflorescence,
        ordered = TRUE
      ),
      flower_symmetry = factor(flower_symmetry),
      inflorescence_symmetry = factor(inflorescence_symmetry),
      flower_sex = factor(flower_sex),
      floral_fragrance = factor(floral_fragrance, ordered = TRUE),
      amount_of_nectar_per_poll_unit_per_day = factor(
        amount_of_nectar_per_poll_unit_per_day,
        ordered = TRUE
      )
    ) |>
    tibble::column_to_rownames("plant")
}

#' Type pollinator trait columns for Gower dissimilarity.
#'
#' Body measurements and foraging preferences stay numeric; the six
#' larval-diet columns are left as 0/1 integers to be passed as asymmetric
#' binary traits to `gowdis()` (co-absence from a diet type is not evidence
#' of similarity).
#'
#' @param pollinator_trait data frame as returned by `clean_coux_data()`.
#'
#' @return data frame of typed traits with `pollinator` as row names.
#' @export
prepare_pollinator_traits <- function(pollinator_trait) {
  pollinator_trait |>
    dplyr::mutate(
      soc_sol = factor(soc_sol),
      Carrying_structure = factor(Carrying_structure),
      season = factor(season),
      daily = factor(daily)
    ) |>
    tibble::column_to_rownames("pollinator")
}

#' Compute Gower dissimilarity and PCoA ordination among species traits.
#'
#' Follows Coux et al. (2016): traits are standardised (z-scores for
#' numeric traits, range-standardisation for the rest) and combined via
#' Gower's (1971) coefficient (`FD::gowdis`), then ordinated with PCoA
#' (`ape::pcoa`), with a Cailliez correction for negative eigenvalues.
#'
#' @param trait_table data frame of traits, one row per species, typed by
#'   `prepare_plant_traits()` or `prepare_pollinator_traits()`.
#' @param asym.bin optional vector of column indices to treat as
#'   asymmetric binary traits (see `FD::gowdis`).
#'
#' @return an `ape::pcoa` object.
#' @export
compute_trait_pcoa <- function(trait_table, asym.bin = NULL) {
  dist <- FD::gowdis(trait_table, asym.bin = asym.bin)
  ape::pcoa(dist, correction = "cailliez")
}

#' Fit each trait to a PCoA ordination, one trait at a time.
#'
#' Wraps `vegan::envfit()`, called separately for each trait column so
#' that a trait's missing values only drop rows for that trait, instead of
#' `envfit()`'s default of dropping any species missing *any* trait (which
#' would needlessly shrink the sample for fully-observed traits). Gives,
#' per trait, its direction (for numeric traits) or level centroids (for
#' factors) plus an R^2 and permutation p-value for how well it fits the
#' ordination.
#'
#' @param pcoa an `ape::pcoa` object, as returned by `compute_trait_pcoa()`.
#' @param trait_table data frame of traits, one row per species matching
#'   `pcoa`, typed by `prepare_plant_traits()` or
#'   `prepare_pollinator_traits()`.
#' @param permutations number of permutations used by `vegan::envfit()`.
#' @param seed random seed, for reproducible permutation p-values.
#'
#' @return named list of `vegan::envfit` objects, one per trait column.
#' @export
fit_trait_ordination <- function(pcoa, trait_table, permutations = 999,
                                 seed = 1) {
  scores <- pcoa$vectors[, 1:2]
  set.seed(seed)
  fits <- lapply(names(trait_table), function(v) {
    vegan::envfit(
      scores, trait_table[, v, drop = FALSE],
      permutations = permutations, na.rm = TRUE
    )
  })
  names(fits) <- names(trait_table)
  fits
}

#' Tidy a list of single-trait `envfit` fits into one row per trait.
#'
#' @param fits named list from `fit_trait_ordination()`.
#'
#' @return data frame with columns `trait`, `r2`, `pval`.
#' @export
trait_importance_table <- function(fits) {
  rows <- lapply(names(fits), function(v) {
    f <- fits[[v]]
    if (!is.null(f$vectors)) {
      data.frame(
        trait = v, r2 = f$vectors$r, pval = f$vectors$pvals, row.names = NULL
      )
    } else {
      data.frame(
        trait = v, r2 = f$factors$r, pval = f$factors$pvals, row.names = NULL
      )
    }
  })
  do.call(rbind, rows)
}

#' Extract biplot-ready trait directions from a list of `envfit` fits.
#'
#' Numeric traits become arrows (direction scaled by fit strength);
#' categorical traits become one centroid per level. Only traits
#' significant at `alpha` are kept, to keep the biplot readable.
#'
#' @param fits named list from `fit_trait_ordination()`.
#' @param alpha significance threshold on the permutation p-value.
#'
#' @return list with `arrows` and `centroids` data frames (`NULL` if no
#'   trait of that type is significant).
#' @export
trait_biplot_data <- function(fits, alpha = 0.05) {
  arrows <- lapply(names(fits), function(v) {
    f <- fits[[v]]
    if (is.null(f$vectors) || f$vectors$pvals > alpha) {
      return(NULL)
    }
    coords <- f$vectors$arrows * sqrt(f$vectors$r)
    data.frame(trait = v, Axis1 = coords[1], Axis2 = coords[2])
  })
  arrows <- do.call(rbind, arrows)

  centroids <- lapply(names(fits), function(v) {
    f <- fits[[v]]
    if (is.null(f$factors) || f$factors$pvals > alpha) {
      return(NULL)
    }
    cent <- as.data.frame(f$factors$centroids)
    names(cent) <- c("Axis1", "Axis2")
    cent$level <- sub(paste0("^", v), "", rownames(cent))
    cent$variable <- v
    cent
  })
  centroids <- do.call(rbind, centroids)

  list(arrows = arrows, centroids = centroids)
}

#' Plot each trait's R^2 against a PCoA ordination, ranked.
#'
#' @param fits named list from `fit_trait_ordination()`.
#'
#' @return a ggplot object.
#' @export
plot_trait_importance <- function(fits) {
  trait_importance_table(fits) |>
    dplyr::mutate(trait = forcats::fct_reorder(trait, r2)) |>
    ggplot2::ggplot(ggplot2::aes(x = trait, y = r2, fill = pval < 0.05)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = "#2a78d6", `FALSE` = "grey70"),
      labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s."), name = NULL
    ) +
    ggplot2::labs(x = NULL, y = expression(R^2 ~ "(fit to PCoA1-2)"))
}

#' Biplot of a trait PCoA with significant trait directions overlaid.
#'
#' @param pcoa an `ape::pcoa` object, as returned by `compute_trait_pcoa()`.
#' @param fits named list from `fit_trait_ordination()`, matching `pcoa`.
#' @param arrow_color color used for numeric-trait arrows and their labels.
#'
#' @return a ggplot object.
#' @export
plot_trait_biplot <- function(pcoa, fits, arrow_color) {
  pct <- round(pcoa$values$Rel_corr_eig[1:2] * 100, 1)
  sp_df <- as.data.frame(pcoa$vectors[, 1:2])
  names(sp_df) <- c("Axis1", "Axis2")

  bp <- trait_biplot_data(fits)
  if (!is.null(bp$arrows)) {
    point_extent <- max(sqrt(sp_df$Axis1^2 + sp_df$Axis2^2))
    arrow_extent <- max(sqrt(bp$arrows$Axis1^2 + bp$arrows$Axis2^2))
    mul <- 0.9 * point_extent / arrow_extent
    bp$arrows$Axis1 <- bp$arrows$Axis1 * mul
    bp$arrows$Axis2 <- bp$arrows$Axis2 * mul
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = sp_df, ggplot2::aes(Axis1, Axis2), color = "grey75", size = 1.5
    )

  if (!is.null(bp$centroids)) {
    p <- p +
      ggplot2::geom_point(
        data = bp$centroids, ggplot2::aes(Axis1, Axis2, color = variable),
        size = 2, shape = 17
      ) +
      ggrepel::geom_text_repel(
        data = bp$centroids,
        ggplot2::aes(Axis1, Axis2, label = level, color = variable),
        size = 3, show.legend = FALSE, seed = 1
      )
  }
  if (!is.null(bp$arrows)) {
    p <- p +
      ggplot2::geom_segment(
        data = bp$arrows,
        ggplot2::aes(x = 0, y = 0, xend = Axis1, yend = Axis2),
        arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
        color = arrow_color
      ) +
      ggrepel::geom_text_repel(
        data = bp$arrows, ggplot2::aes(Axis1, Axis2, label = trait),
        color = arrow_color, size = 3, fontface = "italic", seed = 1
      )
  }
  p +
    ggplot2::labs(
      x = paste0("PCoA1 (", pct[1], "%)"),
      y = paste0("PCoA2 (", pct[2], "%)"),
      color = NULL
    ) +
    ggplot2::coord_equal()
}

#' Compute species uniqueness per species and site
#'
#' @param trait data tibble of species trait
#' @param presence data tibble of species abundance or presence at given site,
#' site should be the first column of the tibble
#'
#' @return data tibble of species uniquess per site
#' @export
compute_uniqueness <- function(trait, presence) {
  distance_matrix <- data.matrix(FD::gowdis(trait))

  site_species_matrix <- data.matrix(presence[, -1]) # Remove site column.
  row.names(site_species_matrix) <- presence$site

  suppressMessages(
    apply(site_species_matrix, 1, function(present) {
      species_here <- present[present > 0 & !is.na(present)]
      funrar::uniqueness(t(as.matrix(species_here)), distance_matrix)
    })
  ) |>
    dplyr::bind_rows(.id = "site") |>
    dplyr::rename(uniqueness = Ui)
}

#' Compute species functional originality per species and site
#'
#' Functional originality is the distance of a species from the trait
#' centroid of the community it belongs to (Laliberte & Legendre 2010;
#' Buisson et al. 2013): species whose traits deviate most from their
#' community's average are the most original. Distances are computed in
#' the Cailliez-corrected PCoA space so that they are proper Euclidean
#' distances. The centroid, across present species, of each PCoA axis is
#' either a plain mean (`weighted = FALSE`) or a mean weighted by each
#' species' value in `presence` (`weighted = TRUE`), e.g. its abundance at
#' that site.
#'
#' @param pcoa an `ape::pcoa` object (with a Cailliez correction, so that
#'   `vectors.cor` is available), as returned by `compute_trait_pcoa()`.
#' @param presence data tibble of species abundance or presence at given site,
#' site should be the first column of the tibble
#' @param weighted if `TRUE`, weight the centroid by `presence` (e.g.
#'   abundance) instead of treating every present species equally.
#'
#' @return data tibble of species originality per site
#' @export
compute_originality <- function(pcoa, presence, weighted = FALSE) {
  coords <- pcoa$vectors.cor |>
    tibble::as_tibble(rownames = "species") |>
    tidyr::pivot_longer(-species, names_to = "axis", values_to = "value")

  presence_long <- presence |>
    tidyr::pivot_longer(-site, names_to = "species", values_to = "present") |>
    dplyr::filter(present > 0, !is.na(present))

  if (!weighted) {
    presence_long <- presence_long |> dplyr::mutate(present = 1)
  }

  presence_long |>
    dplyr::inner_join(
      coords,
      by = "species", relationship = "many-to-many"
    ) |>
    dplyr::group_by(site, axis) |>
    dplyr::mutate(centroid = weighted.mean(value, w = present)) |>
    dplyr::group_by(site, species) |>
    dplyr::summarise(
      originality = sqrt(sum((value - centroid)^2)), .groups = "drop"
    )
}

#' Compute interaction niche size and originality per species and site
#'
#' Following Dehling & Stouffer (2018) and Dehling et al. (2022): a
#' species' interaction niche position is the link-weighted centroid of
#' its partners' traits in the *partner* guild's Cailliez-corrected PCoA
#' space. From that position we
#' compute:
#' - niche size: the link-weighted dispersion of a species' partners
#'   around their own centroid (Laliberte & Legendre 2010-style, as in
#'   `compute_originality()`, rather than a hypervolume; a species with a
#'   single partner has a niche size of 0.
#' - niche originality: the distance of a species' niche position from
#'   the (unweighted) mean niche position across the same guild at that
#'   site.
#'
#' @param web_list named list of interaction matrices, plants as rows and
#'   pollinators as columns, as returned by `get_interaction_matrix()`.
#' @param partner_pcoa an `ape::pcoa` object (with a Cailliez correction,
#'   so that `vectors.cor` is available) for the *partner* guild's
#'   traits, e.g. `plant_trait_pcoa` to compute pollinator niches.
#' @param level `"higher"` to compute niches for pollinators (web
#'   columns) against plant trait space, or `"lower"` for plants (web
#'   rows) against pollinator trait space.
#'
#' @return data tibble of `niche_size` and `niche_originality` per
#'   species and site.
#' @export
compute_interaction_niche <- function(web_list, partner_pcoa,
                                      level = c("higher", "lower")) {
  level <- match.arg(level)

  coords <- partner_pcoa$vectors.cor |>
    tibble::as_tibble(rownames = "partner") |>
    tidyr::pivot_longer(-partner, names_to = "axis", values_to = "value")

  links <- lapply(web_list, function(web) {
    web <- if (level == "higher") t(web) else web
    as.data.frame(web) |>
      tibble::rownames_to_column("species") |>
      tidyr::pivot_longer(
        -species,
        names_to = "partner",
        values_to = "n_links"
      ) |>
      dplyr::filter(n_links > 0)
  }) |>
    dplyr::bind_rows(.id = "site")

  position <- links |>
    dplyr::inner_join(
      coords,
      by = "partner", relationship = "many-to-many"
    ) |>
    dplyr::group_by(site, species, axis) |>
    dplyr::mutate(centroid = weighted.mean(value, w = n_links)) |>
    dplyr::ungroup()

  niche_position <- position |>
    dplyr::distinct(site, species, axis, centroid)

  niche_size <- position |>
    dplyr::group_by(site, species, partner) |>
    dplyr::summarise(
      n_links = dplyr::first(n_links),
      partner_dist_sq = sum((value - centroid)^2),
      .groups = "drop_last"
    ) |>
    dplyr::group_by(site, species) |>
    dplyr::summarise(
      niche_size = weighted.mean(partner_dist_sq, w = n_links), .groups = "drop"
    )

  niche_originality <- niche_position |>
    dplyr::group_by(site, axis) |>
    dplyr::mutate(guild_centroid = mean(centroid)) |>
    dplyr::group_by(site, species) |>
    dplyr::summarise(
      niche_originality = sqrt(sum((centroid - guild_centroid)^2)),
      .groups = "drop"
    )

  dplyr::inner_join(niche_size, niche_originality, by = c("site", "species"))
}
