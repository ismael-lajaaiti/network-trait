#' Count the extinction cascade triggered by removing a pollinator.
#'
#' Topological coextinction (Memmott et al. 2004): removing a pollinator
#' drives extinct every plant left without a pollinator. In a bipartite
#' network the cascade stops there. A pollinator would go extinct only once
#' all its plants are gone, but a plant goes extinct only once all its
#' pollinators are gone, including that surviving pollinator. So the
#' cascade triggered by pollinator \eqn{i} is exactly the number of plants
#' whose only pollinator is \eqn{i}.
#'
#' @param web interaction matrix, plants as rows and pollinators as
#'   columns, as returned by `get_interaction_matrix()`.
#'
#' @return named integer vector: number of secondary plant extinctions
#'   triggered by removing each pollinator (column) of `web`.
#' @export
compute_cascade <- function(web) {
  links <- web > 0
  sole_partner <- rowSums(links) == 1
  colSums(links[sole_partner, , drop = FALSE])
}

#' Plants alive along a pollinator removal order.
#'
#' Pollinators are removed one at a time in `sequence`; after each removal,
#' plants with no remaining pollinator are extinct.
#'
#' @param web interaction matrix, plants as rows and pollinators as
#'   columns.
#' @param sequence column indices of `web`, in removal order.
#'
#' @return logical matrix, one row per plant and one column per step
#'   (0 to `length(sequence)` pollinators removed): is the plant alive?
#' @export
get_plants_alive <- function(web, sequence) {
  links <- web > 0
  vapply(0:length(sequence), function(k) {
    remaining <- setdiff(seq_len(ncol(links)), sequence[seq_len(k)])
    rowSums(links[, remaining, drop = FALSE]) > 0
  }, logical(nrow(links)))
}

#' Area under a curve running from 0 to 1 pollinators removed.
#'
#' @param y values at each step, 0 to n pollinators removed.
#'
#' @return area computed with the trapezoid rule.
#' @export
compute_curve_area <- function(y) {
  x <- seq(0, 1, length.out = length(y))
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

#' Area under an attack tolerance curve for one pollinator removal order.
#'
#' The curve is the fraction of plants still alive against the fraction of
#' pollinators removed, from (0, 1) to (1, 0).
#'
#' @inheritParams get_plants_alive
#'
#' @return robustness, between 0 and 1.
#' @export
compute_atn_area <- function(web, sequence) {
  compute_curve_area(colMeans(get_plants_alive(web, sequence)))
}

#' Build the plant functional dendrogram.
#'
#' UPGMA clustering of Euclidean distances in the Cailliez-corrected PCoA
#' space, the space in which originality and uniqueness are computed
#' (Petchey & Gaston 2002). One tree for all species, so that branch lengths
#' are comparable across sites.
#'
#' @param pcoa an `ape::pcoa` object, as returned by `compute_trait_pcoa()`.
#'
#' @return an `ape::phylo` object.
#' @export
build_fd_tree <- function(pcoa) {
  ape::as.phylo(stats::hclust(stats::dist(pcoa$vectors.cor), "average"))
}

#' Edge-by-species membership of a dendrogram.
#'
#' @param tree an `ape::phylo` object.
#'
#' @return logical matrix, one row per edge and one column per tip: is the
#'   tip below the edge?
#' @export
get_edge_membership <- function(tree) {
  n_tips <- length(tree$tip.label)
  membership <- matrix(
    FALSE, nrow(tree$edge), n_tips,
    dimnames = list(NULL, tree$tip.label)
  )
  # Edges in postorder: a node's children are filled before the node.
  for (e in ape::postorder(tree)) {
    child <- tree$edge[e, 2]
    if (child <= n_tips) {
      membership[e, child] <- TRUE
    } else {
      below <- tree$edge[, 1] == child
      membership[e, ] <- colSums(membership[below, , drop = FALSE]) > 0
    }
  }
  membership
}

#' Dendrogram-based functional diversity of a set of species.
#'
#' FD (Petchey & Gaston 2002) is the total length of the branches that
#' connect the species: an edge counts if it separates at least one present
#' species from another. A single species has no diversity: FD = 0.
#'
#' @param membership edge-by-species matrix from `get_edge_membership()`,
#'   restricted to the species of interest.
#' @param edge_length branch lengths, in the same order as the rows.
#' @param present logical matrix, one row per species (columns of
#'   `membership`) and one column per community.
#'
#' @return FD of each community (column of `present`).
#' @export
compute_dendrogram_fd <- function(membership, edge_length, present) {
  below <- membership %*% present
  n_present <- matrix(colSums(present), nrow(below), ncol(below), byrow = TRUE)
  colSums(edge_length * (below > 0 & below < n_present))
}

#' Area under the FD curve for one pollinator removal order.
#'
#' The curve is the fraction of the site's plant FD retained against the
#' fraction of pollinators removed, from (0, 1) to (1, 0).
#'
#' @inheritParams get_plants_alive
#' @param tree plant functional dendrogram, from `build_fd_tree()`.
#' @param membership edge-by-species matrix of `tree`, precomputed with
#'   `get_edge_membership()` to save time over many orders.
#'
#' @return area, between 0 and 1.
#' @export
compute_fd_area <- function(web, sequence, tree,
                            membership = get_edge_membership(tree)) {
  membership <- membership[, rownames(web), drop = FALSE]
  fd <- compute_dendrogram_fd(
    membership, tree$edge.length, get_plants_alive(web, sequence)
  )
  compute_curve_area(fd / fd[1])
}

#' Contribution of each pollinator to a removal curve area.
#'
#' \eqn{A^{first}_i} is the area under a curve (`area`) averaged over
#' `nrep` random removal orders with \eqn{i} moved to the front, the same
#' orders for every \eqn{i}. A random order starts with any pollinator with
#' probability \eqn{1/n}, so the expected area \eqn{A} is the mean of the
#' \eqn{A^{first}_i}, and we estimate it so. Then
#' \eqn{\Delta A_i = A - A^{first}_i} sums exactly to zero within a site:
#' it compares \eqn{i} with an average pollinator of its site.
#'
#' @param web_list named list of interaction matrices, plants as rows and
#'   pollinators as columns, as returned by `get_interaction_matrix()`.
#' @param area function of `web` and `sequence` returning the area.
#' @param nrep number of random removal orders per site.
#' @param seed random seed. The same seed gives the same orders whatever
#'   `area`, so contributions to different curves are paired.
#'
#' @return data tibble with one row per pollinator and site: `full`,
#'   `first` and `delta`.
#' @export
compute_removal_contributions <- function(web_list, area, nrep = 1000,
                                          seed = 1) {
  set.seed(seed)
  lapply(web_list, function(web) {
    n <- ncol(web)
    orders <- replicate(nrep, sample(n), simplify = FALSE)
    first <- vapply(seq_len(n), function(i) {
      mean(vapply(
        orders, function(o) area(web, c(i, setdiff(o, i))), numeric(1)
      ))
    }, numeric(1))
    full <- mean(first)
    tibble::tibble(
      species = colnames(web), full = full, first = first,
      delta = full - first
    )
  }) |>
    dplyr::bind_rows(.id = "site")
}

#' Compute each pollinator's contribution to network robustness.
#'
#' \eqn{\Delta R_i = R - R^{first}_i}, see
#' `compute_removal_contributions()`, with \eqn{R} the area under the
#' attack tolerance curve (`compute_atn_area()`). It is half the relative
#' cascade, \eqn{(cascade_i - \bar{c}) / 2P} with \eqn{P} plants, exactly
#' when every plant has a single pollinator, and closely otherwise.
#'
#' @inheritParams compute_removal_contributions
#'
#' @return data tibble with one row per pollinator and site: `cascade`,
#'   `r_full`, `r_first` and `delta_r_first`.
#' @export
compute_robustness_contributions <- function(web_list, nrep = 1000,
                                             seed = 1) {
  cascade <- lapply(web_list, compute_cascade) |>
    lapply(tibble::enframe, name = "species", value = "cascade") |>
    dplyr::bind_rows(.id = "site")
  compute_removal_contributions(web_list, compute_atn_area, nrep, seed) |>
    dplyr::rename(r_full = full, r_first = first, delta_r_first = delta) |>
    dplyr::left_join(cascade, by = c("site", "species")) |>
    dplyr::relocate(cascade, .after = species)
}

#' Compute each pollinator's contribution to plant functional diversity.
#'
#' \eqn{\Delta FD_i = FD - FD^{first}_i}, see
#' `compute_removal_contributions()`, with \eqn{FD} the area under the
#' curve of the fraction of plant FD retained (`compute_fd_area()`).
#'
#' @inheritParams compute_removal_contributions
#' @param tree plant functional dendrogram, from `build_fd_tree()`.
#'
#' @return data tibble with one row per pollinator and site: `fd_full`,
#'   `fd_first` and `delta_fd_first`.
#' @export
compute_fd_contributions <- function(web_list, tree, nrep = 1000, seed = 1) {
  membership <- get_edge_membership(tree)
  area <- function(web, sequence) {
    compute_fd_area(web, sequence, tree, membership)
  }
  compute_removal_contributions(web_list, area, nrep, seed) |>
    dplyr::rename(fd_full = full, fd_first = first, delta_fd_first = delta)
}

#' Describe the plants each pollinator supports.
#'
#' Candidate explanations for \eqn{\Delta FD_i - \Delta R_i}, from the
#' plants a pollinator visits and the plants it alone visits (those lost
#' with it, see `compute_cascade()`).
#'
#' @inheritParams compute_removal_contributions
#' @param tree plant functional dendrogram, from `build_fd_tree()`.
#' @param plant_uniqueness plant uniqueness per site, as returned by
#'   `compute_uniqueness()`.
#'
#' @return data tibble with one row per pollinator and site:
#'   `plant_loss` and `fd_loss`, the share of the site's plants and of its
#'   plant FD lost with the plants \eqn{i} alone visits; `excess_fd_loss`,
#'   their difference; `sole_uniqueness` and `partner_uniqueness`, the mean
#'   uniqueness of the plants \eqn{i} alone visits (`NA` if none) and of
#'   all plants it visits.
#' @export
compute_plant_support <- function(web_list, tree, plant_uniqueness) {
  membership <- get_edge_membership(tree)
  fd <- function(species) {
    present <- matrix(tree$tip.label %in% species, ncol = 1)
    compute_dendrogram_fd(membership, tree$edge.length, present)
  }
  mean_uniqueness <- function(s, species) {
    if (length(species) == 0) return(NA_real_)
    keep <- plant_uniqueness$site == s & plant_uniqueness$species %in% species
    mean(plant_uniqueness$uniqueness[keep])
  }
  lapply(names(web_list), function(s) {
    links <- web_list[[s]] > 0
    plants <- rownames(links)
    lapply(colnames(links), function(sp) {
      sole <- plants[rowSums(links) == 1 & links[, sp]]
      tibble::tibble(
        site = s,
        species = sp,
        plant_loss = length(sole) / length(plants),
        fd_loss = 1 - fd(setdiff(plants, sole)) / fd(plants),
        sole_uniqueness = mean_uniqueness(s, sole),
        partner_uniqueness = mean_uniqueness(s, plants[links[, sp]])
      )
    }) |>
      dplyr::bind_rows()
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(excess_fd_loss = fd_loss - plant_loss, .after = fd_loss)
}

#' Plot the plant dendrogram before and after a pollinator is removed.
#'
#' Branches counted in the site's FD are highlighted; plants lost with the
#' pollinator (those it alone visits) are shown in orange.
#'
#' @param tree plant functional dendrogram, from `build_fd_tree()`.
#' @param web interaction matrix of the site, plants as rows and
#'   pollinators as columns.
#' @param pollinator name of the pollinator removed.
#' @param color color of the branches counted in FD.
#' @param lost_color color of the plants lost.
#'
#' @return nothing; draws two base-graphics panels.
#' @export
plot_fd_removal <- function(tree, web, pollinator, color = "#2a78d6",
                            lost_color = "#eb6834") {
  membership <- get_edge_membership(tree)
  links <- web > 0
  plants <- rownames(links)
  lost <- plants[rowSums(links) == 1 & links[, pollinator]]
  counted <- function(species) {
    below <- drop(membership %*% (tree$tip.label %in% species))
    below > 0 & below < length(species)
  }
  fd_all <- sum(tree$edge.length[counted(plants)])
  draw <- function(species, title) {
    fd <- sum(tree$edge.length[counted(species)])
    tip_color <- ifelse(
      tree$tip.label %in% species, "black",
      ifelse(tree$tip.label %in% lost, lost_color, "grey75")
    )
    ape::plot.phylo(
      tree,
      edge.color = ifelse(counted(species), color, "grey85"),
      edge.width = ifelse(counted(species), 3, 1),
      tip.color = tip_color, cex = 0.8, label.offset = 0.005,
      no.margin = FALSE, cex.main = 0.9,
      main = sprintf(
        "%s\nFD = %.2f (%.0f%% of FD, %.0f%% of plants)", title, fd,
        100 * fd / fd_all, 100 * length(species) / length(plants)
      )
    )
  }
  old_par <- graphics::par(mfrow = c(1, 2), mar = c(1, 1, 4, 1))
  on.exit(graphics::par(old_par))
  draw(plants, "All plants")
  draw(setdiff(plants, lost), paste("Without", gsub("_", " ", pollinator)))
}
