#' Read and format ata from Coux et al. 2016.
#'
#' @param dir path where the data is stored.
#'
#' @return list of data frames.
#' @export
read_coux_data <- function(dir) {
  data <- list() # Init.
  data$interaction <- read.csv(here(dir, "interactions.csv")) |>
    select(-X) |>
    rename(site = Site, pollinator = Pol_sp, plant = Plant_sp, links = Links)
  data$plant_trait <- read.csv(here(dir, "plant_traits.csv")) |>
    rename(plant = X)
  data$plant_abundance <- read.csv(here(dir, "plant_abundances_bin.csv")) |>
    rename(site = X)
  data$pollinator_trait <- read.csv(here(dir, "pollinator_traits.csv")) |>
    rename(pollinator = X)
  data$pollinator_abundance <- read.csv(
    here(dir, "pollinator_abundances.csv")
  ) |> rename(site = Site)
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

compute_network_metrics <- function(web_list) {
  lapply(web_list, function(web) {
    metrics <- web |> bipartite::specieslevel()
    rbind(
      metrics[[1]] |> dplyr::mutate(guild = "pollinator"),
      metrics[[2]] |> dplyr::mutate(guild = "plant")
    )
  }) |>
    tibble::enframe() |>
    dplyr::rename(site = name) |>
    tidyr::unnest()
}
