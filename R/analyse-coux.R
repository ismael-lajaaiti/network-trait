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
