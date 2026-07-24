library(targets)
library(tarchetypes)
source(here::here("R", "analyse-coux.R"))
source(here::here("R", "utils.R"))
tar_option_set(
  packages = c(
    "ggplot2",
    "here",
    "dplyr",
    "tidyr",
    "rdryad"
  )
)

list(
  tar_quarto(
    manuscript,
    "output/manuscript.qmd",
    quarto_args = c("--embed-resources")
  ),
  tar_quarto(
    index_main,
    "index.qmd",
    quarto_args = c("--embed-resources")
  ),
  tar_quarto(
    index_reports,
    "reports/index.qmd",
    quarto_args = c("--embed-resources")
  ),
  tar_quarto(
    report_coux,
    "reports/explore-coux2016.qmd",
    quarto_args = c("--embed-resources"),
    quiet = TRUE
  ),
  tar_target(data_coux_folder, "data/raw/coux2016/", format = "file"),
  tar_target(data_coux_raw, read_coux_data(data_coux_folder)),
  tar_target(data_coux, clean_coux_data(data_coux_raw)),
  tar_target(web_list, get_interaction_matrix(data_coux$interaction)),
  tar_target(net_metrics, compute_network_metrics(web_list)),
  tar_target(
    plant_trait_ready,
    prepare_plant_traits(data_coux$plant_trait)
  ),
  tar_target(
    pollinator_trait_ready,
    prepare_pollinator_traits(data_coux$pollinator_trait)
  ),
  tar_target(
    plant_trait_pcoa,
    compute_trait_pcoa(
      plant_trait_ready,
      asym.bin = which(
        names(plant_trait_ready) %in% c("spring", "summer", "fall", "winter")
      )
    )
  ),
  tar_target(
    pollinator_trait_pcoa,
    compute_trait_pcoa(
      pollinator_trait_ready,
      asym.bin = which(grepl("^larv_", names(pollinator_trait_ready)))
    )
  ),
  tar_target(
    plant_trait_fit,
    fit_trait_ordination(plant_trait_pcoa, plant_trait_ready)
  ),
  tar_target(
    pollinator_trait_fit,
    fit_trait_ordination(pollinator_trait_pcoa, pollinator_trait_ready)
  ),
  tar_target(
    plant_uniqueness,
    compute_uniqueness(plant_trait_ready, data_coux$plant_abundance)
  ),
  tar_target(
    pollinator_uniqueness,
    compute_uniqueness(pollinator_trait_ready, data_coux$pollinator_abundance)
  ),
  tar_target(
    plant_originality,
    compute_originality(plant_trait_pcoa, data_coux$plant_abundance)
  ),
  tar_target(
    pollinator_originality,
    compute_originality(pollinator_trait_pcoa, data_coux$pollinator_abundance)
  ),
  tar_target(
    pollinator_originality_w,
    compute_originality(
      pollinator_trait_pcoa, data_coux$pollinator_abundance,
      weighted = TRUE
    )
  ),
  tar_target(
    pollinator_interaction_niche,
    compute_interaction_niche(web_list, plant_trait_pcoa, level = "higher")
  ),
  tar_target(
    plant_interaction_niche,
    compute_interaction_niche(web_list, pollinator_trait_pcoa, level = "lower")
  )
)
