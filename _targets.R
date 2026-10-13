library(targets)
library(tarchetypes)
source(here::here("R", "analyse-coux.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "contributions.R"))
source(here::here("R", "sem.R"))
source(here::here("R", "gbif.R"))
source(here::here("R", "simulation.R"))
tar_option_set(
  packages = c(
    "ggplot2",
    "here",
    "dplyr",
    "tidyr",
    "rdryad",
    "INLA",
    "glmmTMB"
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
  tar_quarto(
    report_contributions,
    "reports/species-contributions.qmd",
    quarto_args = c("--embed-resources"),
    quiet = TRUE
  ),
  tar_quarto(
    report_model,
    "reports/trait-matching-model.qmd",
    quarto_args = c("--embed-resources"),
    quiet = TRUE
  ),
  tar_target(data_coux_folder, "data/raw/coux2016/", format = "file"),
  tar_target(data_coux_raw, read_coux_data(data_coux_folder)),
  tar_target(data_coux, clean_coux_data(data_coux_raw)),
  tar_target(
    pollinator_annotation_file,
    "data/annotation/pollinator_annotation.csv",
    format = "file"
  ),
  tar_target(
    pollinator_annotation,
    read_pollinator_annotation(pollinator_annotation_file)
  ),
  tar_target(
    pollinator_abundance,
    remove_abundance_record(
      data_coux$pollinator_abundance,
      site = "blabrolee", species = "l_sordidum"
    )
  ),
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
      ),
      weights = get_trait_weights(
        plant_trait_ready,
        groups = list(
          flowering_season = c("spring", "summer", "fall", "winter")
        )
      )
    )
  ),
  tar_target(
    pollinator_trait_pcoa,
    compute_trait_pcoa(
      pollinator_trait_ready,
      asym.bin = which(grepl("^larv_", names(pollinator_trait_ready))),
      weights = get_trait_weights(
        pollinator_trait_ready,
        groups = list(
          body_size = c("body_length", "body_weigth"),
          food_preference = c("pref_nec", "pref_pol"),
          larval_diet = grep(
            "^larv_", names(pollinator_trait_ready),
            value = TRUE
          )
        )
      )
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
    compute_uniqueness(plant_trait_pcoa, data_coux$plant_abundance)
  ),
  tar_target(
    pollinator_uniqueness,
    compute_uniqueness(
      pollinator_trait_pcoa, pollinator_abundance
    )
  ),
  tar_target(
    plant_originality,
    compute_originality(plant_trait_pcoa, data_coux$plant_abundance)
  ),
  tar_target(
    pollinator_originality,
    compute_originality(pollinator_trait_pcoa, pollinator_abundance)
  ),
  tar_target(
    pollinator_originality_w,
    compute_originality(
      pollinator_trait_pcoa, pollinator_abundance,
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
  ),
  tar_target(
    pollinator_trait_network,
    prepare_trait_network_data(
      pollinator_originality, pollinator_originality_w, pollinator_uniqueness,
      net_metrics, pollinator_abundance, web_list
    )
  ),
  tar_target(
    trait_network_models,
    compare_trait_network_models(pollinator_trait_network)
  ),
  tar_target(
    robustness_contributions,
    compute_robustness_contributions(web_list)
  ),
  tar_target(plant_fd_tree, build_fd_tree(plant_trait_pcoa)),
  tar_target(
    fd_contributions,
    compute_fd_contributions(web_list, plant_fd_tree)
  ),
  tar_target(
    plant_support,
    compute_plant_support(
      web_list, plant_fd_tree, plant_uniqueness, plant_originality
    )
  ),
  tar_target(
    gbif_months_file,
    "data/raw/gbif/pollinator_months_nz.csv",
    format = "file"
  ),
  tar_target(gbif_months, read.csv(gbif_months_file)),
  tar_target(
    pollinator_breadth,
    compute_phenological_breadth(gbif_months, n = 10)
  ),
  tar_target(
    pollinator_social_breadth,
    get_social_breadth(pollinator_trait_ready)
  ),
  tar_target(plant_floral_breadth, get_floral_breadth(plant_trait_ready)),
  tar_target(
    sem_data,
    prepare_sem_data(
      pollinator_trait_network, robustness_contributions, fd_contributions,
      plant_support, pollinator_social_breadth, plant_floral_breadth, web_list
    )
  ),
  tar_target(contribution_sem, fit_contribution_sem(sem_data)),
  tar_target(
    contribution_sem_summary,
    summarise_contribution_sem(contribution_sem, sem_data)
  ),
  tar_target(
    trait_network_models_richness,
    compare_trait_network_models(
      pollinator_trait_network,
      models = list(
        `weighted originality` = c(
          "originality_w_z", "log_abundance_z", "log_n_plants_z"
        ),
        originality = c("originality_z", "log_abundance_z", "log_n_plants_z")
      )
    )
  ),
  tar_target(
    sim_settings,
    expand.grid(
      n_plants = c(5, 10, 15, 20),
      lambda = c(0, 0.25, 0.5, 0.75, 1),
      rho = c(-0.5, 0, 0.5),
      n_rep = 30
    )
  ),
  tar_target(sim_data, simulate_setting(sim_settings), pattern = map(sim_settings)),
  tar_target(sim_sems, fit_simulated_sems(sim_data)),
  tar_target(
    sim_sem_example,
    fit_simulated_sem(
      dplyr::filter(sim_data, lambda == 0, rho == 0, n_plants == 20, rep <= 15)
    )
  ),
  tar_target(sim_position_metrics, compare_position_metrics(sim_data)),
  tar_target(sim_example, {
    set.seed(3)
    simulate_network(12, 8, lambda = 0, rho = 0)
  })
)
