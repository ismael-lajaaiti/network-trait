library(targets)
library(tarchetypes)
lapply(list.files("R", full.names = TRUE), source)
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
    quarto_args = c("--embed-resources")
  ),
  tar_target(data_coux_folder, "data/raw/coux2016/", format = "file"),
  tar_target(data_coux, read_coux_data(data_coux_folder)),
  tar_target(web_list, get_interaction_matrix(data_coux$interaction)),
  tar_target(net_metrics, compute_network_metrics(web_list))
)
