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
    quarto_args = c("--embed-resources"),
    quiet = FALSE
  )
)
