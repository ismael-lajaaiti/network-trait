#' Download monthly GBIF occurrence counts for pollinators in New Zealand.
#'
#' Not a pipeline target: it queries the GBIF API, so it is run once by hand
#' and its output saved to `data/raw/gbif/`. Each name is matched to the GBIF
#' backbone (accepted name if a synonym), then New Zealand occurrences are
#' counted by month of observation (records without a month are ignored).
#'
#' @param annotation data frame with columns `pollinator` (code) and
#'   `species_name` (binomial), as in `data/annotation/pollinator_annotation.csv`.
#' @param country ISO country code.
#'
#' @return data tibble with one row per pollinator and month: `pollinator`,
#'   `species_name`, `taxon_key`, `match_type`, `month` and `n_records`.
#' @export
fetch_gbif_months <- function(annotation, country = "NZ") {
  api <- "https://api.gbif.org/v1"
  lapply(seq_len(nrow(annotation)), function(k) {
    name <- annotation$species_name[k]
    match <- jsonlite::fromJSON(sprintf(
      "%s/species/match?name=%s&kingdom=Animalia", api, utils::URLencode(name)
    ))
    key <- if (!is.null(match$acceptedUsageKey)) match$acceptedUsageKey else match$usageKey
    occ <- jsonlite::fromJSON(sprintf(
      "%s/occurrence/search?taxonKey=%s&country=%s&facet=month&facetLimit=12&limit=0",
      api, key, country
    ))
    counts <- stats::setNames(rep(0L, 12), 1:12)
    facet <- occ$facets$counts[[1]]
    if (!is.null(facet)) counts[facet$name] <- facet$count
    tibble::tibble(
      pollinator = annotation$pollinator[k], species_name = name,
      taxon_key = key, match_type = match$matchType,
      month = 1:12, n_records = as.integer(counts)
    )
  }) |>
    dplyr::bind_rows()
}

#' Compute phenological breadth from monthly occurrence counts.
#'
#' Breadth is the expected number of distinct months in a random draw of
#' `n` dated records, which corrects for the very unequal number of records
#' between species. Species with fewer than `n` records get `NA`.
#'
#' @param gbif_months data from `fetch_gbif_months()`.
#' @param n number of records drawn.
#'
#' @return data tibble with one row per pollinator: `species`, `n_records`
#'   and `breadth` (months, between 1 and 12).
#' @export
compute_phenological_breadth <- function(gbif_months, n = 10) {
  rarefy <- function(counts) {
    total <- sum(counts)
    if (total < n) return(NA_real_)
    sum(1 - exp(lchoose(total - counts, n) - lchoose(total, n)))
  }
  gbif_months |>
    dplyr::group_by(species = pollinator) |>
    dplyr::summarise(
      breadth = rarefy(n_records),
      n_records = sum(n_records)
    ) |>
    dplyr::relocate(n_records, .after = species)
}
