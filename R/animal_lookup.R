#' animal_lookup
#'
#' @export
#'
animal_lookup <- function(data) {
  livestock <- METAdata %>%
    dplyr::filter(Livestock == "yes") %>%
    dplyr::select(ASSOCIATED_SPECIES) %>%
    unique() %>%
    dplyr::filter(!grepl("^Un", ASSOCIATED_SPECIES),
                  !is.na(ASSOCIATED_SPECIES)) %>%
    dplyr::mutate(dataset = "livestock_dat")

  companion <- METAdata %>%
    dplyr::filter(Companion_animal == "yes") %>%
    dplyr::select(ASSOCIATED_SPECIES) %>%
    unique() %>%
    dplyr::mutate(dataset = "companion_dat")

  wild <- METAdata %>%
    dplyr::filter(Wild_animal == "yes") %>%
    dplyr::select(ASSOCIATED_SPECIES) %>%
    unique() %>%
    dplyr::filter(!grepl("^Un", ASSOCIATED_SPECIES)) %>%
    dplyr::mutate(dataset = "wild_dat")

  lookup <- dplyr::bind_rows(livestock, companion, wild) %>%
    dplyr::rename(species = ASSOCIATED_SPECIES) %>%
    dplyr::mutate(group = gsub("_dat", "", dataset))
}