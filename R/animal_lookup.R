#' animal_lookup
#'
#' @param data data
#'
animal_lookup <- function(data) {
  livestock <- SpARK::METAdata %>%
    dplyr::filter(.data$Livestock == "yes") %>%
    dplyr::select(.data$ASSOCIATED_SPECIES) %>%
    unique() %>%
    dplyr::filter(!grepl("^Un", .data$ASSOCIATED_SPECIES),
                  !is.na(.data$ASSOCIATED_SPECIES)) %>%
    dplyr::mutate(dataset = "livestock_dat")

  companion <- SpARK::METAdata %>%
    dplyr::filter(.data$Companion_animal == "yes") %>%
    dplyr::select(.data$ASSOCIATED_SPECIES) %>%
    unique() %>%
    dplyr::mutate(dataset = "companion_dat")

  wild <- SpARK::METAdata %>%
    dplyr::filter(.data$Wild_animal == "yes") %>%
    dplyr::select(.data$ASSOCIATED_SPECIES) %>%
    unique() %>%
    dplyr::filter(!grepl("^Un", .data$ASSOCIATED_SPECIES)) %>%
    dplyr::mutate(dataset = "wild_dat")

  lookup <- dplyr::bind_rows(livestock, companion, wild) %>%
    dplyr::rename(species = .data$ASSOCIATED_SPECIES) %>%
    dplyr::mutate(group = gsub("_dat", "", .data$dataset))
}
