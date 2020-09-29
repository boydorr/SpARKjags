#' filter_pathogen
#'
#' @export
#'
filter_pathogen <- function(data, path, pathogen, kleb) {

  if (path) {
    data %<>%
      dplyr::filter(Phoenix_Organism %in% pathogen,
                    !grepl("NEG", GUID)) %>%
      dplyr::rename(bacteria = Phoenix_Organism) %>%
      left_join(SpARK::KLEBdata %>% dplyr::select(GUID, ST),
                by = "GUID")
  } else {
    data %<>%
      dplyr::filter(!grepl("NEG", GUID)) %>%
      left_join(SpARK::KLEBdata %>% dplyr::select(GUID, species, ST),
                by = "GUID") %>%
      dplyr::filter(species %in% kleb,
                    Clinical != "dontknow") %>%
      dplyr::rename(bacteria = species) %>%
      dplyr::select(-Phoenix_Organism)
  }
  data
}
