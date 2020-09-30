#' filter_pathogen
#'
#' @param data data
#' @param pathogen pathogen
filter_pathogen <- function(data, pathogen) {

  data %>%
    dplyr::filter(!grepl("NEG", .data$GUID)) %>%
    left_join(SpARK::KLEBdata %>% dplyr::select(.data$GUID, .data$species,
                                                .data$ST),
              by = "GUID") %>%
    dplyr::filter(.data$species %in% pathogen,
                  .data$Clinical != "dontknow") %>%
    dplyr::rename(bacteria = .data$species) %>%
    dplyr::select(-.data$Phoenix_Organism)
}
