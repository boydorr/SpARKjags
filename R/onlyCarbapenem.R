#' onlyCarbapenem
#'
#' @param data data
#'
onlyCarbapenem <- function(data) {
  n <- length(unique(data$GUID))
  keep_these <- data %>%
    dplyr::select(.data$GUID, .data$Classification, .data$Interpretation) %>%
    dplyr::filter(.data$Classification == "Carbapenem") %>%
    dplyr::group_by(.data$GUID, .data$Interpretation) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::filter(.data$Interpretation == "R")
  keep_these <- keep_these$GUID

  data %>% dplyr::filter(.data$GUID %in% keep_these)
}

