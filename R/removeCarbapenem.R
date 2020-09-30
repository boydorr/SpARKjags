#' removeCarbapenem
#'
#' @param data data
#'
removeCarbapenem <- function(data) {
  n <- length(unique(data$GUID))
  remove_these <- data %>%
    dplyr::select(.data$GUID, .data$Classification, .data$Interpretation) %>%
    dplyr::filter(.data$Classification == "Carbapenem") %>%
    dplyr::group_by(.data$GUID, .data$Interpretation) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::filter(.data$Interpretation == "R")
  remove_these <- remove_these$GUID
  data <- data %>% dplyr::filter(!.data$GUID %in% remove_these)

  assertthat::assert_that(
    length(unique(data$GUID)) == n - length(remove_these))

  data
}
