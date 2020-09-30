#' get_parameters
#'
#' @param results results
#'
#' @export
#'
get_parameters <- function(results) {
  results %>%
    coda::as.mcmc.list() %>%
    ggmcmc::ggs() %>%
    .data$Parameter %>%
    unique() %>%
    as.character()
}
