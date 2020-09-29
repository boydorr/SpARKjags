#' get_parameters
#'
#' @export
#'
get_parameters <- function(results) {
  results %>%
    coda::as.mcmc.list() %>%
    ggmcmc::ggs() %$%
    Parameter %>%
    unique() %>%
    as.character()
}