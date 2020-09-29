#' monitored_variables
#'
#' @export
#'
monitored_variables <- function(model) {
  model$summary[[1]] %>%
    rownames() %>%
    .[-which(. == "deviance")] %>%
    .[-which(grepl("^bad.", .))]
}