#' get_vars
#'
#' @export
#'
get_vars <- function() {
  res.a %>%
    monitored_variables() %>%
    paste(collapse = ", ") %>%
    paste0("c(", ., ")")
}