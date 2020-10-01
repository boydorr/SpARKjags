#' get_model
#'
#' @param path path
#'
#' @export
#'
get_model <- function(path) {

  assertthat::assert_that(grepl(".rds$", path))

  # Read results
  readRDS(path)
}
