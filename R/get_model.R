#' get_model
#'
#' @param model_name model_name
#' @param directory directory
#'
#' @export
#'
get_model <- function(model_name, directory) {
  # Find results
  location <- system.file(directory, paste0(model_name, ".rds"),
                          package = "SpARKcarbapenem")

  # Read results
  readRDS(location)
}
