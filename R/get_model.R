#' get_model
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
