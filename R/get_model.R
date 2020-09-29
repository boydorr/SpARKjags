#' get_model
#'
#' @export
#'
get_model <- function(model_name, directory, knit) {
  if (knit) directory <- paste0("../", directory)
  readRDS(paste0(directory, "/", model_name, ".rds"))
}