#' delete_results
#'
#' Delete SpARKjags model results from within the SpARKjags package
#'
#' @param SpARKjags_model a \code{string} specifying the model results to be
#' deleted or the directory in which all model results should be deleted
#'
#' @return Returns TRUE for each file deleted
#'
#' @examples
#' \dontrun{
#' # Delete single model results
#' delete_results("individual_models/test.rds")
#'
#' # Delete all model results in a particular directory
#' delete_results("individual_models")
#' }
#'
delete_results <- function(SpARKjags_model) {
  # Delete model results
  if(grepl("/", SpARKjags_model) & grepl(".rds$", SpARKjags_model)) {
    location <- system.file(SpARKjags_model, package = "SpARKjags")
    if(length(location) == 0) stop("File not found")
    file.remove(location)

  } else {
    # Delete directory results
    location <- system.file(SpARKjags_model, package = "SpARKjags")
    if(length(location) == 0) stop("Directory not found")
    files <- dir(location, ".rds", full.names = TRUE)
    file.remove(files)
  }
}
