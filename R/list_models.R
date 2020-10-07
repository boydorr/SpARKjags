#' list_models
#'
#' List all models contained within model directories
#'
#' @return Returns a \code{data.frame} listing the model directories in the
#' SpARKjags package as well the models scripts and results contained within
#' them
#'
#' @export
#'
#' @examples
#' list_models()
#'
list_models <- function() {
  package_root <- system.file(package = "SpARKjags")
  model_directories <- dir(package_root, "models", full.names = TRUE)

  output <- lapply(model_directories, function(x) {
    models <- sort(dir(x, ".R"))
    results <- sort(dir(x, ".rds"))
    ind <- gsub(".R", "", models) %in% gsub(".rds", "", results)

    tmp <- data.frame(dir = basename(x),
                      model = models,
                      results = ind) %>%
      dplyr::mutate(results = dplyr::case_when(
        results ~ gsub(".R", ".rds", model)))

  })
  do.call(rbind.data.frame, output)
}
