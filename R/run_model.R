#' run_model
#'
#' @param data data
#' @param SpARKjags_model model path e.g. goodbad_models/a.R
#' @param thin thin
#' @param user_model user_model
#' @param save_to save_to
#'
#' @export
#'
run_model <- function(data, SpARKjags_model, thin, user_model, save_to) {

  # Check only one model script has been input
  if(!missing(SpARKjags_model) & !missing(user_model))
    stop("Input SpARKjags_model location or user_model location, not both.")

  # Find model script
  if(missing(SpARKjags_model)) {
    location <- user_model
  } else {
    location <- base::system.file(SpARKjags_model,
                                  package = "SpARKjags")
  }

  # If the user doesn't provide a save_to location, save it in the same
  # directory as the model
  if(missing(save_to)) {
    save_to <- gsub(".r$", "", location)
    save_to <- gsub(".R$", "", location)
    save_to <- paste0(save_to, ".rds")
  }
  assertthat::assert_that(grepl(".rds$", save_to))

  # If file already exists, stop
  if(file.exists(save_to))
    stop(paste0("A file already exists at this location: ", save_to))

  # If directory doesn't exist, create it
  if(!grepl("/", save_to))
    save_to <- paste0("./", save_to)
  directory <- dirname(save_to)
  if(!file.exists(directory))
    dir.create(directory, recursive = TRUE)

  # Run model
  if(missing(thin)) {
    results <- runjags::run.jags(location,
                                 data = data$jags,
                                 n.chains = 2,
                                 silent.jags = T)
  } else {
    results <- runjags::run.jags(location,
                                 data = data$jags,
                                 n.chains = 2,
                                 silent.jags = T,
                                 thin = thin)
  }

  # If a results summary hasn't been generated, generate it
  if(!results$summary.available) {
    results <- runjags::add.summary(results)
  }

  # Save results
  saveRDS(results, save_to)
  save_to
}
