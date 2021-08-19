#' run_SpARKjags_model
#'
#' Run SpARKjags model. If model results already exists, return file location.
#'
#' @param data a \code{list} containing the data input for the runjags model,
#' generated with \code{jags_data()}
#' @param SpARKjags_model a \code{string} specifying a model location within
#' the SpARKjags package. The format should be dir/model (as listed in
#' \code{list_models()}).
#' @param save_to (optional) a \code{string} specifying a custom save location
#' for the model output
#' @param thin (optional; default = 1) an \code{integer} specifying (from
#' \code{runjags::run.jags()}) the thinning interval to be used in JAGS.
#' Increasing the thinning interval may reduce autocorrelation, and therefore
#' reduce the number of samples required, but will increase the time required to
#' run the simulation. Using this option thinning is performed directly in JAGS,
#' rather than on an existing MCMC object.
#'
#' @return Returns code{string} specifying the location of the model output
#' @export
#'
run_SpARKjags_model <- function(data,
                                SpARKjags_model,
                                save_to,
                                thin = 1) {

  if(!grepl(".R$", SpARKjags_model))
    stop("SpARKjags_model must point to an *.R file")

  # Determine where to save results
  filename <- gsub(".r$", "", basename(SpARKjags_model))
  filename <- gsub(".R$", "", filename)
  filename <- paste0(filename, ".rds")
  save_to <- file.path(save_to, filename)

  location <- SpARKjags_model

  # Run model
  # If file already exists, return its path
  if(file.exists(save_to)) {
    message("File already exists")
    return(save_to)

  } else {
    # If directory doesn't exist, create it
    if(!grepl("/", save_to))
      save_to <- paste0("./", save_to)
    directory <- dirname(save_to)
    if(!file.exists(directory))
      dir.create(directory, recursive = TRUE)

    # Run model
    results <- runjags::run.jags(location,
                                 data = data$jags,
                                 n.chains = 2,
                                 silent.jags = T,
                                 thin = thin)

    # If a results summary hasn't been generated, generate it
    if(!results$summary.available) {
      results <- runjags::add.summary(results)
    }

    # Save results
    saveRDS(results, save_to)
    return(save_to)
  }
}
