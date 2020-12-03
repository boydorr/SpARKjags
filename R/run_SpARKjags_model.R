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

  # Find model script
  if(grepl("/", SpARKjags_model)) {
    location <- system.file(SpARKjags_model,
                            package = "SpARKjags")
    if(length(location) == 0)
      stop("SpARKjags_model not found")
    if(length(location) > 1)
      stop("More than one model was found")
  }

  # Determine where to save results
  filename <- gsub(".r$", "", basename(location))
  filename <- gsub(".R$", "", filename)
  filename <- paste0(filename, ".rds")
  save_to <- file.path(save_to, filename)

  # Run model
  run_model(data = data,
            location = location,
            thin = thin,
            save_to = save_to)
}
