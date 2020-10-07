#' run_custom_model
#'
#' Run custom model. If model results already exists, return file location.
#'
#' @param data a \code{list} output from the jags_data function
#' @param custom_model a \code{string} specifying a model location
#' @param thin (optional; default = 1) an \code{integer} specifying (from
#' \code{runjags::run.jags()}) the thinning interval to be used in JAGS.
#' Increasing the thinning interval may reduce autocorrelation, and therefore
#' reduce the number of samples required, but will increase the time required to
#' run the simulation. Using this option thinning is performed directly in JAGS,
#' rather than on an existing MCMC object.
#' @param save_to (optional) a \code{string} specifying a custom save location
#' for the model output. If save_to is missing, the model output will be
#' saved in the same directory as the model script.
#'
#' @return Returns code{string} specifying the location of the model output
#'
#' @export
#'
#' @examples
#' \dontrun{
#' run_model(data, "individual_models/a.R")
#' }
#'
run_custom_model <- function(data,
                             custom_model,
                             thin = 1,
                             save_to) {

  if(!grepl(".R$", custom_model))
    stop("custom_model must point to an *.R file")

  # Find model script
  location <- custom_model
  if(!file.exists(location))
    stop("Model not found")

  # If the user doesn't provide a save_to location, save it in the same
  # directory as the model
  if(missing(save_to)) {
    save_to <- gsub(".r$", "", location)
    save_to <- gsub(".R$", "", location)
    save_to <- paste0(save_to, ".rds")
  }
  assertthat::assert_that(grepl(".rds$", save_to))

  # Run model
  run_model(data = data,
            location = location,
            thin = thin,
            save_to = save_to)
}
