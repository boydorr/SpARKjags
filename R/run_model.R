#' run_model
#'
#' Run model. If model results already exists, return file location.
#'
#' @param data a \code{list} output from the jags_data function
#' @param location a \code{string} specifying a model location
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
run_model <- function(data,
                      location,
                      thin,
                      save_to) {

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
