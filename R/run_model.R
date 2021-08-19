#' #' run_model
#' #'
#' #' Run model. If model results already exists, return file location.
#' #'
#' #' @param data a \code{list} containing the data input for the runjags model,
#' #' generated with \code{jags_data()}
#' #' @param location a \code{string} specifying a model location
#' #' @param thin (optional; default = 1) an \code{integer} specifying (from
#' #' \code{runjags::run.jags()}) the thinning interval to be used in JAGS.
#' #' Increasing the thinning interval may reduce autocorrelation, and therefore
#' #' reduce the number of samples required, but will increase the time required to
#' #' run the simulation. Using this option thinning is performed directly in JAGS,
#' #' rather than on an existing MCMC object.
#' #' @param save_to (optional) a \code{string} specifying a custom save location
#' #' for the model output. If save_to is missing, the model output will be
#' #' saved in the same directory as the model script.
#' #'
#' #' @return Returns \code{string} specifying the location of the model output
#' #'
#' run_model <- function(data,
#'                       location,
#'                       thin,
#'                       save_to) {
#'
#'
#' }
