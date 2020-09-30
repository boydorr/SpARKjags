#' monitored_variables
#'
#' @param model model
#'
monitored_variables <- function(model) {
  tmp <- rownames(model$summary[[1]])

  if(any(grepl("deviance", tmp)))
    tmp <- tmp[-which(tmp == "deviance")]

  if(any(grepl("^bad.", tmp)))
    tmp <- tmp[-which(grepl("^bad.", tmp))]

  tmp
}
