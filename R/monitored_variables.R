#' monitored_variables
#'
#' @export
#'
monitored_variables <- function(model) {
  tmp <- rownames(model$summary[[1]])

  if(any(grepl("deviance", tmp)))
    tmp <- tmp[-which(tmp == "deviance")]

  if(any(grepl("^bad.", tmp)))
    tmp <- tmp[-which(grepl("^bad.", .))]

  tmp
}
