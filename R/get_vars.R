#' get_vars
#'
#' @param result result
#'
#' @export
#'
get_vars <- function(result) {
    tmp <- monitored_variables(result)
    tmp <- gsub(",", "", tmp)
    tmp <- gsub("[[0-9]*]", "", tmp)

    if(any(grepl("sd", tmp)))
        tmp[which(grepl("sd", tmp))] <- "sd"

    if(any(grepl("prob.of.", tmp)))
        tmp[which(grepl("prob.of.", tmp))] <- "prob.of"

    tmp <- tmp %>%
        unique() %>%
        paste(collapse = ")|(")
    paste0("(", tmp, ")")
}
