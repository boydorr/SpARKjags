#' get_vars
#'
#' @export
#'
get_vars <- function(result) {
    tmp <- result %>%
        monitored_variables() %>%
        gsub(",", "", .) %>%
        gsub("[[0-9]*]", "", .)

    if(any(grepl("sd", tmp)))
        tmp[which(grepl("sd", tmp))] <- "sd"

    if(any(grepl("prob.of.", tmp)))
        tmp[which(grepl("sd", tmp))] <- "prob.of"

    tmp %>%
        unique() %>%
        paste(collapse = ")|(") %>%
        paste0("(", ., ")")
}
