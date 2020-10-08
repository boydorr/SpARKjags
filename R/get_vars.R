#' get_vars
#'
#' @param result result
#' @param ignore vector
#'
get_vars <- function(result, ignore) {
    tmp <- monitored_variables(result)
    tmp <- gsub(",", "", tmp)
    tmp <- gsub("[[0-9]*]", "", tmp)

    if(any(grepl("sd", tmp)))
        tmp[which(grepl("sd", tmp))] <- "sd"

    if(any(grepl("prob.of.", tmp)))
        tmp[which(grepl("prob.of.", tmp))] <- "prob.of"

    if(any(grepl("a.prob", tmp)))
        tmp[which(grepl(".prob", tmp))] <- "a.prob"

    tmp <- tmp %>%
        unique() %>%
        paste(collapse = ")|(")

    paste0("(", tmp, ")")
}
