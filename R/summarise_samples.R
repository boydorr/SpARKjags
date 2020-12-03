#' summarise_samples
#'
#' Counts the number of samples identified as being in the good group or the
#' bad group, extracted from posterior probabilities
#'
#' @param model a \code{runjags} object containing model results
#' @param data data
#'
#' @export
#'
summarise_samples <- function(model, data) {

  # data.frame listing SpARK samples, and their resistances to each antibiotic
  # class, as well as the posterior probability of being in the bad group
  # (mean.p.bad), which defines the badgroup (1 if mean.p.bad > 0.5)
  df <- import_data(model, data)

  if(!any(colnames(df) %in% "badgroup"))
    stop(paste("This function counts the number of samples identified as being",
               "in the good group or the bad group, extracted from posterior",
               "probabilities."))

  # Summarise number of samples in original dataset
  df %>%
    dplyr::mutate(name = dplyr::case_when(
      name == "Hospital" & clinical == "yes" ~ "Hospital (Clinical)",
      name == "Hospital" & clinical == "no" ~ "Hospital (Carriage)",
      T ~ name)) %>%
    dplyr::group_by(.data$badgroup, .data$name) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    reshape2::dcast(name ~ badgroup, value.var = "count", fill = 0) %>%
    dplyr::rename(`Bad group` = "1",
                  `Good group` = "0",
                  Category = .data$name) %>%
    flextable::regulartable() %>%
    flextable::autofit()
}
