#' summarise_samples
#'
#' @export
#'
summarise_samples <- function(model, data) {

  df <- import_data(model, data)

  # Summarise number of samples in original dataset
  df %>%
    dplyr::mutate(name = dplyr::case_when(
      name == "Hospital" & clinical == "yes" ~ "Hospital (Clinical)",
      name == "Hospital" & clinical == "no" ~ "Hospital (Carriage)",
      T ~ name)) %>%
    dplyr::group_by(badgroup, name) %>%
    dplyr::summarise(count = n()) %>%
    reshape2::dcast(name ~ badgroup, value.var = "count", fill = 0) %>%
    dplyr::rename(`Bad group` = "1",
                  `Good group` = "0",
                  Category = name) %>%
    flextable::regulartable() %>%
    flextable::autofit()
}