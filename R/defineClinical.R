#' defineClinical
#'
#' @param data data
#'
defineClinical <- function(data) {
  data %>%
    dplyr::mutate(clinical = dplyr::case_when(
      .data$hospital == -2 ~ "no",                              # volunteers
      .data$hospital == 1 &
        .data$ward_type == "Sample_Collection_Center" ~ "yes",  # outpatients
      .data$hospital == -1 ~ "yes",                             # gp
      T ~ .data$clinical),
      sample_type = dplyr::case_when(
        .data$hospital == -2 & .data$sample_type == "feces" ~ "feces_volunteer",
        T ~ .data$sample_type)) # all volunteers are fecal samples
}
