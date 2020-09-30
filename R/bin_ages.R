#' bin_ages
#'
#' @param x data
#' @param N number of bins
#'
bin_ages <- function(x, N) {
  max_age <- ceiling(max(x$age))
  jump <- floor(max_age / (N-1))

  x %>% dplyr::mutate(age_group = ggplot2::cut_number(.data$age, 10),
                      age_group2 = ggplot2::cut_interval(.data$age, 10))
}
