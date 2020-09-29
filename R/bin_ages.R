#' bin_ages
#'
#' @param x data
#' @param N number of bins
#'
#' @export
#'
bin_ages <- function(x, N) {
  max_age <- ceiling(max(x$age))
  jump <- floor(max_age / (N-1))

  x %>% dplyr::mutate(age_group = ggplot2::cut_number(age, 10),
                      age_group2 = ggplot2::cut_interval(age, 10))
}
