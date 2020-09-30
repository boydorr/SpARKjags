#' DICtable
#'
#' @param models models
#'
#' @export
#'
DICtable <- function(models) {
  data.frame(model = models, stringsAsFactors = F) %>%
    dplyr::group_by(.data$model) %>%
    dplyr::mutate(DIC = DIC(get(.data$modelmodel))) %>%
    dplyr::arrange(dplyr::desc(.data$modelDIC)) %>%
    data.frame() %>%
    flextable::regulartable() %>%
    flextable::autofit()
}
