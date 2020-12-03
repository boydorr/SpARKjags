#' DICtable
#'
#' @param models a \code{vector} of \code{strings}, each specifying the name of
#' a \code{runjags} object
#'
#' @return Returns a \code{flextable} containing DIC values for each model
#' output
#' @export
#'
DICtable <- function(models) {
  data.frame(model = models, stringsAsFactors = F) %>%
    dplyr::group_by(.data$model) %>%
    dplyr::mutate(DIC = DIC(get(.data$model))) %>%
    dplyr::arrange(dplyr::desc(.data$DIC)) %>%
    data.frame() %>%
    flextable::regulartable() %>%
    flextable::autofit()
}
