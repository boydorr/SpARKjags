#' DICtable
#'
#' @export
#'
DICtable <- function(models) {
  data.frame(model = models, stringsAsFactors = F) %>%
    dplyr::group_by(model) %>%
    dplyr::mutate(DIC = DIC(get(model))) %>%
    dplyr::arrange(dplyr::desc(DIC)) %>%
    data.frame() %>%
    flextable::regulartable() %>%
    flextable::autofit()
}
