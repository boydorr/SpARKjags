#' true_resistance
#'
#' @param model model
#' @param data data
#'
true_resistance <- function(model, data) {

  # Posterior probability of each sample being in the bad group
  pp.badgroup <- model$summary[[1]] %>%
    as.data.frame()
  pp.badgroup <- cbind.data.frame(pp.badgroup,
                                  index = rownames(pp.badgroup)) %>%
    dplyr::filter(grepl("bad", .data$index),
                  !grepl("prob", .data$index)) %>%
    dplyr::mutate(index = gsub("^[a-z]*.", "", .data$index),
                  index = gsub("]", "", .data$index)) %>%
    tidyr::separate(.data$index, c("cat", "index"), sep = "\\[") %>%
    dplyr::mutate(index = as.numeric(.data$index)) %>%
    data.frame()
  # Original dataset
  p <- cbind.data.frame(GUIDindex = data$all$hospital_dat$GUID, cat = "p",
                        index = seq_along(data$all$hospital_dat$GUID))
  gp <- cbind.data.frame(GUIDindex = data$all$gp_dat$GUID, cat = "gp",
                         index = seq_along(data$all$gp_dat$GUID))
  v <- cbind.data.frame(GUIDindex = data$all$volunteer_dat$GUID, cat = "v",
                        index = seq_along(data$all$volunteer_dat$GUID))
  o <- cbind.data.frame(GUIDindex = data$all$outpatients_dat$GUID, cat = "o",
                        index = seq_along(data$all$outpatients_dat$GUID))
  # Original resistance values
  response <- data$all$response %>%
    as.data.frame()
  response <- cbind.data.frame(response, index = rownames(response)) %>%
    merge(data$all$lookup_tables$GUID, all.x = TRUE) %>%
    select(-.data$index) %>%
    select(.data$GUID, everything())
  # Original dataset and resistance values
  dfGUID <- rbind.data.frame(p, gp, v, o) %>%
    merge(data$all$lookup_tables$GUID %>%
            rename(GUIDindex = .data$index), all.x = TRUE) %>%
    merge(pp.badgroup, all.y = TRUE) %>%
    dplyr::mutate(cat = dplyr::case_when(cat == "gp" ~ "GP",
                                         cat == "p" ~ "Hospital",
                                         cat == "o" ~ "Outpatients",
                                         cat == "v" ~ "Volunteers")) %>%
    dplyr::select(.data$GUID, .data$cat, .data$Mean) %>%
    dplyr::rename(mean.p.bad = .data$Mean) %>%
    dplyr::mutate(badgroup = dplyr::if_else(.data$mean.p.bad > 0.5, 1, 0)) %>%
    merge(data$all$data, all.x = TRUE) %>%
    dplyr::select(.data$GUID, .data$mean.p.bad, .data$hospital, .data$clinical,
                  .data$cat, .data$badgroup) %>%
    merge(response, all.x = TRUE)

}
