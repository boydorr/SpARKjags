#' true_resistance
#'
#' @export
#'
true_resistance <- function(model, data) {

  # Posterior probability of each sample being in the bad group
  pp.badgroup <- model$summary[[1]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("index") %>%
    dplyr::filter(grepl("bad", index),
                  !grepl("prob", index)) %>%
    dplyr::mutate(index = gsub("^[a-z]*.", "", index),
                  index = gsub("]", "", index)) %>%
    tidyr::separate(index, c("cat", "index"), sep = "\\[") %>%
    dplyr::mutate(index = as.numeric(index)) %>%
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
    as.data.frame() %>%
    tibble::rownames_to_column("index") %>%
    merge(data$all$lookup_tables$GUID, all.x = TRUE) %>%
    select(-index) %>%
    select(GUID, everything())
  # Original dataset and resistance values
  dfGUID <- rbind.data.frame(p, gp, v, o) %>%
    merge(data$all$lookup_tables$GUID %>%
            rename(GUIDindex = index), all.x = TRUE) %>%
    merge(pp.badgroup, all.y = TRUE) %>%
    dplyr::mutate(cat = dplyr::case_when(cat == "gp" ~ "GP",
                                         cat == "p" ~ "Hospital",
                                         cat == "o" ~ "Outpatients",
                                         cat == "v" ~ "Volunteers")) %>%
    dplyr::select(GUID, cat, Mean) %>%
    dplyr::rename(mean.p.bad = Mean) %>%
    dplyr::mutate(badgroup = dplyr::if_else(mean.p.bad > 0.5, 1, 0)) %>%
    merge(data$all$data, all.x = TRUE) %>%
    dplyr::select(GUID, mean.p.bad, hospital, clinical, cat, badgroup) %>%
    merge(response, all.x = TRUE)

}