#' human_data
#'
#' @export
#'
human_data <- function(data) {
  # Extract data
  tmp <- data$data.human
  p <- cbind.data.frame(GUIDindex = tmp$hospital_dat$GUID, cat = "p",
                        index = seq_along(tmp$hospital_dat$GUID),
                        stringsAsFactors = F)
  gp <- cbind.data.frame(GUIDindex = tmp$gp_dat$GUID, cat = "gp",
                         index = seq_along(tmp$gp_dat$GUID),
                         stringsAsFactors = F)
  v <- cbind.data.frame(GUIDindex = tmp$volunteer_dat$GUID, cat = "v",
                        index = seq_along(tmp$volunteer_dat$GUID),
                        stringsAsFactors = F)
  o <- cbind.data.frame(GUIDindex = tmp$outpatients_dat$GUID, cat = "o",
                        index = seq_along(tmp$outpatients_dat$GUID),
                        stringsAsFactors = F)

  # Lookup table
  tag <- data.frame(cat = c("gp", "p", "o", "v"),
                    name = c("GP", "Hospital", "Outpatients", "Volunteers"),
                    stringsAsFactors = F)

  # Original resistance values
  response <- tmp$response %>%
    as.data.frame() %>%
    tibble::rownames_to_column("index") %>%
    merge(tmp$lookup_tables$GUID, all.x = TRUE) %>%
    select(-index) %>%
    select(GUID, everything())

  # Combine data
  dfGUID <- dplyr::bind_rows(p, gp, v, o) %>%
    merge(tmp$lookup_tables$GUID %>%
            dplyr::rename(GUIDindex = index), all.x = TRUE) %>%
    merge(tag, all.x = TRUE) %>%
    dplyr::select(GUID, name) %>%
    merge(tmp$data %>%
            dplyr::select(GUID, hospital, clinical), all.x = TRUE) %>%
    merge(response, all.x = TRUE)
}