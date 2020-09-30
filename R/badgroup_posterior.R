#' badgroup_posterior
#'
#' @param model model
#' @param data data
#'
badgroup_posterior <- function(model, data) {

  # Posterior probability of being in the bad group
  pp.badgroup <- model$summary[[1]] %>%
    as.data.frame()

  pp.badgroup <- cbind.data.frame(pp.badgroup, index = rownames(pp.badgroup)) %>%
    dplyr::filter(grepl("bad", .data$index),
                  !grepl("prob", .data$index)) %>%
    dplyr::mutate(index = gsub("^[a-z]*.", "", .data$index),
                  index = gsub("]", "", .data$index)) %>%
    tidyr::separate(.data$index, c("cat", "row_number"), sep = "\\[") %>%
    dplyr::mutate(row_number = as.numeric(.data$row_number)) %>%
    data.frame()


  # Get human GUIDs ---------------------------------------------------------
  GUID.lookup <- data$data.human$lookup_tables$GUID %>%
    dplyr::rename(GUIDindex = .data$index)

  hospital <- data$data.human$hospital_dat %>%
    dplyr::select(.data$GUID) %>%
    dplyr::mutate(cat = "p") %>%
    dplyr::rename(GUIDindex = .data$GUID) %>%
    merge(GUID.lookup, all.x = T)
  hospital <- cbind.data.frame(hospital, row_number = seq_len(nrow(hospital)))

  gp <- data$data.human$gp_dat %>%
    dplyr::select(.data$GUID) %>%
    dplyr::mutate(cat = "gp") %>%
    dplyr::rename(GUIDindex = .data$GUID) %>%
    merge(GUID.lookup, all.x = T)
  gp <- cbind.data.frame(gp, row_number = seq_len(nrow(gp)))

  out <- data$data.human$outpatients_dat %>%
    dplyr::select(.data$GUID) %>%
    dplyr::mutate(cat = "o") %>%
    dplyr::rename(GUIDindex = .data$GUID) %>%
    merge(GUID.lookup, all.x = T )
  out <- cbind.data.frame(out, row_number = seq_len(nrow(out)))

  vol <- data$data.human$volunteer_dat %>%
    dplyr::select(.data$GUID) %>%
    dplyr::mutate(cat = "v") %>%
    dplyr::rename(GUIDindex = .data$GUID) %>%
    merge(GUID.lookup, all.x = T)
  vol <- cbind.data.frame(vol, row_number = seq_len(nrow(vol)))

  humanGUID <- dplyr::bind_rows(hospital, gp, out, vol) %>%
    dplyr::select(-.data$GUIDindex)

  pp.badhuman <- pp.badgroup %>%
    dplyr::filter(.data$cat %in% c("p", "gp", "v", "o")) %>%
    merge(humanGUID, all.x = T)


  # Get animal GUIDs --------------------------------------------------------
  all.pp <- unique(pp.badgroup$cat)
  these_animals <- all.pp[!all.pp %in% c("p", "gp", "v", "o")]

  # Extract data
  tmp <- data$data.animal
  lookup <- animal_lookup(data)

  if(length(these_animals) > 0) {

    pp.badanimal <- lapply(seq_along(these_animals), function(x) {
      if(these_animals[x] %in% lookup$species) {
        this_dataset <- lookup %>%
          dplyr::filter(.data$species == these_animals[x])
        this_dataset <- this_dataset$dataset
      } else if(these_animals[x] %in% lookup$group) {
        this_dataset <- lookup %>%
          dplyr::filter(.data$group == these_animals[x])
        this_dataset <- this_dataset$dataset %>% unique()
      } else stop('animal not found')

      this_group <- tmp[[this_dataset]] %>%
        dplyr::rename(GUIDindex = .data$GUID) %>%
        merge(tmp$lookup_tables$GUID %>%
                dplyr::rename(GUIDindex = .data$index), all.x = TRUE)

      # Filter species index
      if(!these_animals[x] %in% c("livestock", "companion", "wild")) {
        species_index <- tmp$lookup_tables$associated_species %>%
          dplyr::filter(.data$associated_species == these_animals[x])
        species_index <- species_index$index
        # Output
        this_group <- this_group %>%
          dplyr::filter(.data$associated_species %in% species_index) %>%
          dplyr::mutate(cat = these_animals[x])
        # Check
        true_sp <- SpARK::METAdata %>%
          dplyr::filter(.data$ASSOCIATED_SPECIES == these_animals[x])
        true_sp <- true_sp$GUID[grepl("^SP", true_sp$GUID)]
        assertthat::assert_that(all(this_group$GUID %in% true_sp))
      }

      animalGUID <- this_group

      animalGUID <- cbind.data.frame(animalGUID,
                                     row_number = seq_len(nrow(animalGUID))) %>%
        dplyr::select(.data$row_number, .data$GUID)

      pp.badgroup %>%
        dplyr::filter(.data$cat %in% these_animals[x]) %>%
        merge(animalGUID, all.x = TRUE)
    })
    pp.badanimal <- do.call(rbind.data.frame, pp.badanimal)

    output <- dplyr::bind_rows(pp.badhuman, pp.badanimal)

  } else {
    output <- pp.badhuman
  }

  # Lookup table
  tag <- data.frame(cat = c("gp", "p", "o", "v"),
                    name = c("GP", "Hospital", "Outpatients", "Volunteers"),
                    stringsAsFactors = F)

  output %>%
    merge(tag, all.x = T) %>%
    dplyr::mutate(name = dplyr::case_when(is.na(.data$name) ~ .data$cat,
                                          T ~ .data$name)) %>%
    dplyr::select(-.data$cat, -.data$SD, -.data$Naive.SE, -.data$Time.series.SE,
                  -.data$row_number)
}
