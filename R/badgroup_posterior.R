#' badgroup_posterior
#'
#' @export
#'
badgroup_posterior <- function(model, data) {

  # Posterior probability of being in the bad group
  pp.badgroup <- model$summary[[1]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("index") %>%
    dplyr::filter(grepl("bad", index),
                  !grepl("prob", index)) %>%
    dplyr::mutate(index = gsub("^[a-z]*.", "", index),
                  index = gsub("]", "", index)) %>%
    tidyr::separate(index, c("cat", "row_number"), sep = "\\[") %>%
    dplyr::mutate(row_number = as.numeric(row_number)) %>%
    data.frame()


  # Get human GUIDs ---------------------------------------------------------
  GUID.lookup <- data$data.human$lookup_tables$GUID %>%
    dplyr::rename(GUIDindex = index)

  hospital <- data$data.human$hospital_dat %>%
    dplyr::select(GUID) %>%
    dplyr::mutate(cat = "p") %>%
    dplyr::rename(GUIDindex = GUID) %>%
    merge(GUID.lookup, all.x = T) %>%
    tibble::rowid_to_column("row_number")

  gp <- data$data.human$gp_dat %>%
    dplyr::select(GUID) %>%
    dplyr::mutate(cat = "gp") %>%
    dplyr::rename(GUIDindex = GUID) %>%
    merge(GUID.lookup, all.x = T) %>%
    tibble::rowid_to_column("row_number")

  out <- data$data.human$outpatients_dat %>%
    dplyr::select(GUID) %>%
    dplyr::mutate(cat = "o") %>%
    dplyr::rename(GUIDindex = GUID) %>%
    merge(GUID.lookup, all.x = T )%>%
    tibble::rowid_to_column("row_number")

  vol <- data$data.human$volunteer_dat %>%
    dplyr::select(GUID) %>%
    dplyr::mutate(cat = "v") %>%
    dplyr::rename(GUIDindex = GUID) %>%
    merge(GUID.lookup, all.x = T) %>%
    tibble::rowid_to_column("row_number")

  humanGUID <- dplyr::bind_rows(hospital, gp, out, vol) %>%
    dplyr::select(-GUIDindex)

  pp.badhuman <- pp.badgroup %>%
    dplyr::filter(cat %in% c("p", "gp", "v", "o")) %>%
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
          dplyr::filter(species == these_animals[x]) %$% dataset
      } else if(these_animals[x] %in% lookup$group) {
        this_dataset <- lookup %>%
          dplyr::filter(group == these_animals[x]) %$%
          dataset %>% unique()
      } else stop('animal not found')

      this_group <- tmp[[this_dataset]] %>%
        dplyr::rename(GUIDindex = GUID) %>%
        merge(tmp$lookup_tables$GUID %>%
                dplyr::rename(GUIDindex = index), all.x = TRUE)

      # Filter species index
      if(!these_animals[x] %in% c("livestock", "companion", "wild")) {
        species_index <- tmp$lookup_tables$associated_species %>%
          dplyr::filter(associated_species == these_animals[x]) %$% index
        # Output
        this_group %<>%
          dplyr::filter(associated_species %in% species_index) %>%
          dplyr::mutate(cat = these_animals[x])
        # Check
        true_sp <- METAdata %>%
          dplyr::filter(ASSOCIATED_SPECIES == these_animals[x]) %$%
          GUID %>% .[grepl("^SP", .)]
        assertthat::assert_that(all(this_group$GUID %in% true_sp))
      }

      animalGUID <- this_group  %>%
        tibble::rowid_to_column("row_number") %>%
        dplyr::select(row_number, GUID)

      pp.badgroup %>%
        dplyr::filter(cat %in% these_animals[x]) %>%
        merge(animalGUID, all.x = TRUE)
    }) %>%
      do.call(rbind.data.frame, .)

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
    dplyr::mutate(name = dplyr::case_when(is.na(name) ~ cat,
                                          T ~ name)) %>%
    dplyr::select(-cat, -SD, -Naive.SE, -Time.series.SE, -row_number)
}





