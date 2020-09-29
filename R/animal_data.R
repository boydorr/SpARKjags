#' animal_data
#'
#' @export
#'
animal_data <- function(data, these_animals) {

  # Extract GUID for each animal species ------------------------------------
  tmp <- data$data.animal
  lookup <- animal_lookup(data)

  if(any(these_animals %in% c("rabbit", "donkey", "equine", "hare")))
    stop("Some species exist in more than one dataset. Please edit animal_data function.")

  animal_res <- lapply(seq_along(these_animals), function(x) {
    if(these_animals[x] %in% lookup$species) {
      # Find dataset
      this_dataset <- lookup %>%
        dplyr::filter(species == these_animals[x]) %$% dataset
    } else if(these_animals[x] %in% lookup$group) {
      this_dataset <- lookup %>%
        dplyr::filter(group == these_animals[x]) %$%
        dataset %>% unique()
    } else stop('animal not found')

    if(these_animals[x] %in% c("livestock", "companion", "wild")) {
      this_group <- tmp[[this_dataset]]
    } else { # Find species in dataset
      species_index <- tmp$lookup_tables$associated_species %>%
        dplyr::filter(associated_species == these_animals[x]) %$% index
      this_group <- tmp[[this_dataset]] %>%
        dplyr::filter(associated_species == species_index)
    }

    this_group %>%
      dplyr::select(GUID) %>%
      dplyr::rename(GUIDindex = GUID) %>%
      dplyr::mutate(name = these_animals[x]) %>%
      merge(tmp$lookup_tables$GUID %>%
              dplyr::rename(GUIDindex = index), all.x = TRUE) %>%
      dplyr::select(GUID, name)
  }) %>%
    do.call(rbind.data.frame, .)

  # Original resistance values
  response <- tmp$response %>%
    as.data.frame() %>%
    tibble::rownames_to_column("index") %>%
    merge(tmp$lookup_tables$GUID, all.x = TRUE) %>%
    select(-index) %>%
    select(GUID, everything())

  # Combine data
  dfGUID <- animal_res %>%
    merge(response, all.x = TRUE)
}