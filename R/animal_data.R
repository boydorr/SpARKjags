#' animal_data
#'
#' @param data data
#' @param these_animals these_animals
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
        dplyr::filter(.data$species == these_animals[x]) %>% .data$dataset
    } else if(these_animals[x] %in% lookup$group) {
      this_dataset <- lookup %>%
        dplyr::filter(.data$group == these_animals[x]) %>%
        .data$dataset %>% unique()
    } else stop('animal not found')

    if(these_animals[x] %in% c("livestock", "companion", "wild")) {
      this_group <- tmp[[this_dataset]]
    } else { # Find species in dataset
      species_index <- tmp$lookup_tables$associated_species %>%
        dplyr::filter(.data$associated_species == these_animals[x]) %>%
        .data$index
      this_group <- tmp[[this_dataset]] %>%
        dplyr::filter(.data$associated_species == species_index)
    }

    this_group %>%
      dplyr::select(.data$GUID) %>%
      dplyr::rename(GUIDindex = .data$GUID) %>%
      dplyr::mutate(name = these_animals[x]) %>%
      merge(tmp$lookup_tables$GUID %>%
              dplyr::rename(GUIDindex = .data$index), all.x = TRUE) %>%
      dplyr::select(.data$GUID, .data$name)
  })
  animal_res <- do.call(rbind.data.frame, animal_res)

  # Original resistance values
  response <- tmp$response %>%
    as.data.frame()
  response <- cbind.data.frame(response, index = rownames(response)) %>%
    merge(tmp$lookup_tables$GUID, all.x = TRUE) %>%
    select(-.data$index) %>%
    select(.data$GUID, everything())

  # Combine data
  dfGUID <- animal_res %>%
    merge(response, all.x = TRUE)
}
