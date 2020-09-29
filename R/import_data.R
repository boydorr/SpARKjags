#' import_data
#'
#' @export
#'
import_data <- function(model, data) {

  # Animal lookup -----------------------------------------------------------
  species <- SpARK::METAdata %>%
    dplyr::select(ASSOCIATED_SPECIES) %>% unique() %>%
    dplyr::filter(!grepl("^Un", ASSOCIATED_SPECIES),
                  !is.na(ASSOCIATED_SPECIES)) %>%
    unlist()
  names(species) <- NULL
  # species <- data$lookup$associated_species$associated_species
  group <- c("livestock", "companion", "wild")
  animal_lookup <- c(species, group)


  # Original human data -----------------------------------------------------
  human <- human_data(data)


  # Original animal data ----------------------------------------------------

  find_animals <- sapply(animal_lookup, function(x)
    sapply(model$monitor, function(y) grepl(paste0(x, "$"), y))) %>%
    colSums()

  if(sum(find_animals) > 0) {
    these_animals <- names(find_animals[find_animals > 0])
    animal <- animal_data(data, these_animals)

    all_data <- dplyr::bind_rows(human, animal)
  } else {
    all_data <- human
  }


  # Posterior probability of each sample being in the bad group -------------
  if(any(grepl(".bad.", model$monitor))) {
    pp.badgroup <- badgroup_posterior(model, data)

    assertthat::assert_that(nrow(pp.badgroup) == nrow(all_data))

    output <- all_data %>%
      merge(pp.badgroup, all.y = TRUE) %>%
      dplyr::rename(mean.p.bad = Mean) %>%
      dplyr::mutate(badgroup = dplyr::if_else(mean.p.bad > 0.5, 1, 0))

  } else {
    output <- all_data
  }

  output %>%
    dplyr::mutate(label = dplyr::case_when(
      name == "Hospital" & clinical == "yes" ~ "Hospital\n(Clinical)",
      name == "Hospital" & clinical == "no" ~ "Hospital\n(Carriage)",
      T ~ name))
}


