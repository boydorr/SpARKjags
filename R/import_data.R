#' import_data
#'
#' @param model model
#' @param data data
#'
import_data <- function(model, data) {

  # Animal lookup -----------------------------------------------------------

  # List all species in SpARK dataset
  species <- SpARK::METAdata %>%
    dplyr::select(.data$ASSOCIATED_SPECIES) %>% unique() %>%
    dplyr::filter(!grepl("^Un", .data$ASSOCIATED_SPECIES),
                  !is.na(.data$ASSOCIATED_SPECIES)) %>%
    unlist() %>%
    unname()

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

    # If the posterior probability of a sample being in the bad group is more
    # than 0.5, then assign 1, otherwise assign 0
    output <- all_data %>%
      merge(pp.badgroup, all.y = TRUE) %>%
      dplyr::rename(mean.p.bad = .data$Mean) %>%
      dplyr::mutate(badgroup = dplyr::if_else(.data$mean.p.bad > 0.5, 1, 0))

  } else {
    output <- all_data
  }

  output %>%
    dplyr::mutate(label = dplyr::case_when(
      .data$name == "Hospital" & .data$clinical == "yes" ~ "Hospital\n(Clinical)",
      .data$name == "Hospital" & .data$clinical == "no" ~ "Hospital\n(Carriage)",
      T ~ .data$name))
}


