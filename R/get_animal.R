#' get_animal
#'
#' @param classification Antibiotic classification
#' @param pathogen Pathogen name (Kleborate)
#' @param indeterminate c("R", "S"), default is "I"
#' @param removeCarbapenem default is FALSE; remove Carbapenem resitant samples
#' @param onlyCarbapenem default is FALSE
#' @param removeQuinPen default is TRUE
#'
get_animal <- function(classification,
                       pathogen,
                       indeterminate = "I",
                       removeCarbapenem = FALSE,
                       onlyCarbapenem = FALSE,
                       removeQuinPen = TRUE) {

  # Generate Lookup tables --------------------------------------------------

  resistance <- cbind.data.frame(index = c(1, 0),
                                 interpretation = c("R", "S")) %>%
    dplyr::arrange(.data$index)


  # Initialise dataset ------------------------------------------------------

  # Merge METAdata with ward type lookup table
  data <- SpARK::METAdata %>%
    merge(SpARK::ATBdata %>% dplyr::rename(GUID = .data$UNIQUE_SPARK_ID),
          all.x = TRUE)

  # Remove samples that are resistant to Carbapenem
  if(removeCarbapenem) data <- data %>% removeCarbapenem()

  # Remove samples that are not resistant to Carbapenem
  if(onlyCarbapenem) data <- data %>% onlyCarbapenem()

  # Change interpretation of indeterminate results
  data <- data %>% dplyr::mutate(Interpretation = dplyr::case_when(
    .data$Interpretation == "I" ~ indeterminate,
    TRUE ~ Interpretation))

  # Filter by pathogen
  data <- data %>% filter_pathogen(pathogen)

  if(nrow(data) == 0) stop("No animal data remaining")

  if(removeQuinPen)
    data <- data %>% dplyr::filter(.data$Classification != "Quinolone",
                                   .data$Classification != "Penicillin")

  if(nrow(data) == 0) stop("No animal data remaining")

  # Tidy  columns ------------------------------------------------------

  data <- data %>%
    dplyr::filter(.data$used_MIC == "yes",
                  .data$Category %in% "animal")

  if(nrow(data) == 0) stop("No animal data remaining")

  data <- data %>%
    dplyr::rename(associated_species = .data$ASSOCIATED_SPECIES,
                  associated_group = .data$ASSOCIATED_GROUP,
                  type = .data$TYPE,
                  sample_type = .data$SAMPLE_TYPE,
                  interpretation = .data$Interpretation,
                  antibiotic = .data$Antibiotic_name) %>%
    dplyr::mutate(sample_GUID = gsub("_C[1-9]$", "", .data$GUID),
                  sample_month = lubridate::month(
                    lubridate::ymd(.data$SAMPLE_DATE)),
                  sample_season = case_when(
                    sample_month %in% 3:5 ~ "spring",
                    sample_month %in% 6:8 ~ "summer",
                    sample_month %in% 9:11 ~ "autumn",
                    sample_month %in% c(12, 1, 2) ~ "winter")) %>%
    dplyr::select(.data$GUID, .data$interpretation, .data$bacteria, .data$ST,
                  .data$antibiotic, .data$sample_GUID, .data$sample_type,
                  .data$associated_species, .data$associated_group, .data$type,
                  .data$Livestock, .data$Companion_animal, .data$Wild_animal,
                  .data$sample_month, .data$sample_season)


  # Determine class_interpretation ------------------------------------------
  # (resistance to each class of antibiotics)

  if(any(classification == "all")) {
    # Convert antibiotic interpretation to numeric where R = 1, S = 0, and
    # I or ND = NA
    data <- data %>%
      dplyr::mutate(
        interpretation = dplyr::case_when(
          interpretation == "R" ~ 1,
          interpretation == "S" ~ 0)) %>%
      # Transform dataset such that there is only one row per GUID
      dplyr::ungroup() %>%
      tidyr::spread(.data$antibiotic, .data$interpretation) %>%
      dplyr::mutate(class_interpretation = NA)
    assertthat::assert_that(length(unique(data$GUID)) == nrow(data))

    class_tables <- SpARK::ATBdata %>%
      dplyr::select(.data$Antibiotic_name, .data$Classification) %>%
      unique() %>%
      dplyr::filter(.data$Antibiotic_name %in% colnames(data))

    all_classes <- unique(class_tables$Classification) %>%
      as.character() %>%
      sort()

    antibiotics <- unique(class_tables$Antibiotic_name) %>%
      as.character() %>%
      sort()

    # "index", "classification" lookup table
    antibiotic_class <- cbind.data.frame(index = seq_along(all_classes),
                                         classification = all_classes,
                                         stringsAsFactors = FALSE) %>%
      dplyr::arrange(.data$index)


    # Generate a matrix of response variables
    # (column for each class_interpretation)
    response <- lapply(seq_along(all_classes), function(x) {
      these_antibiotics <- class_tables %>%
        dplyr::filter(.data$Classification %in% all_classes[x])
      these_antibiotics <- these_antibiotics$Antibiotic_name %>%
        as.character()
      tmp <- data %>%
        dplyr::select(.data$GUID, one_of(these_antibiotics)) %>%
        dplyr::mutate(class_interpretation =
                        dplyr::case_when(rowSums(. == 1,
                                                 na.rm = TRUE) > 0 ~ 1,
                                         rowSums(. == 0,
                                                 na.rm = TRUE) > 0 ~ 0)) %>%
        dplyr::select(.data$GUID, .data$class_interpretation)
      colnames(tmp)[2] <- all_classes[x]
      tmp
    }) %>%
      purrr::reduce(full_join, by = "GUID") %>%
      dplyr::arrange(.data$GUID) %>%
      dplyr::select(-.data$GUID) %>%
      as.matrix()


  } else {
    class_tables <- SpARK::ATBdata %>%
      dplyr::select(.data$Antibiotic_name, .data$Classification) %>%
      unique() %>%
      dplyr::rename(antibiotic = .data$Antibiotic_name)

    # Filter by antibiotic class
    data <- data %>%
      merge(class_tables, all.x = TRUE) %>%
      dplyr::filter(.data$Classification %in% classification) %>%
      # Determine class_interpretation
      dplyr::mutate(class_interpretation = dplyr::case_when(
        sum(.data$interpretation == "R") > 0 ~ 1,
        sum(.data$interpretation == "S") > 0 ~ 0)) %>%
      # Convert antibiotic interpretation to numeric where R = 1, S = 0, and
      # I or ND = NA
      dplyr::mutate(
        interpretation = dplyr::case_when(
          .data$interpretation == "R" ~ 1,
          .data$interpretation == "S" ~ 0)) %>%
      # Transform dataset such that there is only one row per GUID
      dplyr::ungroup() %>%
      tidyr::spread(.data$antibiotic, .data$interpretation)

    if(length(classification) == 1)
      assertthat::assert_that(length(unique(data$GUID)) == nrow(data))

    antibiotics <- SpARK::ATBdata %>%
      dplyr::filter(.data$Classification %in% classification)
    antibiotics <- antibiotics$Antibiotic_name %>% unique() %>% as.character()
  }


  # Generate lookup table ---------------------------------------------------


  # Combine all lookup tables
  lookup_tables <- lapply(1:ncol(data), function(x) {

    if(colnames(data)[x] %in% c("class_interpretation", antibiotics)) {
      out <- resistance

    }else {
      out <- data[, x, drop = FALSE] %>%
        dplyr::rename(index = colnames(data)[x])
      if(class(out$index) == "character")
        out <- out %>% dplyr::mutate(index = as.factor(.data$index))
      if(class(out$index) == "factor")
        out <- out %>% dplyr::mutate(index = as.numeric(.data$index))
      out <- cbind.data.frame(out, data[, x, drop = F]) %>%
        unique() %>%
        dplyr::arrange(.data$index)
    }
  })
  names(lookup_tables) <- colnames(data)

  if(any(classification == "all")) {
    lookup_tables <- append(lookup_tables, list(antibiotic_class))
    names(lookup_tables)[length(lookup_tables)] <- "antibiotic_class"
  }


  # Convert factors to numeric ----------------------------------------------

  dat_numeric <- data %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    dplyr::mutate_if(is.factor, as.numeric)


  # Generate tables for JAGS models -----------------------------------------

  # Livestock samples
  ind <- lookup_tables$Livestock %>%
    dplyr::filter(.data$Livestock == "yes")
  ind <- ind$index
  livestock_dat <- dat_numeric %>%
    dplyr::filter(.data$Livestock == ind) %>%
    dplyr::select(.data$GUID, .data$sample_GUID, .data$bacteria, .data$ST,
                  .data$class_interpretation, .data$sample_type,
                  .data$associated_species, .data$associated_group,
                  .data$sample_month, .data$sample_season,
                  dplyr::one_of(antibiotics)) %>%
    unique() %>%
    dplyr::arrange(.data$GUID)

  # Companion animal samples
  ind <- lookup_tables$Companion_animal %>%
    dplyr::filter(.data$Companion_animal == "yes")
  ind <- ind$index
  companion_dat <- dat_numeric %>% # see line 116
    dplyr::filter(.data$Companion_animal == ind) %>%
    dplyr::select(.data$GUID, .data$sample_GUID, .data$bacteria, .data$ST,
                  .data$class_interpretation, .data$sample_type,
                  .data$associated_species, .data$associated_group,
                  .data$sample_month, .data$sample_season,
                  dplyr::one_of(antibiotics)) %>%
    unique() %>%
    dplyr::arrange(.data$GUID)

  # Wild animal samples
  ind <- lookup_tables$Wild_animal %>%
    dplyr::filter(.data$Wild_animal == "yes")
  ind <- ind$index
  wild_dat <- dat_numeric %>% # see line 116
    dplyr::filter(.data$Wild_animal == ind) %>%
    dplyr::select(.data$GUID, .data$sample_GUID, .data$bacteria, .data$ST,
                  .data$class_interpretation, .data$sample_type,
                  .data$associated_species, .data$associated_group,
                  .data$sample_month, .data$sample_season,
                  dplyr::one_of(antibiotics)) %>%
    unique() %>%
    dplyr::arrange(.data$GUID)

  # Checks ------------------------------------------------------------------

  # assertthat::assert_that(nrow(dat_numeric) == nrow(livestock_dat))

  # Output ------------------------------------------------------------------

  data <- data %>%
    merge(lookup_tables$GUID, all.x = TRUE) %>%
    dplyr::mutate(dataset = dplyr::case_when(
      .data$index %in% livestock_dat$GUID ~ "livestock",
      .data$index %in% companion_dat$GUID ~ "companion",
      .data$index %in% wild_dat$GUID ~ "wild"))


  df <- list(livestock_dat = livestock_dat,
             companion_dat = companion_dat,
             wild_dat = wild_dat,
             lookup_tables = lookup_tables,
             data = data)

  if(any(classification == "all")) {
    df$response <- response
    df$antibiotic_classes <- as.numeric(as.factor(all_classes))
  }
  df

}
