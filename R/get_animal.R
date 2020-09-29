#' get_animal
#'
#' @param classification Antibiotic classification
#' @param pathogen Pathogen name (Phoenix_Organism)
#' @param kleb Pathogen name (Kleborate)
#' @param indeterminate c("R", "S"), default is "I"
#' @param removeCarbapenem default is FALSE; remove Carbapenem resitant samples
#' @param onlyCarbapenem default is FALSE
#' @param removeQuinPen default is TRUE
#'
#' @export
#'
get_animal <- function(classification,
                        pathogen,
                        kleb,
                        indeterminate = "I",
                        removeCarbapenem = F,
                        onlyCarbapenem = F,
                        removeQuinPen = T) {

  if ((missing(pathogen) && missing(kleb)) ||
      ((!missing(pathogen)) && (pathogen == "all"))) {
    pathogen <- c("Klebsiella_pneumoniae",
                  "Klebsiella_aerogenes",
                  "Klebsiella_oxytoca",
                  "Klebsiella_variicola",
                  "Raoultella_planticola",
                  "Raoultella_ornithinolytica",
                  "Klebsiella_spp.",
                  "Raoultella_spp.",
                  "Klebsiella_oxytoca/Raoultella_ornithinolytica",
                  "Enterobacter_aerogenes",
                  "Raoultella_terrigena",
                  "Klebsiella_ozaenae")
    path <- TRUE

  } else if (!missing(kleb) && (kleb == "all"))  {
    kleb <- c("Klebsiella aerogenes",
              "Klebsiella grimontii",
              "Klebsiella huaxiensis",
              "Klebsiella michiganensis",
              "Klebsiella oxytoca",
              "Klebsiella pneumoniae",
              "Klebsiella quasipneumoniae subsp. quasipneumoniae",
              "Klebsiella quasipneumoniae subsp. similipneumoniae",
              "Klebsiella quasivariicola",
              "Klebsiella variicola",
              "Raoultella ornithinolytica",
              "Raoultella planticola",
              "Raoultella terrigena")
    path <- FALSE
  } else
    path <- missing(kleb)

  # Generate Lookup tables --------------------------------------------------


  resistance <- cbind.data.frame(index = c(1, 0),
                                 interpretation = c("R", "S")) %>%
    dplyr::arrange(index)


  # Initialise dataset ------------------------------------------------------


  # Merge METAdata with ward type lookup table
  data <- SpARK::METAdata %>%
    merge(SpARK::ATBdata %>% dplyr::rename(GUID = UNIQUE_SPARK_ID) , all.x = T)

  # Remove samples that are resistant to Carbapenem
  if(removeCarbapenem) data %<>% removeCarbapenem()

  # Remove samples that are not resistant to Carbapenem
  if(onlyCarbapenem) data %<>% onlyCarbapenem()

  # Change interpretation of indeterminate results
  data %<>% dplyr::mutate(Interpretation = dplyr::case_when(
    Interpretation == "I" ~ indeterminate,
    T ~ Interpretation))

  # Filter by pathogen
  data %<>% filter_pathogen(path, pathogen, kleb)

  if(removeQuinPen)
    data %<>% dplyr::filter(Classification != "Quinolone",
                            Classification != "Penicillin")


  # Tidy  columns ------------------------------------------------------


  data %<>%
    dplyr::filter(used_MIC == "yes",
                  Category %in% "animal") %>%
    dplyr::rename(associated_species = ASSOCIATED_SPECIES,
                  associated_group = ASSOCIATED_GROUP,
                  type = TYPE,
                  sample_type = SAMPLE_TYPE,
                  interpretation = Interpretation,
                  antibiotic = Antibiotic_name) %>%
    dplyr::mutate(sample_GUID = gsub("_C[1-9]$", "", GUID),
                  sample_month = lubridate::month(
                    lubridate::ymd(SAMPLE_DATE)),
                  sample_season = case_when(
                    sample_month %in% 3:5 ~ "spring",
                    sample_month %in% 6:8 ~ "summer",
                    sample_month %in% 9:11 ~ "autumn",
                    sample_month %in% c(12, 1, 2) ~ "winter")) %>%
    dplyr::select(GUID, interpretation, bacteria, ST, antibiotic,
                  sample_GUID, sample_type, associated_species,
                  associated_group, type, Livestock, Companion_animal,
                  Wild_animal, sample_month, sample_season)


  # Determine class_interpretation ------------------------------------------
  # (resistance to each class of antibiotics)

  if(classification == "all") {
    # Convert antibiotic interpretation to numeric where R = 1, S = 0, and
    # I or ND = NA
    data %<>%
      dplyr::mutate(
        interpretation = dplyr::case_when(
          interpretation == "R" ~ 1,
          interpretation == "S" ~ 0)) %>%
      # Transform dataset such that there is only one row per GUID
      dplyr::ungroup() %>%
      tidyr::spread(antibiotic, interpretation) %>%
      dplyr::mutate(class_interpretation = NA)
    assertthat::assert_that(length(unique(data$GUID)) == nrow(data))

    class_tables <- SpARK::ATBdata %>%
      dplyr::select(Antibiotic_name, Classification) %>%
      unique() %>%
      dplyr::filter(Antibiotic_name %in% colnames(data))

    all_classes <- unique(class_tables$Classification) %>%
      as.character() %>%
      sort()

    antibiotics <- unique(class_tables$Antibiotic_name) %>%
      as.character() %>%
      sort()

    # "index", "classification" lookup table
    antibiotic_class <- cbind.data.frame(index = seq_along(all_classes),
                                         classification = all_classes,
                                         stringsAsFactors = F) %>%
      dplyr::arrange(index)


    # Generate a matrix of response variables
    # (column for each class_interpretation)
    response <- lapply(seq_along(all_classes), function(x) {
      these_antibiotics <- class_tables %>%
        dplyr::filter(Classification %in% all_classes[x]) %$%
        Antibiotic_name %>%
        as.character()
      tmp <- data %>%
        dplyr::select(GUID, one_of(these_antibiotics)) %>%
        dplyr::mutate(class_interpretation =
                        dplyr::case_when(rowSums(. == 1,
                                                 na.rm = T) > 0 ~ 1,
                                         rowSums(. == 0,
                                                 na.rm = T) > 0 ~ 0)) %>%
        dplyr::select(GUID, class_interpretation)
      colnames(tmp)[2] <- all_classes[x]
      tmp
    }) %>%
      purrr::reduce(full_join, by = "GUID") %>%
      dplyr::arrange(GUID) %>%
      dplyr::select(-GUID) %>%
      as.matrix()


  } else {
    # Filter by antibiotic class
    data %<>%
      dplyr::filter(Classification %in% classification) %>%
      # Determine class_interpretation
      dplyr::mutate(class_interpretation = dplyr::case_when(
        sum(interpretation == "R") > 0 ~ 1,
        sum(interpretation == "S") > 0 ~ 0)) %>%
      # Convert antibiotic interpretation to numeric where R = 1, S = 0, and
      # I or ND = NA
      dplyr::mutate(
        interpretation = dplyr::case_when(
          interpretation == "R" ~ 1,
          interpretation == "S" ~ 0)) %>%
      # Transform dataset such that there is only one row per GUID
      dplyr::ungroup() %>%
      tidyr::spread(antibiotic, interpretation)
    assertthat::assert_that(length(unique(data$GUID)) == nrow(data))

    antibiotics <- SpARK::ATBdata %>%
      dplyr::filter(Classification %in% classification) %$%
      Antibiotic_name %>% unique() %>% as.character()
  }


  # Generate lookup table ---------------------------------------------------


  # Combine all lookup tables
  lookup_tables <- lapply(1:ncol(data), function(x) {

    if(colnames(data)[x] %in% c("class_interpretation", antibiotics)) {
      out <- resistance

    }else {
      out <- data[, x, drop = F] %>%
        dplyr::rename(index = colnames(data)[x])
      if(class(out$index) == "character")
        out %<>% dplyr::mutate(index = as.factor(index))
      if(class(out$index) == "factor")
        out %<>% dplyr::mutate(index = as.numeric(index))
      out <- cbind.data.frame(out, data[, x, drop = F]) %>%
        unique() %>%
        dplyr::arrange(index)
    }
  })
  names(lookup_tables) <- colnames(data)

  if(classification == "all") {
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
    filter(Livestock == "yes") %$% index
  livestock_dat <- dat_numeric %>% # see line 116
    dplyr::filter(Livestock == ind) %>%
    dplyr::select(GUID, sample_GUID, bacteria, ST, class_interpretation,
                  sample_type, associated_species, associated_group,
                  sample_month, sample_season,
                  dplyr::one_of(antibiotics)) %>%
    unique() %>%
    dplyr::arrange(GUID) %>%
    tibble::as_tibble()

  # Companion animal samples
  ind <- lookup_tables$Companion_animal %>%
    filter(Companion_animal == "yes") %$% index
  companion_dat <- dat_numeric %>% # see line 116
    dplyr::filter(Companion_animal == ind) %>%
    dplyr::select(GUID, sample_GUID, bacteria, ST, class_interpretation,
                  sample_type, associated_species, associated_group,
                  sample_month, sample_season,
                  dplyr::one_of(antibiotics)) %>%
    unique() %>%
    dplyr::arrange(GUID) %>%
    tibble::as_tibble()

  # Wild animal samples
  ind <- lookup_tables$Wild_animal %>%
    filter(Wild_animal == "yes") %$% index
  wild_dat <- dat_numeric %>% # see line 116
    dplyr::filter(Wild_animal == ind) %>%
    dplyr::select(GUID, sample_GUID, bacteria, ST, class_interpretation,
                  sample_type, associated_species, associated_group,
                  sample_month, sample_season,
                  dplyr::one_of(antibiotics)) %>%
    unique() %>%
    dplyr::arrange(GUID) %>%
    tibble::as_tibble()


  # Checks ------------------------------------------------------------------


  # assertthat::assert_that(nrow(dat_numeric) == nrow(livestock_dat))


  # Output ------------------------------------------------------------------


  data %<>%
    merge(lookup_tables$GUID, all.x = T) %>%
    dplyr::mutate(dataset = dplyr::case_when(
      index %in% livestock_dat$GUID ~ "livestock",
      index %in% companion_dat$GUID ~ "companion",
      index %in% wild_dat$GUID ~ "wild"))


  df <- list(livestock_dat = livestock_dat,
             companion_dat = companion_dat,
             wild_dat = wild_dat,
             lookup_tables = lookup_tables,
             data = data)

  if(classification == "all") {
    df$response <- response
    df$antibiotic_classes <- as.numeric(as.factor(all_classes))
  }
  df

}
