#' get_human
#'
#' @param classification Antibiotic classification
#' @param pathogen Pathogen name (Kleborate)
#' @param indeterminate c("R", "S"), default is "I"
#' @param removeCarbapenem default is FALSE; remove Carbapenem resitant samples
#' @param onlyCarbapenem default is FALSE
#' @param removeQuinPen default is TRUE
#'
get_human <- function(classification,
                      pathogen,
                      indeterminate = "I",
                      removeCarbapenem = FALSE,
                      onlyCarbapenem = FALSE,
                      removeQuinPen = TRUE) {

  # Generate Lookup tables --------------------------------------------------

  ward_types <- SpARKjags::wards %>%
    dplyr::select(.data$HOSPITAL_WARD, .data$Area) %>%
    dplyr::rename(ward_type = .data$Area) %>%
    dplyr::arrange(.data$ward_type) %>%
    dplyr::mutate(ward_type = dplyr::case_when(
      .data$ward_type == "" ~ paste0("system", dplyr::row_number()),
      T ~ .data$ward_type))

  resistance <- cbind.data.frame(index = c(1, 0),
                                 interpretation = c("R", "S")) %>%
    dplyr::arrange(.data$index)

  hospitals <- cbind.data.frame(index = c(1, 2, 3, 4, -1, -2),
                                hospital = c("San_Matteo", "Montescano",
                                             "Maugeri", "Santa_Margherita",
                                             "Belgioioso", "Pavia_citizen")) %>%
    dplyr::arrange(.data$index)


  # Initialise dataset ------------------------------------------------------

  # Merge METAdata with ward type lookup table
  data <- SpARK::METAdata %>%
    merge(ward_types, all.x = TRUE) %>%
    merge(SpARK::ATBdata %>% dplyr::rename(GUID = .data$UNIQUE_SPARK_ID),
          all.x = TRUE)

  # Remove samples that are resistant to Carbapenem
  if(removeCarbapenem) data <- data %>% removeCarbapenem()

  # Remove samples that are not resistant to Carbapenem
  if(onlyCarbapenem) data <- data %>% onlyCarbapenem()

  # Change interpretation of indeterminate results
  data <- data %>% dplyr::mutate(Interpretation = dplyr::case_when(
    Interpretation == "I" ~ indeterminate,
    TRUE ~ Interpretation))

  # Filter by pathogen
  data <- data %>% filter_pathogen(pathogen)

  if(nrow(data) == 0) stop("No human data remaining")

  # Remove Quinolone and Penicillin
  if(removeQuinPen)
    data <- data %>% dplyr::filter(.data$Classification != "Quinolone",
                                   .data$Classification != "Penicillin")

  if(nrow(data) == 0) stop("No human data remaining")

  # Tidy dataset ------------------------------------------------------


  data <- data %>%
    dplyr::filter(.data$used_MIC == "yes",
                  .data$Category %in% "human",
                  .data$DATE_OF_BIRTH != "XXXX")

  if(nrow(data) == 0) stop("No human data remaining")

  data <- data %>%
    dplyr::rename(ward = .data$HOSPITAL_WARD,
                  sample_type = .data$SAMPLE_TYPE,
                  hospital = .data$SPECIFIC_GROUP,
                  gender = .data$SEX,
                  clinical = .data$Clinical,
                  ward_type = .data$ward_type,
                  interpretation = .data$Interpretation,
                  antibiotic = .data$Antibiotic_name) %>%
    merge(hospitals, by = "hospital") %>%
    dplyr::mutate(hospital = .data$index,
                  sample_GUID = gsub("_C[1-9]$", "", .data$GUID),
                  age = as.numeric(lubridate::interval(
                    lubridate::ymd(.data$DATE_OF_BIRTH),
                    lubridate::ymd(.data$SAMPLE_DATE)),
                    unit = "years"),
                  sample_month = lubridate::month(
                    lubridate::ymd(.data$SAMPLE_DATE)),
                  sample_season = case_when(
                    .data$sample_month %in% 3:5 ~ "spring",
                    .data$sample_month %in% 6:8 ~ "summer",
                    .data$sample_month %in% 9:11 ~ "autumn",
                    .data$sample_month %in% c(12, 1, 2) ~ "winter")) %>%
    bin_ages(10) %>%
    dplyr::select(.data$GUID, .data$interpretation, .data$bacteria,
                  .data$ST, .data$antibiotic, .data$sample_GUID,
                  .data$sample_type, .data$ward, .data$hospital,
                  .data$gender, .data$clinical, .data$age_group,
                  .data$age_group2, .data$ward_type, .data$sample_month,
                  .data$sample_season, .data$age)

  # Define clinical status
  data <- data %>% defineClinical()

  # "Remove the sample taken from the Microbiology_and_Virology_Laboratory
  # ward in San Matteo
  delete_sample <- data %>%
    dplyr::filter(.data$hospital != -2,
                  .data$ward_type == "Volunteer")
  delete_sample <- unique(delete_sample$GUID)

  if(length(delete_sample) > 0) data <- data %>%
    dplyr::filter(.data$GUID != delete_sample)

  # "index", "ward" lookup table
  w <- data %>%
    dplyr::select(.data$hospital, .data$ward) %>%
    unique() %>%
    dplyr::arrange(desc(.data$hospital)) %>%
    dplyr::mutate(index = dplyr::case_when(
      .data$hospital == -1 ~ -1,
      .data$hospital == -2 ~ -2,
      TRUE ~ as.numeric(dplyr::row_number()))) %>%
    dplyr::select(.data$index, .data$ward, -.data$hospital)

  # Merge data with ward lookup table
  data <- data %>%
    merge(w, by = "ward") %>%
    dplyr::select(-.data$index)


  # Determine class_interpretation ------------------------------------------
  # (resistance to each class of antibiotics)

  if(any(classification == "all")) {
    # Convert antibiotic interpretation to numeric where R = 1, S = 0, and
    # I or ND = NA
    data <- data %>%
      dplyr::mutate(
        interpretation = dplyr::case_when(
          .data$interpretation == "R" ~ 1,
          .data$interpretation == "S" ~ 0)) %>%
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
  lookup_tables <- lapply(seq_len(ncol(data)), function(x) {

    if(colnames(data)[x] %in% c("class_interpretation", antibiotics)) {
      out <- resistance

    }else if(colnames(data)[x] == "hospital") {
      out <- hospitals

    }else if(colnames(data)[x] == "ward") {
      out <- w

    }else if(colnames(data)[x] == "age") {
      out <- cbind.data.frame(index = data$age, age = data$age)

    }else {
      out <- data[, x, drop = FALSE] %>%
        dplyr::rename(index = colnames(data)[x])
      if(class(out$index) == "character")
        out <- out %>% dplyr::mutate(index = as.factor(.data$index))
      if(class(out$index) == "factor")
        out <- out %>% dplyr::mutate(index = as.numeric(.data$index))
      out <- cbind.data.frame(out, data[, x, drop = FALSE]) %>%
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

  tmp <- data %>%
    dplyr::select(.data$sample_type, .data$clinical) %>%
    dplyr::filter(.data$sample_type == "rectal_swab") %>%
    unique() %>% nrow()

  if(tmp == 2) {
    data <- data %>%
      dplyr::mutate(sample_type = dplyr::case_when(
        .data$sample_type == "rectal_swab" &
          .data$clinical == "yes" ~ "rectal_swab_clinical",
        .data$sample_type == "rectal_swab" &
          .data$clinical == "no" ~ "rectal_swab_carriage",
        TRUE ~ .data$sample_type)) %>%
      dplyr::mutate(sample_type = dplyr::case_when(
        .data$sample_type == "skin_swab" &
          .data$clinical == "yes" ~ "skin_swab_clinical",
        .data$sample_type == "skin_swab" &
          .data$clinical == "no" ~ "skin_swab_carriage",
        TRUE ~ .data$sample_type))
  }

  dat_numeric <- data %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    dplyr::mutate_if(is.factor, as.numeric)


  # Generate tables for JAGS models -----------------------------------------


  # Sample_collection centre index
  ind <- lookup_tables$ward_type %>%
    dplyr::filter(.data$ward_type == "Sample_Collection_Center")
  ind <- ind$index

  # Hospital samples
  hospital_dat <- dat_numeric %>%
    dplyr::filter(!.data$hospital %in% c(-2, -1),
                  .data$ward_type != ind) %>%
    dplyr::select(.data$GUID, .data$sample_GUID, .data$bacteria, .data$ST,
                  .data$class_interpretation,
                  dplyr::one_of(antibiotics)) %>%
    unique() %>%
    dplyr::arrange(.data$GUID)

  # GP samples (sample collection centre, Belgioioso)
  gp_dat <- dat_numeric %>%
    dplyr::filter(.data$hospital == -1) %>%
    dplyr::select(.data$GUID, .data$sample_GUID, .data$bacteria, .data$ST,
                  .data$class_interpretation,
                  dplyr::one_of(antibiotics)) %>%
    unique() %>%
    dplyr::arrange(.data$GUID)

  # Outpatients (sample collection centre, San Matteo)
  outpatients_dat <- dat_numeric %>%
    dplyr::filter(.data$hospital == 1,
                  .data$ward_type == ind) %>%
    dplyr::select(.data$GUID, .data$sample_GUID, .data$bacteria, .data$ST,
                  .data$class_interpretation,
                  dplyr::one_of(antibiotics)) %>%
    unique() %>%
    dplyr::arrange(.data$GUID)

  # Volunteer samples
  volunteer_dat <- dat_numeric %>%
    dplyr::filter(.data$hospital == -2) %>%
    dplyr::select(.data$GUID, .data$sample_GUID, .data$bacteria, .data$ST,
                  .data$class_interpretation,
                  dplyr::one_of(antibiotics)) %>%
    unique() %>%
    dplyr::arrange(.data$GUID)

  # Ward information
  ward_dat <- dat_numeric %>%
    dplyr::select(.data$ward, .data$hospital, .data$ward_type) %>%
    dplyr::filter(!.data$ward %in% -1:-2) %>%
    unique() %>%
    dplyr::arrange(.data$ward)

  # Sample type information
  sample_type_dat <- dat_numeric %>%
    dplyr::select(.data$sample_type, .data$clinical) %>%
    unique() %>%
    dplyr::arrange(.data$sample_type)

  # General information
  sample_dat <- dat_numeric %>%
    dplyr::select(.data$sample_GUID, .data$sample_month, .data$sample_season,
                  .data$sample_type, .data$ward, .data$age, .data$age_group,
                  .data$age_group2, .data$gender) %>%
    unique() %>%
    dplyr::arrange(.data$sample_GUID)


  # Checks ------------------------------------------------------------------

  assertthat::assert_that(nrow(dat_numeric) == nrow(hospital_dat) +
                            nrow(gp_dat) + nrow(outpatients_dat) +
                            nrow(volunteer_dat))
  assertthat::assert_that(all(sample_dat$sample_GUID ==
                                seq_len(nrow(sample_dat))))
  assertthat::assert_that(all(sample_type_dat$sample_type ==
                                seq_len(nrow(sample_type_dat))))
  assertthat::assert_that(all(ward_dat$ward == seq_len(nrow(ward_dat))))


  # Output ------------------------------------------------------------------

  data <- data %>%
    merge(lookup_tables$GUID, all.x = TRUE) %>%
    dplyr::mutate(dataset = dplyr::case_when(
      .data$index %in% hospital_dat$GUID ~ "hosp",
      .data$index %in% gp_dat$GUID ~ "gp",
      .data$index %in% volunteer_dat$GUID ~ "vol",
      .data$index %in% outpatients_dat$GUID ~ "out"))

  df <- list(hospital_dat = hospital_dat,
             gp_dat = gp_dat,
             outpatients_dat = outpatients_dat,
             volunteer_dat = volunteer_dat,
             sample_dat = sample_dat,
             ward_dat = ward_dat,
             sample_type_dat = sample_type_dat,
             lookup_tables = lookup_tables,
             data = data)

  if(any(classification == "all")) {
    df$response <- response
    df$antibiotic_classes <- as.numeric(as.factor(all_classes))
  }
  df

}
