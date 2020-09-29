#' get_human
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
get_human <- function(classification,
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


  ward_types <- wards %>%
    dplyr::select(HOSPITAL_WARD, Area) %>%
    dplyr::rename(ward_type = Area) %>%
    dplyr::arrange(ward_type) %>%
    dplyr::mutate(ward_type = dplyr::case_when(
      ward_type == "" ~ paste0("system", dplyr::row_number()),
      T ~ ward_type))

  resistance <- cbind.data.frame(index = c(1, 0),
                                 interpretation = c("R", "S")) %>%
    dplyr::arrange(index)

  hospitals <- cbind.data.frame(index = c(1, 2, 3, 4, -1, -2),
                                hospital = c("San_Matteo", "Montescano",
                                             "Maugeri", "Santa_Margherita",
                                             "Belgioioso", "Pavia_citizen")) %>%
    dplyr::arrange(index)


  # Initialise dataset ------------------------------------------------------


  # Merge METAdata with ward type lookup table
  data <- merge(SpARK::METAdata, ward_types, all.x = T) %>%
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


  # Tidy dataset ------------------------------------------------------


  data %<>%
    dplyr::filter(used_MIC == "yes",
                  Category %in% "human",
                  DATE_OF_BIRTH != "XXXX") %>%
    dplyr::rename(ward = HOSPITAL_WARD,
                  sample_type = SAMPLE_TYPE,
                  hospital = SPECIFIC_GROUP,
                  gender = SEX,
                  clinical = Clinical,
                  ward_type = ward_type,
                  interpretation = Interpretation,
                  antibiotic = Antibiotic_name) %>%
    merge(hospitals, by = "hospital") %>%
    dplyr::mutate(hospital = index,
                  sample_GUID = gsub("_C[1-9]$", "", GUID),
                  age = as.numeric(lubridate::interval(
                    lubridate::ymd(DATE_OF_BIRTH),
                    lubridate::ymd(SAMPLE_DATE)),
                    unit = "years"),
                  sample_month = lubridate::month(
                    lubridate::ymd(SAMPLE_DATE)),
                  sample_season = case_when(
                    sample_month %in% 3:5 ~ "spring",
                    sample_month %in% 6:8 ~ "summer",
                    sample_month %in% 9:11 ~ "autumn",
                    sample_month %in% c(12, 1, 2) ~ "winter")) %>%
    bin_ages(10) %>%
    dplyr::select(GUID, interpretation, bacteria, ST, antibiotic,
                  sample_GUID, sample_type, ward, hospital,
                  gender, clinical, age_group, age_group2,
                  ward_type, sample_month, sample_season, age)

  # Define clinical status
  data %<>% defineClinical()

  # "Remove the sample taken from the Microbiology_and_Virology_Laboratory
  # ward in San Matteo
  delete_sample <- data %>%
    dplyr::filter(hospital != -2, ward_type == "Volunteer") %$%
    GUID %>%
    unique()
  if(length(delete_sample) > 0) data %<>%
    dplyr::filter(GUID != delete_sample)

  # "index", "ward" lookup table
  w <- data %>%
    dplyr::select(hospital, ward) %>%
    unique() %>%
    dplyr::arrange(desc(hospital)) %>%
    dplyr::mutate(index = dplyr::case_when(
      hospital == -1 ~ -1,
      hospital == -2 ~ -2,
      T ~ as.numeric(dplyr::row_number()))) %>%
    dplyr::select(index, ward, -hospital)

  # Merge data with ward lookup table
  data %<>% merge(w, by = "ward") %>%
    dplyr::select(-index)


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
    class_tables <- SpARK::ATBdata %>%
      dplyr::select(Antibiotic_name, Classification) %>%
      unique() %>%
      dplyr::rename(antibiotic = Antibiotic_name)

    # Filter by antibiotic class
    data %>%
      merge(class_tables, all.x = TRUE) %>%
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

    }else if(colnames(data)[x] == "hospital") {
      out <- hospitals

    }else if(colnames(data)[x] == "ward") {
      out <- w

    }else if(colnames(data)[x] == "age") {
      out <- cbind.data.frame(index = data$age, age = data$age)

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


  # Sample_collection centre index
  ind <- lookup_tables$ward_type %>%
    dplyr::filter(ward_type == "Sample_Collection_Center") %$%
    index

  # Hospital samples
  hospital_dat <- dat_numeric %>%
    dplyr::filter(!hospital %in% c(-2, -1),
                  ward_type != ind) %>%
    dplyr::select(GUID, sample_GUID, bacteria, ST, class_interpretation,
                  dplyr::one_of(antibiotics)) %>%
    unique() %>%
    dplyr::arrange(GUID) %>%
    tibble::as_tibble()

  # GP samples (sample collection centre, Belgioioso)
  gp_dat <- dat_numeric %>%
    dplyr::filter(hospital == -1) %>%
    dplyr::select(GUID, sample_GUID, bacteria, ST, class_interpretation,
                  dplyr::one_of(antibiotics)) %>%
    unique() %>%
    dplyr::arrange(GUID) %>%
    tibble::as_tibble()

  # Outpatients (sample collection centre, San Matteo)
  outpatients_dat <- dat_numeric %>%
    dplyr::filter(hospital == 1,
                  ward_type == ind) %>%
    dplyr::select(GUID, sample_GUID, bacteria, ST, class_interpretation,
                  dplyr::one_of(antibiotics)) %>%
    unique() %>%
    dplyr::arrange(GUID) %>%
    tibble::as_tibble()

  # Volunteer samples
  volunteer_dat <- dat_numeric %>%
    dplyr::filter(hospital == -2) %>%
    dplyr::select(GUID, sample_GUID, bacteria, ST, class_interpretation,
                  dplyr::one_of(antibiotics)) %>%
    unique() %>%
    dplyr::arrange(GUID) %>%
    tibble::as_tibble()

  # Ward information
  ward_dat <- dat_numeric %>%
    dplyr::select(ward, hospital, ward_type) %>%
    dplyr::filter(!ward %in% -1:-2) %>%
    unique() %>%
    dplyr::arrange(ward) %>%
    tibble::as_tibble()

  # Sample type information
  sample_type_dat <- dat_numeric %>%
    dplyr::select(sample_type, clinical) %>%
    unique() %>%
    dplyr::arrange(sample_type) %>%
    tibble::as_tibble()

  # General information
  sample_dat <- dat_numeric %>%
    dplyr::select(sample_GUID, sample_month, sample_season, sample_type,
                  ward, age, age_group, age_group2, gender) %>%
    unique() %>%
    dplyr::arrange(sample_GUID) %>%
    tibble::as_tibble()


  # Checks ------------------------------------------------------------------


  assertthat::assert_that(nrow(dat_numeric) == nrow(hospital_dat) +
                            nrow(gp_dat) + nrow(outpatients_dat) +
                            nrow(volunteer_dat))
  assertthat::assert_that(all(sample_dat$sample_GUID == 1:nrow(sample_dat)))
  assertthat::assert_that(all(sample_type_dat$sample_type ==
                                1:nrow(sample_type_dat)))
  assertthat::assert_that(all(ward_dat$ward == 1:nrow(ward_dat)))


  # Output ------------------------------------------------------------------


  data %<>%
    merge(lookup_tables$GUID, all.x = T) %>%
    dplyr::mutate(dataset = dplyr::case_when(
      index %in% hospital_dat$GUID ~ "hosp",
      index %in% gp_dat$GUID ~ "gp",
      index %in% volunteer_dat$GUID ~ "vol",
      index %in% outpatients_dat$GUID ~ "out"))

  df <- list(hospital_dat = hospital_dat,
             gp_dat = gp_dat,
             outpatients_dat = outpatients_dat,
             volunteer_dat = volunteer_dat,
             sample_dat = sample_dat,
             ward_dat = ward_dat,
             sample_type_dat = sample_type_dat,
             lookup_tables = lookup_tables,
             data = data)

  if(classification == "all") {
    df$response <- response
    df$antibiotic_classes <- as.numeric(as.factor(all_classes))
  }
  df

}



#' removeCarbapenem
#'
removeCarbapenem <- function(data) {
  n <- length(unique(data$GUID))
  remove_these <- data %>%
    dplyr::select(GUID, Classification, Interpretation) %>%
    dplyr::filter(Classification == "Carbapenem") %>%
    dplyr::group_by(GUID, Interpretation) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::filter(Interpretation == "R") %$%
    GUID
  data %<>% dplyr::filter(!GUID %in% remove_these)
  assertthat::assert_that(
    length(unique(data$GUID)) == n - length(remove_these))
  data
}


#' onlyCarbapenem
#'
onlyCarbapenem <- function(data) {
  n <- length(unique(data$GUID))
  keep_these <- data %>%
    dplyr::select(GUID, Classification, Interpretation) %>%
    dplyr::filter(Classification == "Carbapenem") %>%
    dplyr::group_by(GUID, Interpretation) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::filter(Interpretation == "R") %$%
    GUID
  data %<>% dplyr::filter(GUID %in% keep_these)
}


#' defineClinical
#'
defineClinical <- function(data) {
  data %>%
    dplyr::mutate(clinical = dplyr::case_when(
      hospital == -2 ~ "no",                              # volunteers
      hospital == 1 &&
        ward_type == "Sample_Collection_Center" ~ "yes",  # outpatients
      hospital == -1 ~ "yes",                             # gp
      T ~ clinical),
      sample_type = dplyr::case_when(
        hospital == -2 & sample_type == "feces" ~ "feces_volunteer",
        T ~ sample_type)) # all volunteers are fecal samples
}








