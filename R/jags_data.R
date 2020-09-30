#' jags_data
#'
#' Function to generate data input for jags models
#'
#' @param classification a \code{string} specifying the antibiotic
#' classification, options are:
#' \itemize{
#' \item{"all"} {}
#' }
#' or any combination of the following, in a \code{vector}:
#' \itemize{
#' \item{"Aminoglycoside"} {}
#' \item{"Penicillin Combination"} {}
#' \item{"Penicillin"} {}
#' \item{"Monobactam"} {}
#' \item{"Cephalosporin"} {}
#' \item{"Fluoroquinolone"} {}
#' \item{"Colistin"} {}
#' \item{"Carbapenem"} {}
#' \item{"Fosfomycin"} {}
#' \item{"Penicillin (Penams)"} {}
#' \item{"Quinolone"} {}
#' \item{"Nitrofurantoin"} {}
#' \item{"Tetracycline"} {}
#' \item{"Trimethoprim"} {}
#' \item{"Trimethoprim/Sulfamethoxazole"} {}
#' }
#' @param categories a \code{string} specifying data category, options are:
#' \itemize{
#' \item{"human"}
#' \item{c("human", "animal")}
#' }
#' @param pathogen a \code{string} specifying pathogen name (taken from
#' Kleborate data), options are:
#' \itemize{
#' \item{"all"} {}
#' }
#' or any combination of the following, in a \code{vector}:
#' \itemize{
#' \item{"Klebsiella pneumoniae} {}
#' \item{"Klebsiella quasipneumoniae subsp. quasipneumoniae"} {}
#' \item{"Klebsiella michiganensis"} {}
#' \item{"Klebsiella (Raoultella) planticola"} {}
#' \item{"Klebsiella variicola subsp. variicola"} {}
#' \item{"Klebsiella quasipneumoniae subsp. similipneumoniae"} {}
#' \item{"Klebsiella grimontii"} {}
#' \item{"Klebsiella oxytoca"} {}
#' \item{"Klebsiella (Raoultella) ornithinolytica"} {}
#' \item{"Klebsiella aerogenes"} {}
#' \item{"Klebsiella pasteurii"} {}
#' \item{"Klebsiella spallanzanii"} {}
#' \item{"Klebsiella quasivariicola"} {}
#' \item{"Klebsiella (Raoultella) terrigena"} {}
#' \item{"Klebsiella huaxiensis"} {}
#' \item{"unknown"} - {Removed from possible options}
#' \item{"Raoultella quasiterrigena"} {}
#' }
#' @param ... Additional parameters may be included:
#' \itemize{
#'  \item{indeterminate}: {a \code{string} specifying how indeterminate samples
#'  should be interpreted; options are "R", "S", or "I"; default is "I"}
#'  \item{removeCarbapenem}: {a boolean flag to indicate removal of
#'  Carbapenem resitant samples (\code{TRUE}); default is \code{FALSE}}
#'  \item{onlyCarbapenem}: {a boolean flag to indicate removal of everything
#'  but Carbapenem resitant samples (\code{TRUE}); default is \code{FALSE}}
#'  \item{removeQuinPen}: {a boolean flag to indicate removal of Quinolone and
#'  Penicillin samples (\code{TRUE}); default is \code{TRUE}}
#' }
#'
#' @export
#'
jags_data <- function(classification,
                      categories,
                      pathogen,
                      ...) {

  if(isTRUE(removeCarbapenem) & isTRUE(onlyCarbapenem))
    stop("removeCarbapenem and onlyCarbapenem can't both be true")

  # Ensure classification is valid
  if(any(classification == "all")) {
    assertthat::assert_that(length(classification) == 1)

  } else {
    tmp <- SpARK::ATBdata %>%
      dplyr::select(.data$Classification) %>%
      unlist() %>%
      unique()
    assertthat::assert_that(all(classification %in% tmp))
  }

  # If pathogen == "all" define them, otherwise ensure pathogen is valid
  if(any(pathogen == "all")) {
    assertthat::assert_that(length(pathogen) == 1)

    tmp <- SpARK::KLEBdata %>%
      dplyr::select(.data$species) %>%
      unlist() %>%
      unique()
    pathogen <- tmp[-which(tmp == "unknown")]

  } else {

    tmp <- SpARK::KLEBdata %>%
      dplyr::select(.data$species) %>%
      unlist() %>%
      unique()
    assertthat::assert_that(all(pathogen %in% tmp))
  }

  # human -------------------------------------------------------------------

  if(any(categories %in% c("human", "animal"))) {
    data.human <- get_human(classification, pathogen, ...)

    lookup.human <- data.human$lookup_tables

    data.jags <- list(h_resist = data.human$hospital_dat$class_interpretation,
                      ward = data.human$sample_dat$ward,
                      ward_type = data.human$ward_dat$ward_type,
                      hospital = data.human$ward_dat$hospital,
                      h_GUID = data.human$hospital_dat$GUID,
                      h_sample_GUID = data.human$hospital_dat$sample_GUID,
                      h_bacteria = data.human$hospital_dat$bacteria,
                      gp_resist = data.human$gp_dat$class_interpretation,
                      gp_GUID = data.human$gp_dat$GUID,
                      gp_sample_GUID = data.human$gp_dat$sample_GUID,
                      gp_bacteria = data.human$gp_dat$bacteria,
                      v_resist = data.human$volunteer_dat$class_interpretation,
                      v_GUID = data.human$volunteer_dat$GUID,
                      v_sample_GUID = data.human$volunteer_dat$sample_GUID,
                      v_bacteria = data.human$volunteer_dat$bacteria,
                      o_resist = data.human$outpatients_dat$class_interpretation,
                      o_GUID = data.human$outpatients_dat$GUID,
                      o_sample_GUID = data.human$outpatients_dat$sample_GUID,
                      o_bacteria = data.human$outpatients_dat$bacteria,
                      clinical = data.human$sample_type_dat$clinical,
                      sample_type = data.human$sample_dat$sample_type,
                      N_patients = nrow(data.human$hospital_dat),
                      N_gp = nrow(data.human$gp_dat),
                      N_volunteers = nrow(data.human$volunteer_dat),
                      N_outpatients = nrow(data.human$outpatients_dat),
                      N_sample = nrow(data.human$sample_dat),
                      N_ward = nrow(data.human$ward_dat),
                      gender = data.human$sample_dat$gender,
                      age = data.human$sample_dat$age,
                      age_group = data.human$sample_dat$age_group,
                      age_group2 = data.human$sample_dat$age_group2,
                      N_age_group = length(unique(data.human$sample_dat$age_group)),
                      N_sample_type = nrow(data.human$sample_type_dat),
                      sample_month = data.human$sample_dat$sample_month,
                      sample_season = data.human$sample_dat$sample_season,
                      N_sample_month = length(unique(
                        data.human$sample_dat$sample_month)),
                      N_sample_season = length(unique(
                        data.human$sample_dat$sample_season))
    )

    l_clin <- data.human$lookup_tables$clinical
    data.jags$ncarr       <- l_clin$index[l_clin$clinical=="no"]
    data.jags$nclin       <- l_clin$index[l_clin$clinical=="yes"]
    data.jags$v_clinical  <- data.jags$ncarr # not
    data.jags$gp_clinical <- data.jags$nclin # yes
    data.jags$o_clinical  <- data.jags$nclin # yes

    l_gen <- data.human$lookup_tables$gender
    data.jags$female      <- l_gen$index[l_gen$gender=="F"]
    data.jags$male        <- l_gen$index[l_gen$gender=="M"]
    data.jags$genders     <- c(data.jags$female, data.jags$male)

    data.jags$N_hosp      <- length(unique(data.human$ward_dat$hospital[
      data.human$sample_dat$ward[data.human$hospital_dat$sample_GUID]]))

    data.jags$hosp_wards <- sort(unique(data.jags$ward[data.jags$h_sample_GUID]))
    data.jags$hosp_wardtypes <- sort(unique(data.jags$ward_type[
      data.jags$ward[data.jags$h_sample_GUID]]))
    data.jags$agegroups      <- sort(unique(data.jags$age_group2))
    data.jags$bact_species   <- sort(unique(c(data.jags$h_bacteria,
                                              data.jags$gp_bacteria,
                                              data.jags$v_bacteria,
                                              data.jags$o_bacteria)))
    data.jags$sampletypes    <- sort(unique(data.jags$sample_type))

    if(any(classification == "all")) {
      data.jags$response           <- data.human$response
      data.jags$antibiotic_classes <- length(data.human$antibiotic_classes)
      data.jags$carbapenem         <- data.jags$response[, "Carbapenem"]
      data.jags$carb               <- sort(unique(data.jags$carbapenem))
      data.jags$carb_sus           <- which(data.jags$carb == 0)
      data.jags$carb_res           <- which(data.jags$carb == 1)
      N_classes                    <- length(data.jags$response)
    }

  }


  # animal ------------------------------------------------------------------

  if(any(categories == "animal")) {
    data.animal <- get_animal(classification, pathogen, ...)

    # animal
    lookup.animal <- data.animal$lookup_tables
    data.jags$response_animal <- data.animal$response
    data.jags$associated_species <- data.animal$livestock_dat$associated_species

    # livestock
    data.jags$N_livestock <- nrow(data.animal$livestock_dat)
    data.jags$livestock_GUID <- data.animal$livestock_dat$GUID
    data.jags$livestock_sample_GUID <- data.animal$livestock_dat$sample_GUID

    # - cattle
    cattle_index <- data.animal$lookup_tables$associated_species %>%
      dplyr::filter(.data$associated_species == "cattle")
    cattle_index <- cattle_index$index
    cattle <- data.animal$livestock_dat %>%
      dplyr::filter(.data$associated_species == cattle_index)
    data.jags$N_cattle <- nrow(cattle)
    data.jags$cattle_GUID <- cattle$GUID
    data.jags$cattle_sample_GUID <- cattle$sample_GUID

    # - pigs
    pig_index <- data.animal$lookup_tables$associated_species %>%
      dplyr::filter(.data$associated_species == "pig")
    pig_index <- pig_index$index
    pig <- data.animal$livestock_dat %>%
      dplyr::filter(.data$associated_species == pig_index)
    data.jags$N_pig <- nrow(pig)
    data.jags$pig_GUID <- pig$GUID
    data.jags$pig_sample_GUID <- pig$sample_GUID

    # - chicken
    chicken_index <- data.animal$lookup_tables$associated_species %>%
      dplyr::filter(.data$associated_species == "chicken")
    chicken_index <- chicken_index$index
    chicken <- data.animal$livestock_dat %>%
      dplyr::filter(.data$associated_species == chicken_index)
    data.jags$N_chicken <- nrow(chicken)
    data.jags$chicken_GUID <- chicken$GUID
    data.jags$chicken_sample_GUID <- chicken$sample_GUID


    # companion animal
    data.jags$N_companion <- nrow(data.animal$companion_dat)
    data.jags$companion_GUID <- data.animal$companion_dat$GUID
    data.jags$companion_sample_GUID <- data.animal$companion_dat$sample_GUID

    # - cat
    cat_index <- data.animal$lookup_tables$associated_species %>%
      dplyr::filter(.data$associated_species == "cat")
    cat_index <- cat_index$index
    cat <- data.animal$companion_dat %>%
      dplyr::filter(.data$associated_species == cat_index)
    data.jags$N_cat <- nrow(cat)
    data.jags$cat_GUID <- cat$GUID
    data.jags$cat_sample_GUID <- cat$sample_GUID

    # - dog
    dog_index <- data.animal$lookup_tables$associated_species %>%
      dplyr::filter(.data$associated_species == "dog")
    dog_index <- dog_index$index
    dog <- data.animal$companion_dat %>%
      dplyr::filter(.data$associated_species == dog_index)
    data.jags$N_dog <- nrow(dog)
    data.jags$dog_GUID <- dog$GUID
    data.jags$dog_sample_GUID <- dog$sample_GUID


    # wild animal
    data.jags$N_wild <- nrow(data.animal$wild_dat)
    data.jags$wild_GUID <- data.animal$wild_dat$GUID
    data.jags$wild_sample_GUID <- data.animal$wild_dat$sample_GUID

    # - fly
    fly_index <- data.animal$lookup_tables$associated_species %>%
      dplyr::filter(.data$associated_species == "fly")
    fly_index <- fly_index$index
    fly <- data.animal$wild_dat %>%
      dplyr::filter(.data$associated_species == fly_index)
    data.jags$N_fly <- nrow(fly)
    data.jags$fly_GUID <- fly$GUID
    data.jags$fly_sample_GUID <- fly$sample_GUID

    # - turtle
    turtle_index <- data.animal$lookup_tables$associated_species %>%
      dplyr::filter(.data$associated_species == "turtle")
    turtle_index <- turtle_index$index
    turtle <- data.animal$wild_dat %>%
      dplyr::filter(.data$associated_species == turtle_index)
    data.jags$N_turtle <- nrow(turtle)
    data.jags$turtle_GUID <- turtle$GUID
    data.jags$turtle_sample_GUID <- turtle$sample_GUID

    # - crow
    crow_index <- data.animal$lookup_tables$associated_species %>%
      dplyr::filter(.data$associated_species == "crow")
    crow_index <- crow_index$index
    crow <- data.animal$wild_dat %>%
      dplyr::filter(.data$associated_species == crow_index)
    data.jags$N_crow <- nrow(crow)
    data.jags$crow_GUID <- crow$GUID
    data.jags$crow_sample_GUID <- crow$sample_GUID

    return(list(data.animal = data.animal,
                data.human = data.human,
                lookup =  c(lookup.human, lookup.animal),
                jags = data.jags))

  } else {
    lookup <- lookup.human

    return(list(data.human = data.human,
                lookup = lookup,
                jags = data.jags))
  }
}
