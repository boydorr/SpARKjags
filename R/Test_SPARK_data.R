#Test data for Sonia
# library(dplyr)
# library(ggplot2)
# library(codetools)

################################################################################

#' add_resistance
#'
#' @param df argument
#' @param N argument
#' @param ac.prob.list argument
#'
add_resistance <- function(df, N, ac.prob.list){
  ab_vec <- c()
  for (ab in c(1,2,3)){
    for (i in 1:N){
      if (df$clinical[i] == 1 & df$bad[i] == 1) {
        ab_vec[i] <- stats::rbinom(1,1, ac.prob.list[["clinical"]][ab, "bad"])
      }
      if (df$clinical[i] == 1 & df$bad[i] == 0) {
        ab_vec[i] <- stats::rbinom(1,1, ac.prob.list[["clinical"]][ab, "good"])
      }
      if (df$clinical[i] == 0 & df$bad[i] == 1) {
        ab_vec[i] <- stats::rbinom(1,1, ac.prob.list[["carriage"]][ab, "bad"])
      }
      if (df$clinical[i] == 0 & df$bad[i] == 0) {
        ab_vec[i] <- stats::rbinom(1,1, ac.prob.list[["carriage"]][ab, "good"])
      }
    }
    df <- cbind(df, ab_vec)
  }
  colnames(df) <- c("id", "source", "clinical", "bad", "ab1", "ab2", "ab3")
  df
}




#' add_bad_column_hosp
#'
#' @param df argument
#' @param N argument
#' @param prob.of.bad.hosp argument
#'
add_bad_column_hosp <- function(df, N, prob.of.bad.hosp){
  bad <- c()
  for (i in 1:N){
    if (df$clinical[i] == 1){
      bad[i] <- stats::rbinom(1, 1, prob.of.bad.hosp[2])
    }
    if (df$clinical[i] == 0){
      bad[i] <- stats::rbinom(1, 1, prob.of.bad.hosp[1])
    }
  }
  df$bad <- bad
  df
}

add_bad_column_other <- function(df, N, prob.of.bad){
  bad <- c()
  for (i in 1:N){
    bad[i] <- stats::rbinom(1, 1, prob.of.bad)
  }
  df$bad <- bad
  df
}




#' get_test_data
#'
#' @export
#'
get_test_data <- function() {

  ################################################################################
  #Set parameters
  intercept <- -1.5
  diff <- 2 #difference between good and bad
  intercept.plus <- intercept + diff

  sd.class <- 0.5
  sd.clin <- 0.1

  prob.of.bad.hosp <- c(0.6, 0.4) #Carriage, Clinical
  prob.of.bad.gp <- 0.2
  prob.of.bad.vol <- 0.05
  prob.of.bad.out <- 0.3

  ################################################################################
  ac.prob.list <- list()
  antibiotic_classes <- 3
  for (c in c(1,2)){

    my_dimnames = list(c("ab1", "ab2", "ab3"), c("good", "bad"))
    antibiotic.class.effect <- matrix(nrow = antibiotic_classes, ncol = 2,
                                      dimnames = my_dimnames)
    ac.effect <- matrix(nrow = antibiotic_classes, ncol = 2,
                        dimnames = my_dimnames)
    ac.prob <- matrix(nrow = antibiotic_classes, ncol = 2,
                      dimnames = my_dimnames)
    for (a in 1:antibiotic_classes){

      antibiotic.class.effect[a, 1] <- stats::rnorm(1,intercept,
                                                    sd.class) #good group
      antibiotic.class.effect[a, 2] <- stats::rnorm(1,intercept.plus,
                                                    sd.class) #bad group

      ac.effect[a,1] <- stats::rnorm(1,antibiotic.class.effect[a, 1], sd.clin)
      ac.effect[a,2] <- stats::rnorm(1,antibiotic.class.effect[a, 2], sd.clin)
      ac.prob[a,1] <- exp(ac.effect[a,1])/(1+exp(ac.effect[a,1]))
      ac.prob[a,2] <- exp(ac.effect[a,2])/(1+exp(ac.effect[a,2]))

      ac.prob.list[[c]] <- ac.prob
    }
  }
  names(ac.prob.list) <- c("carriage", "clinical")
  ac.prob.list

  ################################################################################
  #Create test data
  #hospital
  hospital_data <- data.frame()
  N_patients <- 874
  proportion_clinical <- 0.8

  hospital_data <- data.frame(id = seq(1:N_patients))
  hospital_data$source <- "hospital"
  hospital_data$clinical <- stats::rbinom(N_patients, 1, proportion_clinical)
  hospital_data <- add_bad_column_hosp(hospital_data, N_patients, prob.of.bad.hosp)
  hospital_data <- add_resistance(hospital_data, N_patients, ac.prob.list)
  #remove bad column
  hospital_data <- hospital_data %>% dplyr::select(-bad)
  # hospital_data

  #gp
  gp_data <- data.frame()
  N_gp <- 17
  proportion_clinical <- 0.0

  gp_data <- data.frame(id = seq(1:N_gp))
  gp_data$source <- "gp"
  gp_data$clinical <- stats::rbinom(N_gp, 1, proportion_clinical)
  gp_data <- add_bad_column_other(gp_data, N_gp, prob.of.bad.gp)
  gp_data <- add_resistance(gp_data, N_gp, ac.prob.list)
  #remove bad column
  gp_data <- gp_data %>% dplyr::select(-bad)
  # gp_data

  #out
  out_data <- data.frame()
  N_out <- 106
  proportion_clinical <- 1.0

  out_data <- data.frame(id = seq(1:N_out))
  out_data$source <- "out"
  out_data$clinical <- stats::rbinom(N_out, 1, proportion_clinical)
  out_data <- add_bad_column_other(out_data, N_out, prob.of.bad.out)
  out_data <- add_resistance(out_data, N_out, ac.prob.list)
  #remove bad column
  out_data <- out_data %>% dplyr::select(-bad)
  # out_data

  #vol
  vol_data <- data.frame()
  N_vol <- 22
  proportion_clinical <- 0.0

  vol_data <- data.frame(id = seq(1:N_vol))
  vol_data$source <- "vol"
  vol_data$clinical <- stats::rbinom(N_vol, 1, proportion_clinical)
  vol_data <- add_bad_column_other(vol_data, N_vol, prob.of.bad.vol)
  vol_data <- add_resistance(vol_data, N_vol, ac.prob.list)
  #remove bad column
  vol_data <- vol_data %>% dplyr::select(-bad)
  # vol_data

  # head(hospital_data)
  # head(gp_data)
  # head(out_data)
  # head(vol_data)

  #Combine dataframes
  all_test_data <- rbind(hospital_data, gp_data, out_data, vol_data)
  #Remove group specific id and replace with unique id
  all_test_data <- all_test_data %>% dplyr::select(-id)
  all_test_data$unique_id <- seq(1:nrow(all_test_data))
  # View(all_test_data)


  ################################################################################
  # findGlobals(add_resistance, merge=FALSE)$variables
  # findGlobals(add_bad_column_other, merge=FALSE)$variables
  # findGlobals(add_bad_column_hosp, merge=FALSE)$variables


  # -------------------------------------------------------------------------

  data <- list()
  data$data.human <- list()
  data$lookup <- list()
  data$jags <- list()

  data$data.human$hospital_dat <- all_test_data  %>%
    dplyr::filter(source == "hospital") %>%
    dplyr::select(.data$unique_id, .data$ab1, .data$ab2, .data$ab3) %>%
    dplyr::rename(GUID = .data$unique_id) %>%
    dplyr::mutate(sample_GUID = .data$GUID,
           bacteria = 1,
           ST = NA,
           class_interpretation = NA) %>%
    dplyr::rename(Aminoglycoside = .data$ab1,
           Carbapenem = .data$ab2,
           Cephalosporin = .data$ab3)

  data$data.human$gp_dat <- all_test_data  %>%
    dplyr::filter(source == "gp") %>%
    dplyr::select(.data$unique_id, .data$ab1, .data$ab2, .data$ab3) %>%
    dplyr::rename(GUID = .data$unique_id) %>%
    dplyr::mutate(sample_GUID = .data$GUID,
           bacteria = 1,
           ST = NA,
           class_interpretation = NA) %>%
    dplyr::rename(Aminoglycoside = .data$ab1,
           Carbapenem = .data$ab2,
           Cephalosporin = .data$ab3)

  data$data.human$outpatients_dat <- all_test_data  %>%
    dplyr::filter(source == "out") %>%
    dplyr::select(.data$unique_id, .data$ab1, .data$ab2, .data$ab3) %>%
    dplyr::rename(GUID = .data$unique_id) %>%
    dplyr::mutate(sample_GUID = .data$GUID,
           bacteria = 1,
           ST = NA,
           class_interpretation = NA) %>%
    dplyr::rename(Aminoglycoside = .data$ab1,
           Carbapenem = .data$ab2,
           Cephalosporin = .data$ab3)

  data$data.human$volunteer_dat <- all_test_data  %>%
    dplyr::filter(source == "vol") %>%
    dplyr::select(.data$unique_id, .data$ab1, .data$ab2, .data$ab3) %>%
    dplyr::rename(GUID = .data$unique_id) %>%
    dplyr::mutate(sample_GUID = .data$GUID,
           bacteria = 1,
           ST = NA,
           class_interpretation = NA) %>%
    dplyr::rename(Aminoglycoside = .data$ab1,
           Carbapenem = .data$ab2,
           Cephalosporin = .data$ab3)

  data$data.human$response <- all_test_data %>%
    dplyr::select(.data$ab1, .data$ab2, .data$ab3) %>%
    dplyr::rename(Aminoglycoside = .data$ab1,
           Carbapenem = .data$ab2,
           Cephalosporin = .data$ab3)

  data$data.human$data <- all_test_data %>%
    dplyr::mutate(GUID = paste0("SPARK_", .data$unique_id),
           ward = NA,
           bacteria = NA,
           ST = NA,
           sample_GUID = .data$GUID,
           sample_type = NA,
           hospital = NA,
           gender = NA,
           clinical = dplyr::case_when(clinical == 1 ~ "yes",
                                       clinical == 0 ~ "no"),
           age_group = NA,
           age_group2 = NA,
           ward_type = NA,
           sample_month = NA,
           sample_season = NA,
           age = NA,
           class_interpretation = NA,
           dataset = NA) %>%
    dplyr::rename(index = .data$unique_id,
           Aminoglycoside = .data$ab1,
           Carbapenem = .data$ab2,
           Cephalosporin = .data$ab3)

  data$data.human$lookup_tables$GUID <- all_test_data %>%
    dplyr::rename(GUID = .data$unique_id) %>%
    dplyr::mutate(GUID =  paste0("SPARK_", .data$GUID),
           index = 1:nrow(all_test_data)) %>%
    dplyr::rename(Aminoglycoside = .data$ab1,
           Carbapenem = .data$ab2,
           Cephalosporin = .data$ab3) %>%
    dplyr::mutate(source = dplyr::case_when(source == "hospital" ~ "Hospital",
                                     source == "gp" ~ "GP",
                                     source == "out" ~ "Outpatients",
                                     source == "vol" ~ "Volunteers"),
           clinical = dplyr::case_when(clinical == 1 ~ "yes",
                                       clinical == 0 ~ "no"),
           hospital = NA)

  data$lookup$clinical <- data.frame(clinical = c("no", "yes"), index = 1:2)
  data$lookup$antibiotic_class <- data.frame(
    index = 1:3,
    classification = c("Aminoglycoside", "Carbapenem", "Cephalosporin"))

  data$jags$antibiotic_classes <- 3

  data$jags$hospital_response <- all_test_data  %>%
    dplyr::filter(source == "hospital") %>%
    dplyr::select(.data$ab1, .data$ab2, .data$ab3) %>%
    as.matrix()

  data$jags$N_patients <- all_test_data %>%
    dplyr::filter(source == "hospital") %>%
    nrow()

  data$jags$hospital_clinical <- all_test_data %>%
    dplyr::filter(source == "hospital") %>%
    dplyr::mutate(clinical = as.numeric(as.factor(clinical))) %>%
    dplyr::select(clinical) %>%
    unlist()

  data$jags$ncarr <- 1
  data$jags$nclin <- 2
  data$jags$o_clinical <- 2
  data$jags$gp_clinical <- 2
  data$jags$v_clinical <- 1

  data$jags$gp_response <- all_test_data  %>%
    dplyr::filter(source == "gp") %>%
    dplyr::select(.data$ab1, .data$ab2, .data$ab3) %>%
    as.matrix()

  data$jags$N_gp <- all_test_data %>%
    dplyr::filter(source == "gp") %>%
    nrow()

  data$jags$vol_response <- all_test_data  %>%
    dplyr::filter(source == "vol") %>%
    dplyr::select(.data$ab1, .data$ab2, .data$ab3) %>%
    as.matrix()

  data$jags$N_volunteers <- all_test_data %>%
    dplyr::filter(source == "vol") %>%
    nrow()

  data$jags$out_response <- all_test_data  %>%
    dplyr::filter(source == "out") %>%
    dplyr::select(.data$ab1, .data$ab2, .data$ab3) %>%
    as.matrix()

  data$jags$N_outpatients <- all_test_data %>%
    dplyr::filter(source == "out") %>%
    nrow()

  data
}

