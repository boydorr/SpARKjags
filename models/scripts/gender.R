kleb <- "Klebsiella pneumoniae"

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


antibiotics <- unique(ATBdata$Antibiotic_name) %>% as.character()

data <- merge(METAdata, ward_types) %>%
  merge(ATBdata %>% dplyr::rename(GUID = UNIQUE_SPARK_ID))

data %<>% dplyr::mutate(Interpretation = dplyr::case_when(
  Interpretation == "I" ~ indeterminate,
  T ~ Interpretation))

data %<>%
  dplyr::filter(!grepl("NEG", GUID)) %>%
  left_join(KLEBdata %>% dplyr::select(GUID, species, ST),
            by = "GUID") %>%
  dplyr::filter(species %in% kleb,
                Clinical != "dontknow") %>%
  dplyr::rename(bacteria = species) %>%
  dplyr::select(-Phoenix_Organism)

data %<>%
  dplyr::filter(used_MIC == "yes",
                Category == "human",
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
                age = as.numeric(lubridate::interval(lubridate::ymd(DATE_OF_BIRTH),
                                                     lubridate::ymd(SAMPLE_DATE)),
                                 unit = "years"),
                sample_month = lubridate::month(lubridate::ymd(SAMPLE_DATE))) %>%
  bin_ages(10) %>%
  dplyr::select(GUID, interpretation, bacteria, ST, antibiotic,
                sample_GUID, sample_type, ward, hospital, gender,
                clinical, age_group, age_group2, ward_type,
                sample_month, age) %>%
  dplyr::group_by(GUID)

# Define clinical status
data %<>%
  dplyr::mutate(clinical = dplyr::case_when(
    hospital == -2 ~ "no",                              # volunteers
    hospital == 1 &&
      ward_type == "Sample_Collection_Center" ~ "yes",  # outpatients
    hospital == -1 ~ "yes",                             # gp
    T ~ clinical),
    sample_type = dplyr::case_when(
      hospital == -2 & sample_type == "feces" ~ "feces_volunteer",
      T ~ sample_type)) # all volunteers are fecal samples

data %<>%
  dplyr::mutate(class_interpretation = NA)

data %<>% dplyr::mutate(
  interpretation = dplyr::case_when(
    interpretation == "R" ~ 1,
    interpretation == "S" ~ 0))

data %<>%
  dplyr::ungroup() %>%
  tidyr::spread(antibiotic, interpretation)

# "Remove the sample taken from the Microbiology_and_Virology_Laboratory
# ward in San Matteo
delete_sample <- data %>% dplyr::filter(hospital != -2,
                                        ward_type == "Volunteer") %$%
  GUID %>%
  unique()
if(length(delete_sample) > 0)
  data %<>% dplyr::filter(GUID != delete_sample)

# Determine class_interpretation (resistance to each class of antibiotics)
class_tables <- ATBdata %>%
  dplyr::select(Antibiotic_name, Classification) %>%
  unique() %>%
  filter(Antibiotic_name %in% colnames(data))
all_classes <- unique(class_tables$Classification) %>% as.character()
antibiotics <- unique(class_tables$Antibiotic_name) %>% as.character()

# Generate a matrix of response variables (column for each antibiotic class)
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
  dplyr::arrange(GUID)
# %>%
# dplyr::select(-GUID) %>%
# as.matrix()
# colnames(tmp) <- NULL

w <- data %>%
  dplyr::select(hospital, ward) %>%
  unique() %>%
  dplyr::arrange(desc(hospital)) %>%
  dplyr::mutate(index = dplyr::case_when(
    hospital == -1 ~ -1,
    hospital == -2 ~ -2,
    T ~ as.numeric(dplyr::row_number()))) %>%
  dplyr::select(index, ward, -hospital)

data %<>% merge(w, by = "ward") %>%
  dplyr::select(-index)

response <- reshape2::melt(response, id.vars = "GUID",
                           variable.name = "antibiotic",
                           value.name = "resistance")

data %<>% select(-one_of(antibiotics)) %>%
  merge(response, by = "GUID")


# -------------------------------------------------------------------------


g <- data %>%
  mutate(antibiotic = as.character(antibiotic),
         antibiotic = case_when(
           antibiotic == "Trimethoprim/Sulfamethoxazole" ~ "Tri/sul",
           antibiotic == "Penicillin Combination" ~ "Pen (comb)",
           T ~ antibiotic)) %>%
  group_by(antibiotic, resistance, gender) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  ggplot2::ggplot() + ggplot2::theme_minimal() +
  ggplot2::facet_grid(resistance~antibiotic) +
  ggplot2::geom_bar(ggplot2::aes(x = gender, y = count, fill = gender),
                    colour = "black", stat = "identity") +
  ggplot2::scale_fill_manual(values = c("#d01c8b", "#2c7bb6")) +
  ggplot2::theme(strip.text.x = ggplot2::element_text(angle = 90, hjust = 0)) +
  ggplot2::labs(y = "Number of samples", x = "Gender", fill = "Gender")
ggsave("gender.pdf", g, height = 6, width = 8)


tmp <- data %>%
  mutate(antibiotic = as.character(antibiotic),
         antibiotic = case_when(
           antibiotic == "Trimethoprim/Sulfamethoxazole" ~ "Tri/sul",
           antibiotic == "Penicillin Combination" ~ "Pen (comb)",
           T ~ antibiotic)) %>%
  group_by(antibiotic, resistance, ward_type) %>%
  summarise(count = n()) %>%
  arrange(desc(count))


  ggplot2::ggplot() + ggplot2::theme_minimal() +
  ggplot2::facet_grid(resistance~antibiotic) +
  ggplot2::geom_bar(ggplot2::aes(x = ward_type, y = count, fill = ward_type),
                    colour = "black", stat = "identity") +
  # ggplot2::scale_fill_manual(values = c("#d01c8b", "#2c7bb6")) +
  ggplot2::theme(strip.text.x = ggplot2::element_text(angle = 90, hjust = 0)) +
  ggplot2::labs(y = "Number of samples", x = "ward_type", fill = "ward_type")















