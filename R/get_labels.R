#' get_labels
#'
#' @export
#'
get_labels <- function(data, ind, cat) {

  tmp <- data$lookup$antibiotic_class %>%
    mutate(index = paste0("a.prob[", 1:13, ",1]"))

  tmp2 <- paste0("prob.of.bad.", c("vol", "out", "gp", "hosp")) %>%
    data.frame(index = ., stringsAsFactors = F) %>%
    mutate(category = c("Volunteer", "Outpatient", "GP", "Hospital"))

  clinical <- data$lookup$clinical %>%
    dplyr::filter(clinical == "yes") %$% index
  carriage <- data$lookup$clinical %>%
    dplyr::filter(clinical == "no") %$% index

  tmp.clin <- data.frame(index = paste0("prob.of.bad.hosp[", clinical, "]"),
                         category = "Hospital (clinical)", stringsAsFactors = F)
  tmp.car <- data.frame(index = paste0("prob.of.bad.hosp[", carriage, "]"),
                        category = "Hospital (carriage)", stringsAsFactors = F)

  tmp2 %<>% dplyr::bind_rows(tmp.clin, tmp.car)

  tmp.ainmals <- data.frame(index = c("prob.of.bad.livestock",
                                      "prob.of.bad.cattle",
                                      "prob.of.bad.pig",
                                      "prob.of.bad.chicken",
                                      "prob.of.bad.companion",
                                      "prob.of.bad.dog",
                                      "prob.of.bad.cat",
                                      "prob.of.bad.wild",
                                      "prob.of.bad.fly",
                                      "prob.of.bad.turtle",
                                      "prob.of.bad.crow"),
                            category = c("Livestock",
                                         "Cattle",
                                         "Pig",
                                         "Chicken",
                                         "Companion",
                                         "Dog",
                                         "Cat",
                                         "Wild",
                                         "Fly",
                                         "Turtle",
                                         "Crow"))

  tmp2 %<>% dplyr::bind_rows(tmp.ainmals)

  if(!missing(ind) && !missing(cat)) {
    tmp3 <- data.frame(index = ind, category = cat, stringsAsFactors = F)
    tmp2 %<>% dplyr::bind_rows(tmp3)
  }


  list(tmp2, tmp, tmp %>% mutate(index = gsub(",1", ",2", index)),
       NA, NA)
}