#' get_labels
#'
#' @param data data
#' @param ind ind
#' @param cat cat
#'
#' @export
#'
get_labels <- function(data, ind, cat) {

  tmp <- data$lookup$antibiotic_class

  tmp2 <- data.frame(index = paste0("prob.of.bad.",
                                    c("vol", "out", "gp", "hosp")),
                     stringsAsFactors = F) %>%
    mutate(category = c("Volunteer", "Outpatient", "GP", "Hospital"))

  clinical <- clinical(data, "index")
  carriage <- carriage(data, 'index')

  tmp.clin <- data.frame(index = paste0("prob.of.bad.hosp[", clinical, "]"),
                         category = "Hospital (clinical)", stringsAsFactors = F)
  tmp.car <- data.frame(index = paste0("prob.of.bad.hosp[", carriage, "]"),
                        category = "Hospital (carriage)", stringsAsFactors = F)

  tmp2 <- tmp2 %>% dplyr::bind_rows(tmp.clin, tmp.car)

  tmp.animals <- data.frame(index = c("prob.of.bad.livestock",
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

  tmp2 <- tmp2 %>% dplyr::bind_rows(tmp.animals)

  if(!missing(ind) && !missing(cat)) {
    tmp3 <- data.frame(index = ind, category = cat, stringsAsFactors = F)
    tmp2 <- tmp2 %>% dplyr::bind_rows(tmp3)
  }

  list(category = tmp2,
       good = tmp,
       bad = tmp %>% mutate(index = gsub(",1", ",2", .data$index)),
       NA, NA)
}
