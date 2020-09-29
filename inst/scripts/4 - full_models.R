#' ---
#' title: Seemingly unrelated regressions
#' date: Feb 2020
#' output:
#'   html_document:
#'     theme: paper
#'     format: html_clean
#'     code_folding: hide
#'     highlight: pygments
#' ---
#' Last updated: `r Sys.time()`

#+ r setup, include=FALSE
library(SpARK)
library(SpARKcarbapenem)
library(runjags)
library(dplyr)
set.seed(1234)

knitr::opts_chunk$set(warning = FALSE,
                      echo = FALSE)

directory <- "full_models"

data <- jags_data(classification = "all",
                  categories = "human",
                  kleb = "Klebsiella pneumoniae",
                  removeQuinPen = T)

plotnull <- function(model, data) {
  facet.titles.and.contents <- list(`probability of resistance` = "^a.prob",
                                    "intercept", "sd")
  axis.labels <- list(data$lookup$antibiotic_class %>%
                        mutate(index = paste0("a.prob[", 1:13, "]")), NA, NA)

  densityplot(model = model,
              data = data,
              var.regex = get_vars(model),
              params = facet.titles.and.contents,
              labels = axis.labels)
}

#+


#' ### null

# run_model(data, "null", directory)
res.null <- get_model("null", directory)

res.null %>% DIC() # 13030.99
res.null %>% testSSEF()
res.null %>% testPSRF()



#'
#' # 1. Control for the structure of the data
#'
#' ## 1a. Experimental structure
#' Find the best model with antibiotic class and sampling date
#'

#' ### ~ antibiotic class

# run_model(data, "a", directory)
res.a <- get_model("a", directory)

res.a %>% DIC() # 12262.57
res.a %>% testSSEF()
res.a %>% testPSRF()
res.a %>% traceplot(get_vars(.))
res.a %>% plotAC(get_vars(.))


#' ### ~ antibiotic_class + sample_month

# run_model(data, "asm", directory, thin = 20)
res.asm <- get_model("asm", directory)

res.asm %>% DIC() # 12011.55
res.asm %>% testSSEF()
res.asm %>% testPSRF()
res.asm %>% traceplot(get_vars(.))
res.asm %>% plotAC(get_vars(.))


#' ### ~ antibiotic_class + sample_season

# run_model(data, "ass", directory, thin = 20)
res.ass <- get_model("ass", directory)

res.ass %>% DIC() # 12097.86
res.ass %>% testSSEF()
res.ass %>% testPSRF()
res.ass %>% traceplot(get_vars(.))
res.ass %>% plotAC(get_vars(.))

#' ### Summary of DIC results
DICtable(c("res.null", "res.a", "res.asm", "res.ass"))



#'
#' ## 1b. Demographic
#' Find the best model with gender and age
#'

#' ### ~ antibiotic_class + sample_month + gender

# run_model(data, "asmg", directory, thin = 20)
res.asmg <- get_model("asmg", directory)

res.asmg %>% DIC() # 11862.36
res.asmg %>% testSSEF()
res.asmg %>% testPSRF()
res.asmg %>% traceplot(get_vars(.))
res.asmg %>% plotAC(get_vars(.))

#' ### ~ antibiotic_class + sample_month + agegroup

# run_model(data, "asmag", directory, thin = 20)
res.asmag <- get_model("asmag", directory)

res.asmag %>% DIC() # 11786.2
res.asmag %>% testSSEF()
res.asmag %>% testPSRF()
res.asmag %>% traceplot(get_vars(.))
res.asmag %>% plotAC(get_vars(.))

#' ### ~ antibiotic_class + sample_month + agegroup2

# run_model(data, "asmag2", directory)
res.asmag2 <- get_model("asmag2", directory)

res.asmag2 %>% DIC() # 11788.27
res.asmag2 %>% testSSEF()
res.asmag2 %>% testPSRF()
res.asmag2 %>% traceplot(get_vars(.))
res.asmag2 %>% plotAC(get_vars(.))


#' ### ~ antibiotic_class + sample_month + age

# run_model(data, "asmage", directory, thin = 20)
res.asmage <- get_model("asmage", directory)

res.asmage %>% DIC() # 11993.5
res.asmage %>% testSSEF()
res.asmage %>% testPSRF()
res.asmage %>% traceplot(get_vars(.))


#' ### ~ antibiotic_class + sample_month + age^2 + age

# run_model(data, "asmagesq", directory, thin = 20)
res.asmagesq <- get_model("asmagesq", directory)

res.asmagesq %>% DIC() # 11836.08
res.asmagesq %>% testSSEF()
res.asmagesq %>% testPSRF()
res.asmagesq %>% traceplot(get_vars(.))


#' ### ~ antibiotic_class + sample_month + agegroup + gender

# run_model(data, "asmagg", directory, thin = 20)
res.asmagg <- get_model("asmagg", directory)

res.asmagg %>% DIC() # 11647.46
res.asmagg %>% testSSEF()
res.asmagg %>% testPSRF()
res.asmagg %>% traceplot(get_vars(.))

#' ### Summary of DIC results
DICtable(c("res.asmg", "res.asmag", "res.asmag2",
           "res.asmage", "res.asmagesq", "res.asmagg"))



#'
#' ## 1c. Hospital structure
#' Find the best model with hospital, wardtype, and ward
#'

#' ### ~ antibiotic_class + sample_month + agegroup + gender + hospital

# run_model(data, "asmaggh", directory, thin = 20)
res.asmaggh <- get_model("asmaggh", directory)

res.asmaggh %>% DIC() # 10973.54
res.asmaggh %>% testSSEF()
res.asmaggh %>% testPSRF()
res.asmaggh %>% traceplot(get_vars(.))
res.asmaggh %>% plotAC(get_vars(.))


#' ### ~ antibiotic_class + sample_month + agegroup + gender + wardtype

# run_model(data, "asmaggwt", directory, thin = 10)
res.asmaggwt <- get_model("asmaggwt", directory)

res.asmaggwt %>% DIC() # 11057.68
res.asmaggwt %>% testSSEF()
res.asmaggwt %>% testPSRF()
res.asmaggwt %>% traceplot(get_vars(.))


#' ### ~ antibiotic_class + sample_month + agegroup + gender + ward

# run_model(data, "asmaggw", directory, thin = 10)
res.asmaggw <- get_model("asmaggw", directory)

res.asmaggw %>% DIC() # 9914.729
res.asmaggw %>% testSSEF()
res.asmaggw %>% testPSRF()
res.asmaggw %>% traceplot(get_vars(.))


#' ### ~ antibiotic_class + sample_month + agegroup + gender + wardtype_{ward}

#+ asmaggwt_w
# run_model(data, "asmaggwt_w", directory, thin = 10)
res.asmaggwt_w <- get_model("asmaggwt_w", directory)

res.asmaggwt_w %>% DIC() # 9915.182
res.asmaggwt_w %>% testSSEF()
res.asmaggwt_w %>% testPSRF()
res.asmaggwt_w %>% traceplot(get_vars(.))


#' ### ~ antibiotic_class + sample_month + agegroup + gender + hospital + wardtype + ward

#+ asmagghwtw
# run_model(data, "asmagghwtw", directory, thin = 10)
res.asmagghwtw <- get_model("asmagghwtw", directory)

res.asmagghwtw %>% DIC() # 9916.447
res.asmagghwtw %>% testSSEF()
res.asmagghwtw %>% testPSRF()
res.asmagghwtw %>% traceplot(get_vars(.))
res.asmagghwtw %>% plotAC(get_vars(.))


#' ### ~ antibiotic_class + sample_month + agegroup + gender + hospital_{wardtype_{ward}}

#+ asmaggh_wt_w
# run_model(data, "asmaggh_wt_w", directory, thin = 10)
res.asmaggh_wt_w <- get_model("asmaggh_wt_w", directory)

res.asmaggh_wt_w %>% DIC() # 9916.827
res.asmaggh_wt_w %>% testSSEF()
res.asmaggh_wt_w %>% testPSRF()
res.asmaggh_wt_w %>% traceplot(get_vars(.))
res.asmaggh_wt_w %>% plotAC(get_vars(.))


#' ### Summary of DIC results

DICtable(c("res.asmg", "res.asmag", "res.asmag2", "res.asmage",
           "res.asmagesq", "res.asmagg"))
DICtable(c("res.asmaggh", "res.asmaggwt", "res.asmaggw", "res.asmaggwt_w",
           "res.asmagghwtw", "res.asmaggh_wt_w"))


#'
#' # 2. Check that null components are important
#'

#' I think I already did this.
#'



#' ### This is the best null model
#+ best null
res.asmaggw <- get_model("asmaggw", directory)

plotnull(res.asmaggw, data)


#' ### Apply best null model to Carbapenem positive samples only
data.carb.res <- jags_data(classification = "all",
                           categories = "human",
                           kleb = "Klebsiella pneumoniae",
                           onlyCarbapenem = T,
                           removeQuinPen = T)
# run_model(data.carb.res, "asmaggw_CarbRes", directory, thin = 10)
#+ best null carb res
res.asmaggw_CarbRes <- get_model("asmaggw_CarbRes", directory)

plotnull(res.asmaggw_CarbRes, data)


#' ### Apply best null model to Carbapenem susceptible samples only
data.carb.sus <- jags_data(classification = "all",
                           categories = "human",
                           kleb = "Klebsiella pneumoniae",
                           removeCarbapenem = T,
                           removeQuinPen = T)

# run_model(data.carb.sus, "asmaggw_CarbSus", directory, thin = 10)
#+ best null carb sus
res.asmaggw_CarbSus <- get_model("asmaggw_CarbSus", directory)

plotnull(res.asmaggw_CarbSus, data)




#' # 3. Check that any null components vary with antibiotic class


#'
#' # 4. Add clinical, sampletype, do they improve the model
#'
#' ### ~ antibioticclass + sample_month + agegroup + gender + ward + clinical

#+ asmaggwc
# run_model(data, "asmaggwc", directory, thin = 10)
res.asmaggwc <- get_model("asmaggwc", directory)

res.asmaggwc %>% DIC() # 9911.505
res.asmaggwc %>% testSSEF()
res.asmaggwc %>% testPSRF()
res.asmaggwc %>% traceplot(get_vars(.))
res.asmaggwc %>% plotAC(get_vars(.))


#' ### ~ antibioticclass + sample_month + agegroup + gender + ward + clinical + sampletype

#+ asmaggwcst
# run_model(data, "asmaggwcst", directory, thin = 10)
res.asmaggwcst <- get_model("asmaggwcst", directory)

res.asmaggwcst %>% DIC() # 9764.723
res.asmaggwcst %>% testSSEF()
res.asmaggwcst %>% testPSRF()
res.asmaggwcst %>% traceplot(get_vars(.))
res.asmaggwcst %>% plotAC(get_vars(.))


#' ### ~ antibioticclass + sample_month + agegroup + gender + ward + clinical_{sampletype}

#+ asmaggwc_st
# run_model(data, "asmaggwc_st", directory, thin = 10)
res.asmaggwc_st <- get_model("asmaggwc_st", directory)

res.asmaggwc_st %>% DIC() # 9748.901
res.asmaggwc_st %>% testSSEF()
res.asmaggwc_st %>% testPSRF()
res.asmaggwc_st %>% traceplot(get_vars(.))
res.asmaggwc_st %>% plotAC(get_vars(.))


#' ### Summary of DIC results

DICtable(c("res.asmg", "res.asmag", "res.asmag2", "res.asmage",
           "res.asmagesq", "res.asmagg"))
DICtable(c("res.asmaggh", "res.asmaggwt", "res.asmaggw", "res.asmaggwt_w",
           "res.asmagghwtw", "res.asmaggh_wt_w"))
DICtable(c("res.asmaggwc", "res.asmaggwcst", "res.asmaggwc_st"))





#'
#' # 4b. Add another variable with levels hospital, non-hospital
#'

#' ### ~ antibioticclass + sample_month + agegroup + gender + ward + hosp + clinical_{sampletype}

# run_model(data, "asmaggwhc_st", directory, thin = 30)
# res.asmaggwhc_st <- get_model("asmaggwhc_st", directory)
#
# res.asmaggwhc_st %>% DIC() #
# res.asmaggwhc_st %>% testSSEF()
# res.asmaggwhc_st %>% testPSRF()
# res.asmaggwhc_st %>% traceplot()
# res.asmaggwhc_st %>% plotAC()

#' Can't get this to work. At thin = 20, SSeff and PRSF are terrible.
#'

#' ### antibioticclass + agegroup + gender + hosp + clinical_{sampletype}
# results.aagghc_st  <- run.jags("full_models/aagghc_st.R",
#                                 data = data.jags,
#                                 n.chains = 2,
#                                 silent.jags = T,
#                                 thin = 10)
# saveRDS(results.aagghc_st, "full_models/aagghc_st.rds")
# results.aagghc_st <- readRDS(paste0(filename, "aagghc_st.rds"))
#
# DIC(results.aagghc_st) #
# testSSEF(results.aagghc_st)
# testPSRF(results.aagghc_st)
# traceplot(results.aagghc_st)
# plotAC(results.aagghc_st)

#' Can't get this to work either. At thin = 10, SSeff and PRSF are terrible.



#'
#' # 4c. Add another variable with levels carbapenem, not-carbapenem
#'

#' #' ### antibioticclass + carbapenem + agegroup + gender + ward + clinical_{sampletype}
# results.acaggwc_st  <- run.jags("full_models/acaggwc_st.R",
#                                 data = data.jags,
#                                 n.chains = 2,
#                                 silent.jags = T,
#                                 thin = 10)
# saveRDS(results.acaggwc_st, "full_models/acaggwc_st.rds")
# results.acaggwc_st <- readRDS(paste0(filename, "acaggwc_st.rds"))
#
# DIC(results.acaggwc_st) #
# testSSEF(results.acaggwc_st)
# testPSRF(results.acaggwc_st)
# traceplot(results.acaggwc_st)
# plotAC(results.acaggwc_st)

#' Can't get this to work. There's something wrong with indexing the carbapenem effect in outpatient data.




#'
#' # 5. Do they vary with antibiotic class
#'
#
#
# model_name <- "aaggwc_st.rds"
# get_model(model_name, directory)
#
#
# labs <- cbind.data.frame(Parameter = paste0("a.prob[", 1:15, "]"),
#                          Label = colnames(data.jags$response))
# densityplot(results.aaggwc_st, labs)
#







