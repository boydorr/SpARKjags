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

#+ setup, include=F

library(SpARKjags)
library(runjags)
library(dplyr)
set.seed(1234)

knitr::opts_chunk$set(warning = FALSE)

directory <- "full_models"

data <- jags_data(classification = "all",
                  categories = "human",
                  pathogen = "Klebsiella pneumoniae",
                  removeQuinPen = T)

plotnull <- function(model, data) {
  facet.titles.and.contents <- list(`probability of resistance` = "^a.prob",
                                    "intercept", "sd")
  axis.labels <- list(data$lookup$antibiotic_class %>%
                        mutate(index = paste0("a.prob[", 1:13, "]")), NA, NA)

  plot_density(model = model,
               data = data,
               var.regex = get_vars(model),
               params = facet.titles.and.contents,
               labels = axis.labels)
}

#+


#' ### null {.tabset}

# run_model(data, file.path(directory, "null.R"))
res.null <- get_model(file.path("..", directory, "null.rds"))

#' #### Diagnostics
res.null %>% DIC() # 13030.99
res.null %>% testSSEF()
res.null %>% testPSRF()



#'
#' # 1. Control for the structure of the data
#'
#' ## 1a. Experimental structure
#' Find the best model with antibiotic class and sampling date
#'

#' ### ~ antibiotic class {.tabset}

# run_model(data, file.path(directory, "a.R"))
res.a <- get_model(file.path("..", directory, "a.rds"))

#' #### Diagnostics
res.a %>% DIC() # 12262.57
res.a %>% testSSEF()
res.a %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.a %>% plot_caterpillar(get_vars(.))

#' #### Autocorrelation
#+ fig.height = 6
res.a %>% plot_autocorr(get_vars(.))

#' #### Model
res.a$model

#' #### Results
res.a


#' ### ~ antibiotic_class + sample_month {.tabset}

# run_model(data, file.path(directory, "asm.R"), thin = 20)
res.asm <- get_model(file.path("..", directory, "asm.rds"))

#' #### Diagnostics
res.asm %>% DIC() # 12011.55
res.asm %>% testSSEF()
res.asm %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asm %>% plot_caterpillar(get_vars(.))

#' #### Autocorrelation
#+ fig.height = 6
res.asm %>% plot_autocorr(get_vars(.))

#' #### Model
res.asm$model

#' #### Results
res.asm


#' ### ~ antibiotic_class + sample_season {.tabset}

# run_model(data, file.path(directory, "ass.R"), thin = 20)
res.ass <- get_model(file.path("..", directory, "ass.rds"))

#' #### Diagnostics
res.ass %>% DIC() # 12097.86
res.ass %>% testSSEF()
res.ass %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.ass %>% plot_caterpillar(get_vars(.))

#' #### Autocorrelation
#+ fig.height = 6
res.ass %>% plot_autocorr(get_vars(.))

#' #### Model
res.ass$model

#' #### Results
res.ass


#' ### Summary of DIC results
DICtable(c("res.null", "res.a", "res.asm", "res.ass"))



#'
#' ## 1b. Demographic
#' Find the best model with gender and age
#'

#' ### ~ antibiotic_class + sample_month + gender {.tabset}

# run_model(data, file.path(directory, "asmg.R"), thin = 20)
res.asmg <- get_model(file.path("..", directory, "asmg.rds"))

#' #### Diagnostics
res.asmg %>% DIC() # 11862.36
res.asmg %>% testSSEF()
res.asmg %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asmg %>% plot_caterpillar(get_vars(.))

#' #### Autocorrelation
#+ fig.height = 6
res.asmg %>% plot_autocorr(get_vars(.))

#' ### ~ antibiotic_class + sample_month + agegroup {.tabset}

# run_model(data, file.path(directory, "asmag.R"), thin = 20)
res.asmag <- get_model(file.path("..", directory, "asmag.rds"))

#' #### Diagnostics
res.asmag %>% DIC() # 11786.2
res.asmag %>% testSSEF()
res.asmag %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asmag %>% plot_caterpillar(get_vars(.))

#' #### Autocorrelation
#+ fig.height = 6
res.asmag %>% plot_autocorr(get_vars(.))

#' #### Model
res.asmag$model

#' #### Results
res.asmag


#' ### ~ antibiotic_class + sample_month + agegroup2 {.tabset}

# run_model(data, file.path(directory, "asmag2.R"))
res.asmag2 <- get_model(file.path("..", directory, "asmag2.rds"))

#' #### Diagnostics
res.asmag2 %>% DIC() # 11788.27
res.asmag2 %>% testSSEF()
res.asmag2 %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asmag2 %>% plot_caterpillar(get_vars(.))

#' #### Autocorrelation
#+ fig.height = 6
res.asmag2 %>% plot_autocorr(get_vars(.))

#' #### Model
res.asmag2$model

#' #### Results
res.asmag2


#' ### ~ antibiotic_class + sample_month + age {.tabset}

# run_model(data, file.path(directory, "asmage.R"), thin = 20)
res.asmage <- get_model(file.path("..", directory, "asmage.rds"))

#' #### Diagnostics
res.asmage %>% DIC() # 11993.5
res.asmage %>% testSSEF()
res.asmage %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asmage %>% plot_caterpillar(get_vars(.))

#' #### Model
res.asmage$model

#' #### Results
res.asmage


#' ### ~ antibiotic_class + sample_month + age^2 + age {.tabset}

# run_model(data, file.path(directory, "asmagesq.R"), thin = 20)
res.asmagesq <- get_model(file.path("..", directory, "asmagesq.rds"))

#' #### Diagnostics
res.asmagesq %>% DIC() # 11836.08
res.asmagesq %>% testSSEF()
res.asmagesq %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asmagesq %>% plot_caterpillar(get_vars(.))

#' #### Model
res.asmagesq$model

#' #### Results
res.asmagesq


#' ### ~ antibiotic_class + sample_month + agegroup + gender {.tabset}

# run_model(data, file.path(directory, "asmagg.R"), thin = 20)
res.asmagg <- get_model(file.path("..", directory, "asmagg.rds"))

#' #### Diagnostics
res.asmagg %>% DIC() # 11647.46
res.asmagg %>% testSSEF()
res.asmagg %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asmagg %>% plot_caterpillar(get_vars(.))

#' #### Model
res.asmagg$model

#' #### Results
res.asmagg


#' ### Summary of DIC results
DICtable(c("res.asmg", "res.asmag", "res.asmag2",
           "res.asmage", "res.asmagesq", "res.asmagg"))



#'
#' ## 1c. Hospital structure
#' Find the best model with hospital, wardtype, and ward
#'

#' ### ~ antibiotic_class + sample_month + agegroup + gender + hospital {.tabset}

# run_model(data, file.path(directory, "asmaggh.R"), thin = 20)
res.asmaggh <- get_model(file.path("..", directory, "asmaggh.rds"))

#' #### Diagnostics
res.asmaggh %>% DIC() # 10973.54
res.asmaggh %>% testSSEF()
res.asmaggh %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asmaggh %>% plot_caterpillar(get_vars(.))

#' #### Autocorrelation
#+ fig.height = 6
res.asmaggh %>% plot_autocorr(get_vars(.))

#' #### Model
res.asmaggh$model

#' #### Results
res.asmaggh



#' ### ~ antibiotic_class + sample_month + agegroup + gender + wardtype {.tabset}

# run_model(data, file.path(directory, "asmaggwt.R"), thin = 10)
res.asmaggwt <- get_model(file.path("..", directory, "asmaggwt.rds"))

#' #### Diagnostics
res.asmaggwt %>% DIC() # 11057.68
res.asmaggwt %>% testSSEF()
res.asmaggwt %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asmaggwt %>% plot_caterpillar(get_vars(.))

#' #### Model
res.asmaggwt$model

#' #### Results
res.asmaggwt


#' ### ~ antibiotic_class + sample_month + agegroup + gender + ward {.tabset}

# run_model(data, file.path(directory, "asmaggw.R"), thin = 10)
res.asmaggw <- get_model(file.path("..", directory, "asmaggw.rds"))

#' #### Diagnostics
res.asmaggw %>% DIC() # 9914.729
res.asmaggw %>% testSSEF()
res.asmaggw %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asmaggw %>% plot_caterpillar(get_vars(.))

#' #### Model
res.asmaggw$model

#' #### Results
res.asmaggw


#' ### ~ antibiotic_class + sample_month + agegroup + gender + wardtype_{ward} {.tabset}

#+ asmaggwt_w
# run_model(data, file.path(directory, "asmaggwt_w.R"), thin = 10)
res.asmaggwt_w <- get_model(file.path("..", directory, "asmaggwt_w.rds"))

#' #### Diagnostics
res.asmaggwt_w %>% DIC() # 9915.182
res.asmaggwt_w %>% testSSEF()
res.asmaggwt_w %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asmaggwt_w %>% plot_caterpillar(get_vars(.))

#' #### Model
res.asmaggwt_w$model

#' #### Results
res.asmaggwt_w


#' ### ~ antibiotic_class + sample_month + agegroup + gender + hospital + wardtype + ward {.tabset}

#+ asmagghwtw
# run_model(data, file.path(directory, "asmagghwtw.R"), thin = 10)
res.asmagghwtw <- get_model(file.path("..", directory, "asmagghwtw.rds"))

#' #### Diagnostics
res.asmagghwtw %>% DIC() # 9916.447
res.asmagghwtw %>% testSSEF()
res.asmagghwtw %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asmagghwtw %>% plot_caterpillar(get_vars(.))

#' #### Autocorrelation
#+ fig.height = 6
res.asmagghwtw %>% plot_autocorr(get_vars(.))

#' #### Model
res.asmagghwtw$model

#' #### Results
res.asmagghwtw


#' ### ~ antibiotic_class + sample_month + agegroup + gender + hospital_{wardtype_{ward}} {.tabset}

#+ asmaggh_wt_w
# run_model(data, file.path(directory, "asmaggh_wt_w.R"), thin = 10)
res.asmaggh_wt_w <- get_model(file.path("..", directory, "asmaggh_wt_w.rds"))

#' #### Diagnostics
res.asmaggh_wt_w %>% DIC() # 9916.827
res.asmaggh_wt_w %>% testSSEF()
res.asmaggh_wt_w %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asmaggh_wt_w %>% plot_caterpillar(get_vars(.))

#' #### Autocorrelation
#+ fig.height = 6
res.asmaggh_wt_w %>% plot_autocorr(get_vars(.))

#' #### Model
res.asmaggh_wt_w$model

#' #### Results
res.asmaggh_wt_w


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
res.asmaggw <- get_model(file.path("..", directory, "asmaggw.rds"))

plotnull(res.asmaggw, data)


#' ### Apply best null model to Carbapenem positive samples only
data.carb.res <- jags_data(classification = "all",
                           categories = "human",
                           pathogen = "Klebsiella pneumoniae",
                           onlyCarbapenem = T,
                           removeQuinPen = T)
# run_model(data.carb.res, "asmaggw_CarbRes", directory, thin = 10)
#+ best null carb res
res.asmaggw_CarbRes <- get_model(file.path("..", directory, "asmaggw_CarbRes.rds"))

plotnull(res.asmaggw_CarbRes, data)


#' ### Apply best null model to Carbapenem susceptible samples only
data.carb.sus <- jags_data(classification = "all",
                           categories = "human",
                           pathogen = "Klebsiella pneumoniae",
                           removeCarbapenem = T,
                           removeQuinPen = T)

# run_model(data.carb.sus, "asmaggw_CarbSus", directory, thin = 10)
#+ best null carb sus
res.asmaggw_CarbSus <- get_model(file.path("..", directory, "asmaggw_CarbSus.rds"))

plotnull(res.asmaggw_CarbSus, data)




#' # 3. Check that any null components vary with antibiotic class


#'
#' # 4. Add clinical, sampletype, do they improve the model
#'
#' ### ~ antibioticclass + sample_month + agegroup + gender + ward + clinical {.tabset}

#+ asmaggwc
# run_model(data, file.path(directory, "asmaggwc.R"), thin = 10)
res.asmaggwc <- get_model(file.path("..", directory, "asmaggwc.rds"))

#' #### Diagnostics
res.asmaggwc %>% DIC() # 9911.505
res.asmaggwc %>% testSSEF()
res.asmaggwc %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asmaggwc %>% plot_caterpillar(get_vars(.))

#' #### Autocorrelation
#+ fig.height = 6
res.asmaggwc %>% plot_autocorr(get_vars(.))

#' #### Model
res.asmaggwc$model

#' #### Results
res.asmaggwc


#' ### ~ antibioticclass + sample_month + agegroup + gender + ward + clinical + sampletype {.tabset}

#+ asmaggwcst
# run_model(data, file.path(directory, "asmaggwcst.R"), thin = 10)
res.asmaggwcst <- get_model(file.path("..", directory, "asmaggwcst.rds"))

#' #### Diagnostics
res.asmaggwcst %>% DIC() # 9764.723
res.asmaggwcst %>% testSSEF()
res.asmaggwcst %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asmaggwcst %>% plot_caterpillar(get_vars(.))

#' #### Autocorrelation
#+ fig.height = 6
res.asmaggwcst %>% plot_autocorr(get_vars(.))

#' #### Model
res.asmaggwcst$model

#' #### Results
res.asmaggwcst


#' ### ~ antibioticclass + sample_month + agegroup + gender + ward + clinical_{sampletype} {.tabset}

#+ asmaggwc_st
# run_model(data, file.path(directory, "asmaggwc_st.R"), thin = 10)
res.asmaggwc_st <- get_model(file.path("..", directory, "asmaggwc_st.rds"))

#' #### Diagnostics
res.asmaggwc_st %>% DIC() # 9748.901
res.asmaggwc_st %>% testSSEF()
res.asmaggwc_st %>% testPSRF()

#' #### plot_caterpillar
#+ fig.height = 6
res.asmaggwc_st %>% plot_caterpillar(get_vars(.))

#' #### Autocorrelation
#+ fig.height = 6
res.asmaggwc_st %>% plot_autocorr(get_vars(.))

#' #### Model
res.asmaggwc_st$model

#' #### Results
res.asmaggwc_st

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
# res.asmaggwhc_st
#' Can't get this to work. At thin = 20, SSeff and PRSF are terrible.
#'

#' ### antibioticclass + agegroup + gender + hosp + clinical_{sampletype}
# results.aagghc_st

#' Can't get this to work either. At thin = 10, SSeff and PRSF are terrible.

#'
#' # 4c. Add another variable with levels carbapenem, not-carbapenem
#'

#' #' ### antibioticclass + carbapenem + agegroup + gender + ward + clinical_{sampletype}
# results.acaggwc_st
#' Can't get this to work. There's something wrong with indexing the carbapenem effect in outpatient data.

#'
#' # 5. Do they vary with antibiotic class
#'
# res.aaggwc_st
