#' ---
#' title: Full models
#' date: Feb 2020
#' output:
#'   html_document:
#'     theme: paper
#'     format: html_clean
#'     code_folding: hide
#'     highlight: pygments
#'     toc: TRUE
#'     toc_float: TRUE
#'     toc_collapsed: true
#' ---
#' Last updated: `r Sys.time()`

#+ setup, include=F

library(SpARKjags)
library(runjags)
library(dplyr)
set.seed(1234)

knitr::opts_chunk$set(warning = FALSE)

# Find models here (installed with the package)
directory <- "full_models"
# Save results here (use absolute path so that results can be found when
# RMarkdown report is generated)
res_dir <- file.path("", "Users", "Soniam", "Desktop", "git", "SpARK",
                     "SpARKjags", "results", directory)

data <- jags_data(classification = "all",
                  categories = "human",
                  pathogen = "Klebsiella pneumoniae",
                  removeQuinPen = T)
#+


#' ### null {.tabset}

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "null.R"),
                            save_to = res_dir)
res.null <- get_model(path)
model_name <- "res.null"

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))




#'
#' # 1. Control for the structure of the data
#'
#' ## 1a. Experimental structure
#' Find the best model with antibiotic class and sampling date
#'

#' ### res.a {.tabset}
#' response ~ antibiotic class
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "a.R"),
                            save_to = res_dir)
res.a <- get_model(path)
model_name <- "res.a"

#' #### Posterior
#+ res.a_naive, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)



#' ### res.asm {.tabset}
#' response ~ antibiotic_class + sample_month
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "asm.R"),
                            thin = 20,
                            save_to = res_dir)
res.asm <- get_model(path)
model_name <- "res.asm"

#' #### Posterior
#+ res.asm, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




#' ### res.ass {.tabset}
#' response ~ antibiotic_class + sample_season
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "ass.R"),
                            thin = 20,
                            save_to = res_dir)
res.ass <- get_model(path)
model_name <- "res.ass"

#' #### Posterior
#+ res.ass, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)


#' ### Summary of DIC results
DICtable(c("res.null", "res.a", "res.asm", "res.ass"))




#'
#' ## 1b. Demographic
#' Find the best model with gender and age
#'

#' ### res.asmg {.tabset}
#' response ~ antibiotic_class + sample_month + gender
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "asmg.R"),
                            thin = 20,
                            save_to = res_dir)
res.asmg <- get_model(path)
model_name <- "res.asmg"

#' #### Posterior
#+ res.asmg, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




#' ### res.asmag {.tabset}
#' response ~ antibiotic_class + sample_month + agegroup
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "asmag.R"),
                            thin = 20,
                            save_to = res_dir)
res.asmag <- get_model(path)
model_name <- "res.asmag"

#' #### Posterior
#+ res.asmag, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




#' ### res.asmag2 {.tabset}
#' response ~ antibiotic_class + sample_month + agegroup2
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "asmag2.R"),
                            thin = 20,
                            save_to = res_dir)
res.asmag2 <- get_model(path)
model_name <- "res.asmag2"

#' #### Posterior
#+ res.asmag2, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




#' ### res.asmage {.tabset}
#' response ~ antibiotic_class + sample_month + age
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "asmage.R"),
                            thin = 20,
                            save_to = res_dir)
res.asmage <- get_model(path)
model_name <- "res.asmage"

#' #### Posterior
#+ res.asmage, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




#' ### res.asmagesq {.tabset}
#' response ~ antibiotic_class + sample_month + age^2 + age
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory,
                                                        "asmagesq.R"),
                            thin = 20,
                            save_to = res_dir)
res.asmagesq <- get_model(path)
model_name <- "res.asmagesq"

#' #### Posterior
#+ res.asmagesq, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




#' ### res.asmagg {.tabset}
#' response ~ antibiotic_class + sample_month + agegroup + gender
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "asmagg.R"),
                            thin = 20,
                            save_to = res_dir)
res.asmagg <- get_model(path)
model_name <- "res.asmagg"

#' #### Posterior
#+ res.asmagg, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)


#' ### Summary of DIC results
DICtable(c("res.asmg", "res.asmag", "res.asmag2",
           "res.asmage", "res.asmagesq", "res.asmagg"))




#'
#' ## 1c. Hospital structure
#' Find the best model with hospital, wardtype, and ward
#'

#' ### res.asmaggh {.tabset}
#' response ~ antibiotic_class + sample_month + agegroup + gender + hospital
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "asmaggh.R"),
                            thin = 20,
                            save_to = res_dir)
res.asmaggh <- get_model(path)
model_name <- "res.asmaggh"

#' #### Posterior
#+ res.asmaggh, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




#' ### res.asmaggwt {.tabset}
#' response ~ antibiotic_class + sample_month + agegroup + gender + wardtype
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "asmaggwt.R"),
                            thin = 10,
                            save_to = res_dir)
res.asmaggwt <- get_model(path)
model_name <- "res.asmaggwt"

#' #### Posterior
#+ res.asmaggwt, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




#' ### res.asmaggw {.tabset}
#' response ~ antibiotic_class + sample_month + agegroup + gender + ward
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "asmaggw.R"),
                            thin = 10,
                            save_to = res_dir)
res.asmaggw <- get_model(path)
model_name <- "res.asmaggw"

#' #### Posterior
#+ res.asmaggw, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)



#' ### res.asmaggwt_w {.tabset}
#' response ~ antibiotic_class + sample_month + agegroup + gender +
#' wardtype_{ward}
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory,
                                                        "asmaggwt_w.R"),
                            thin = 10,
                            save_to = res_dir)
res.asmaggwt_w <- get_model(path)
model_name <- "res.asmaggwt_w"

#' #### Posterior
#+ res.asmaggwt_w, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)



#' res.asmagghwtw {.tabset}
#' response ~ antibiotic_class + sample_month + agegroup + gender + hospital +
#' wardtype + ward
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory,
                                                        "asmagghwtw.R"),
                            thin = 10,
                            save_to = res_dir)
res.asmagghwtw <- get_model(path)
model_name <- "res.asmagghwtw"

#' #### Posterior
#+ res.asmagghwtw, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




#' ### res.asmaggh_wt_w {.tabset}
#' response ~ antibiotic_class + sample_month + agegroup + gender +
#' hospital_{wardtype_{ward}}
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory,
                                                        "asmaggh_wt_w.R"),
                            thin = 10,
                            save_to = res_dir)
res.asmaggh_wt_w <- get_model(path)
model_name <- "res.asmaggh_wt_w"

#' #### Posterior
#+ res.asmaggh_wt_w, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)


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
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "asmaggw.R"),
                            thin = 10,
                            save_to = res_dir)
res.asmaggw <- get_model(path)
model_name <- "res.asmaggw"

plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)


#' ### Apply best null model to Carbapenem positive samples only
data.carb.res <- jags_data(classification = "all",
                           categories = "human",
                           pathogen = "Klebsiella pneumoniae",
                           onlyCarbapenem = T,
                           removeQuinPen = T)

#+ best null carb res
path <- run_SpARKjags_model(data.carb.res,
                            file.path(directory, "asmaggw_CarbRes.R"),
                            thin = 10)
res.asmaggw_CarbRes <- get_model(path)
model_name <- "res.asmaggw_CarbRes"

plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' ### Apply best null model to Carbapenem susceptible samples only
data.carb.sus <- jags_data(classification = "all",
                           categories = "human",
                           pathogen = "Klebsiella pneumoniae",
                           removeCarbapenem = T,
                           removeQuinPen = T)

#+ best null carb sus
path <- run_SpARKjags_model(data.carb.sus,
                            file.path(directory, "asmaggw_CarbSus.R"),
                            thin = 10)
res.asmaggw_CarbSus <- get_model(path)
model_name <- "res.asmaggw_CarbSus"

plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)




#' # 3. Check that any null components vary with antibiotic class


#'
#' # 4. Add clinical, sampletype, do they improve the model
#'
#' ### res.asmaggwc {.tabset}
#' response ~ antibioticclass + sample_month + agegroup + gender + ward +
#' clinical
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "asmaggwc.R"),
                            thin = 10,
                            save_to = res_dir)
res.asmaggwc <- get_model(path)
model_name <- "res.asmaggwc"

#' #### Posterior
#+ res.asmaggwc, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




#' ### res.asmaggwcst {.tabset}
#' response ~ antibioticclass + sample_month + agegroup + gender + ward +
#' clinical + sampletype
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory,
                                                        "asmaggwcst.R"),
                            thin = 10,
                            save_to = res_dir)
res.asmaggwcst <- get_model(path)
model_name <- "res.asmaggwcst"

#' #### Posterior
#+ res.asmaggwcst, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




#' ### res.asmaggwc_st {.tabset}
#' response ~ antibioticclass + sample_month + agegroup + gender + ward +
#' clinical_{sampletype}
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory,
                                                        "asmaggwc_st.R"),
                            thin = 10,
                            save_to = res_dir)
res.asmaggwc_st <- get_model(path)
model_name <- "res.asmaggwc_st"

#' #### Posterior
#+ res.asmaggwc_st, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir,
                 model_name = model_name)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)

#' ### Summary of DIC results

DICtable(c("res.asmg", "res.asmag", "res.asmag2", "res.asmage",
           "res.asmagesq", "res.asmagg"))
DICtable(c("res.asmaggh", "res.asmaggwt", "res.asmaggw", "res.asmaggwt_w",
           "res.asmagghwtw", "res.asmaggh_wt_w"))
DICtable(c("res.asmaggwc", "res.asmaggwcst", "res.asmaggwc_st"))




#'
#' # 4b. Add another variable with levels hospital, non-hospital
#'

#' response ~ antibioticclass + sample_month + agegroup + gender + ward + hosp + clinical_{sampletype}
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
