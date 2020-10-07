# SpARKjags

[![test-build](https://github.com/soniamitchell/SpARKjags/workflows/build/badge.svg)](https://github.com/soniamitchell/SpARKjags/actions)
[![codecov](https://codecov.io/gh/soniamitchell/SpARKjags/branch/master/graph/badge.svg?=1)](https://codecov.io/gh/soniamitchell/SpARKjags)

## Table of contents
* [General info](#general-info)
* [Setup](#setup)
* [How to](#how-to)
* [Function map](#function-map)

## General info
This package contains JAGS models and associated functionality for use with SpARK project data.

## Setup
To install this package in R, run:
```
install.packages("devtools")
devtools::install_github("soniamitchell/SpARKjags")
```
Note that you must have access to the SpARK project datasets for any of this code to work.

## How to
* **List all of the models installed as part of this package**
  ```
  list_models()
  ```
  If you have already run any of the models, outputs will also be listed.

* **Run one of the models included in the SpARKjags package**
  ```
  # Generate JAGS input dataset
  data <- jags_data(classification = "Carbapenem",
                    categories = "human",
                    pathogen = "Klebsiella pneumoniae",
                    removeQuinPen = T)

  # Run the model                  
  path <- run_SpARKjags_model(data = data,
                              SpARKjags_model = "individual_models/h.R")
  ```
  This will save the model output to the same directory as the model script, within the SpARKjags package itself. The `run_SpARKjags_model()` function will return the path of the model output.

  * **Choose your own model output save location**
    ```
    path <- run_SpARKjags_model(data = data,
                                SpARKjags_model = "individual_models/h.R",
                                saveto = "mylocation/myresults.rds)
    ```

  * **Run your own model**
    ```
    path <- run_custom_model()
    ```
* **Read model output into R**
  ```
  results <- get_model(path)
  ```


## Function map

| Function            | Description                  | Called by | Checked                | Tested                 |
| ------------------- | ---------------------------- | --------- | ---------------------- | ---------------------- |
| animal_data         | -                          | import_data | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| animal_lookup       | -     | badgroup_posterior & import_data | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| badgroup_posterior  | -                          | import_data | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| bin_ages            | -                            | get_human | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| check_fit           |  |                                       | <ul><li>[ ] </li></ul> | <ul><li>[ ] </li></ul> |
| data                | -                            | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| defineClinical      | -                            | get_human | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| DIC                 | Get DIC                      | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| DICtable            |  |                                       | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| filter_pathogen     | -               | get_animal & get_human | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| get_animal          | -               | jags_data              | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| get_human           | -               | jags_data              | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| get_labels          | Generate labels argument |               | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| get_model           | Load model results           | -         | <ul><li>[x] </li></ul> | <ul><li>[x] </li></ul> |
| get_parameters      |  |                                       | <ul><li>[ ] </li></ul> | <ul><li>[ ] </li></ul> |
| get_params          | Generate param argument |                | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| get_vars            | Generate var.regex argument |            | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| human_data          | -                          | import_data | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| import_data         | - | plot_antibiotics, _density & _correlation & summarise_samples | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| jags_data           | Get jags data                | -         | <ul><li>[x] </li></ul> | <ul><li>[x] </li></ul> |
| monitored_variables |  | get_vars                              | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| onlyCarbapenem      | -               | get_animal & get_human | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| plot_antibiotics    | Generate antibiotics heatmap | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| plot_autocorr       | Plot autocorrelation         | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| plot_caterpillar    | Plot caterpillars |                      | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| plot_correlation    | Generate correlation plot    | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| plot_density        | Generate density plot        | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| plot_jags           |  |                                       | <ul><li>[ ] </li></ul> | <ul><li>[ ] </li></ul> |
| removeCarbapenem    | -               | get_animal & get_human | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| run_model           | Run jags model               | -         | <ul><li>[x] </li></ul> | <ul><li>[x] </li></ul> |
| summarise_samples   | Summarise samples            | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| testPSRF            | Test PSRF                    | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| testSSEF            | Test SSEF                    | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| true_resistance     |  |                                       | <ul><li>[ ] </li></ul> | <ul><li>[ ] </li></ul> |
