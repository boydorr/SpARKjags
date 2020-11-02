# SpARKjags

[![test-build](https://github.com/soniamitchell/SpARKjags/workflows/build/badge.svg)](https://github.com/soniamitchell/SpARKjags/actions)
[![codecov](https://codecov.io/gh/soniamitchell/SpARKjags/branch/master/graph/badge.svg?token=0oOkymQf8t)](https://codecov.io/gh/soniamitchell/SpARKjags)
[![CodeFactor](https://www.codefactor.io/repository/github/soniamitchell/SpARKjags/badge)](https://www.codefactor.io/repository/github/soniamitchell/SpARKjags)

## Table of contents
* [General info](#general-info)
* [Setup](#setup)
* [How to](#how-to)
  * [List all of the models included in the SpARKjags package](#list-all-of-the-models-included-in-the-sparkjags-package)
  * [Run one of the SpARKjags models](#run-one-of-the-SpARKjags-models)
  * [Run one of your own models](#run-one-of-your-own-models)
  * [Read model output into R](#read-model-output-into-R)
  * [Delete SpARKjags model output](#delete-SpARKjags-model-output)
* [Function map](#function-map)

## General info
This package contains JAGS models and associated functionality for use with SpARK project data.

## Setup
To install this package in R, run:
```R
install.packages("devtools")
devtools::install_github("soniamitchell/SpARKjags")
```
Note that you must have access to the SpARK project datasets for any of this code to work.

## How to
### List all of the models included in the SpARKjags package
```R
list_models()
```

### Run one of the SpARKjags models
```R
# Generate JAGS input dataset
data <- jags_data(classification = "Carbapenem",
                  categories = "human",
                  pathogen = "Klebsiella pneumoniae",
                  removeQuinPen = T)

# Run the model                  
path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = "individual_models/h.R",
                            save_to = "results/individual_models")
```
The `run_SpARKjags_model()` function will return the path of the model output.

    
### Read model output into R
```R
results <- get_model(path = path)
```
