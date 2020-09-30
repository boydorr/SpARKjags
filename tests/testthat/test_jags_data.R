library(SpARKjags)
library(SpARK)

testthat::context("Test jags_data() can generate different datasets")

testthat::test_that("all antibiotic classes, human data, Klebsiella pneumoniae", {
  testthat::expect_silent(
    tmp <- jags_data(classification = "all",
                     categories = "human",
                     pathogen = "Klebsiella pneumoniae")
  )
})

testthat::test_that("Carbapenem, human data, Klebsiella pneumoniae", {
  testthat::expect_silent(
    jags_data(classification = "Carbapenem",
              categories = "human",
              pathogen = "Klebsiella pneumoniae")
  )
})

testthat::test_that("Carbapenem and Monobactam, human data, Klebsiella pneumoniae", {
  testthat::expect_silent(
    tmp <- jags_data(classification = c("Carbapenem", "Monobactam"),
                     categories = "human",
                     pathogen = "Klebsiella pneumoniae")
  )
})

testthat::test_that("all antibiotic classes, human and animal data, Klebsiella pneumoniae", {
  testthat::expect_silent(
    tmp <- jags_data(classification = "all",
                     categories = c("human", "animal"),
                     pathogen = "Klebsiella pneumoniae")
  )
})

testthat::test_that("Carbapenem, human and animal data, Klebsiella pneumoniae", {
  testthat::expect_silent(
    tmp <- jags_data(classification = "Carbapenem",
                     categories = c("human", "animal"),
                     pathogen = "Klebsiella pneumoniae")
  )
})

testthat::test_that("Carbapenem and Monobactam, human and animal data, Klebsiella pneumoniae", {
  testthat::expect_silent(
    tmp <- jags_data(classification = c("Carbapenem", "Monobactam"),
                     categories = c("human", "animal"),
                     pathogen = "Klebsiella pneumoniae")
  )
})

testthat::test_that("all antibiotic classes, human and animal data, Klebsiella pneumoniae and Klebsiella pasteurii", {
  testthat::expect_silent(
    tmp <- jags_data(classification = "all",
                     categories = c("human", "animal"),
                     pathogen = c("Klebsiella pneumoniae",
                                  "Klebsiella pasteurii"))
  )
})

testthat::test_that("all antibiotic classes, human and animal data, all pathogens", {
  testthat::expect_silent(
    tmp <- jags_data(classification = "all",
                     categories = c("human", "animal"),
                     pathogen = "all")
  )
})

testthat::test_that("all antibiotic classes, no Carbapenem samples", {
  testthat::expect_silent(
    tmp <- jags_data(classification = "all",
                     categories = c("human", "animal"),
                     pathogen = "Klebsiella pneumoniae",
                     removeCarbapenem = TRUE)
  )
})

testthat::test_that("all antibiotic classes, human data, only Carbapenem samples", {
  testthat::expect_silent(
    tmp <- jags_data(classification = "all",
                     categories = "human",
                     pathogen = "Klebsiella pneumoniae",
                     onlyCarbapenem = TRUE)
  )
})

testthat::test_that("all antibiotic classes, human and animal data, only Carbapenem samples", {
  testthat::expect_error(
    tmp <- jags_data(classification = "all",
                     categories = c("human", "animal"),
                     pathogen = "Klebsiella pneumoniae",
                     onlyCarbapenem = TRUE)
  )
})

testthat::test_that("all antibiotic classes, remove Quinolone and Penicillin samples", {
  testthat::expect_silent(
    tmp <- jags_data(classification = "all",
                     categories = c("human", "animal"),
                     pathogen = "Klebsiella pneumoniae",
                     removeQuinPen = TRUE))
})

testthat::test_that("all antibiotic classes, indeterminate is resistant", {
  testthat::expect_silent(
    tmp <- jags_data(classification = "all",
                     categories = c("human", "animal"),
                     pathogen = "Klebsiella pneumoniae",
                     indeterminate = "R")
  )
})

testthat::test_that("all antibiotic classes, indeterminate is susceptible", {
  testthat::expect_silent(
    tmp <- jags_data(classification = "all",
                     categories = c("human", "animal"),
                     pathogen = "Klebsiella pneumoniae",
                     indeterminate = "S")
  )
})

