library(SpARKjags)

context("Test run_model()")

# Load in data
data <- jags_data(classification = "Carbapenem",
                  categories = "human",
                  pathogen = "Klebsiella pneumoniae",
                  removeQuinPen = T)

data$jags <- data$jags[c("h_resist", "gp_resist", "v_resist", "o_resist",
                         "N_patients", "N_gp", "N_volunteers", "N_outpatients",
                         "N_hosp", "hospital", "ward", "h_sample_GUID")]

save_to <- "tmp_test_dir/h.rds"

test_that("a runjags object is generated", {
  # Run model and save output to tmp_test_dir/a.rds
    expect_warning(
      run_model(data, "individual_models/h.R",
                save_to = save_to),
      "No initial values were provided - JAGS will use the same initial values for all chains"
    )

  # Read output file
  tmp <- readRDS(save_to)

  expect_equal(class(tmp), "runjags")
})

# Delete tmp directory and contents
unlink("tmp_test_dir", recursive = TRUE)
