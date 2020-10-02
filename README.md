# SpARKjags

[![test-build](https://github.com/soniamitchell/SpARKjags/workflows/build/badge.svg)](https://github.com/soniamitchell/SpARKjags/actions)


| Function            | Description                  | Called by | Checked                | Tested                 |
| ------------------- | ---------------------------- | --------- | ---------------------- | ---------------------- |
| animal_data         | -                          | import_data | <ul><li>[ ] </li></ul> | <ul><li>[ ] </li></ul> |
| animal_lookup       |  |                                       | <ul><li>[ ] </li></ul> | <ul><li>[ ] </li></ul> |
| antibioticsplot     | Generate antibiotics heatmap | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| badgroup_posterior  | -                          | import_data | <ul><li>[ ] </li></ul> | <ul><li>[ ] </li></ul> |
| bin_ages            | -                            | get_human | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| check_fit           |  |                                       | <ul><li>[ ] </li></ul> | <ul><li>[ ] </li></ul> |
| data                | -                            | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| defineClinical      | -                            | get_human | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| densityplot         | Generate density plot        | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
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
| human_data          | -                          | import_data | <ul><li>[ ] </li></ul> | <ul><li>[ ] </li></ul> |
| import_data         | -        | antibioticsplot & densityplot & plotcorrelation & summarise_samples | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| jags_data           | Get jags data                | -         | <ul><li>[x] </li></ul> | <ul><li>[x] </li></ul> |
| monitored_variables |  | get_vars                              | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| onlyCarbapenem      | -               | get_animal & get_human | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| plot_jags           |  |                                       | <ul><li>[ ] </li></ul> | <ul><li>[ ] </li></ul> |
| plotAC              | Plot autocorrelation         | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| plotcorrelation     | Generate correlation plot    | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| removeCarbapenem    | -               | get_animal & get_human | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| run_model           | Run jags model               | -         | <ul><li>[x] </li></ul> | <ul><li>[x] </li></ul> |
| summarise_samples   | Summarise samples            | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| testPSRF            | Test PSRF                    | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| testSSEF            | Test SSEF                    | -         | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| traceplot           | Plot caterpillars |                      | <ul><li>[x] </li></ul> | <ul><li>[ ] </li></ul> |
| true_resistance     |  |                                       | <ul><li>[ ] </li></ul> | <ul><li>[ ] </li></ul> |
