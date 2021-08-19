model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability g.prob
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(g.prob[gender[h_sample_GUID[p]]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(g.prob[gender[gp_sample_GUID[gp]]])
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(g.prob[gender[v_sample_GUID[v]]])
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(g.prob[gender[o_sample_GUID[o]]])
  }

  # ------------------------

  # Define the priors:
  # Prior distribution for gender.effect (log-odds for each gender).
  gender.effect[female] ~ dnorm(0, 0.001)
  gender.effect[male] <- -gender.effect[female]

  for (g in genders)
  {
    logit(g.prob[g]) <- intercept + gender.effect[g]
  }

  # ------------------------

  # Calculate odds
  g.diff <- gender.effect[male] - gender.effect[female]
  odds.g <- exp(g.diff)

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  #monitor# full.pd, dic, deviance, intercept, gender.effect, g.prob, g.diff, odds.g, sd
}
