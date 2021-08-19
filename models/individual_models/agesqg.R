model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability wc.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(genage.prob[gender[h_sample_GUID[p]],
                                    h_sample_GUID[p]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(genage.prob[gender[gp_sample_GUID[gp]],
                                             gp_sample_GUID[gp]])
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(genage.prob[gender[v_sample_GUID[v]],
                                           v_sample_GUID[v]])
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(genage.prob[gender[o_sample_GUID[o]],
                                           o_sample_GUID[o]])
  }

  # ------------------------

  # Define the priors:

  # Prior distribution for w.effect
  gender.effect[female] ~ dnorm(0, 0.001)
  gender.effect[male] <- -gender.effect[female]

  for (s in 1:N_sample)
  {
    mu.age[s] <- age.effect * age[s] + age.sq.effect * pow(age[s], 2)
  }

  for (g in genders)
  {
    for (s in 1:N_sample)
    {
      logit(genage.prob[g,s]) <- gender.effect[g] + mu.age[s] + intercept
    }
  }

  # ------------------------

  age.effect ~ dnorm(0, 0.001)
  age.sq.effect ~ dnorm(0, 0.001)

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  #monitor# full.pd, dic, deviance, gender.effect, age.effect, intercept, age.sq.effect
}
