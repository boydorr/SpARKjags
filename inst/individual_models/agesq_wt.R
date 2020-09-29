model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability wc.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(wtage2.prob[clinical[sample_type[h_sample_GUID[p]]],
                                    h_sample_GUID[p]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(wtage2.prob[gp_clinical,
                                      gp_sample_GUID[gp]])
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(wtage2.prob[v_clinical,
                                    v_sample_GUID[v]])
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(wtage2.prob[o_clinical,
                                    o_sample_GUID[o]])
  }

  # ------------------------

  # Define the priors:
  age.effect ~ dnorm(0, 0.001)
  age.sq.effect ~ dnorm(0, 0.001)

  for (wt in hosp_wardtypes)
  {
    wt.age.effect[wt] ~ dnorm(age.effect, tau.wt.age)
    wt.age.sq.effect[wt] ~ dnorm(age.sq.effect, tau.wt.age.sq)

    for (s in 1:N_sample)
    {
      wtage2.effect[wt,s] <- intercept + wt.age.effect[wt] * age[wt] +
        wt.age.sq.effect[wt] * pow(age[s], 2)
      logit(wtage2.prob[wt,s]) <- wtage2.effect[wt,s]
    }
  }

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  tau.wt.age ~ dgamma(0.001, 0.001)
  tau.wt.age.sq ~ dgamma(0.001, 0.001)

  sd.wt.age <- sqrt(1/tau.wt.age)
  sd.wt.age.sq <- sqrt(1/tau.wt.age.sq)

  #monitor# full.pd, dic, deviance, clin.effect, age.effect, intercept, age.sq.effect, sd.wt.age, sd.wt.age.sq
}
