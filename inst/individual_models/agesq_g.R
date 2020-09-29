model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability wc.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(gage2.prob[gender[h_sample_GUID[p]],
                                    h_sample_GUID[p]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(gage2.prob[gender[gp_sample_GUID[gp]],
                                      gp_sample_GUID[gp]])
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(gage2.prob[gender[v_sample_GUID[v]],
                                    v_sample_GUID[v]])
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(gage2.prob[gender[o_sample_GUID[o]],
                                    o_sample_GUID[o]])
  }

  # ------------------------

  # # Define the priors:
  # age.effect ~ dnorm(0, 0.001)
  # age.sq.effect ~ dnorm(0, 0.001)

  for (g in genders)
  {
    age.effect[g] ~ dnorm(0, 0.001)
    age.sq.effect[g] ~ dnorm(0, 0.001)

    # g.age.effect[g] ~ dnorm(age.effect, 0.001)
    # g.age.sq.effect[g] ~ dnorm(age.sq.effect, 0.001)

    for (s in 1:N_sample)
    {
      gage2.effect[g,s] <- intercept + age.effect[g] * age[s] +
        age.sq.effect[g] * pow(age[s], 2)
      logit(gage2.prob[g,s]) <- gage2.effect[g,s]
    }
  }

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  #monitor# full.pd, dic, deviance, gender.effect, age.effect, intercept, age.sq.effect
}
