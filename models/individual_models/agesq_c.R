model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability wc.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(cage2.prob[clinical[sample_type[h_sample_GUID[p]]],
                                   h_sample_GUID[p]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(cage2.prob[gp_clinical,
                                     gp_sample_GUID[gp]])
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(cage2.prob[v_clinical,
                                   v_sample_GUID[v]])
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(cage2.prob[o_clinical,
                                   o_sample_GUID[o]])
  }

  # ------------------------

  # Define the priors:


  for (c in c(ncarr, nclin))
  {
    age.effect[c] ~ dnorm(0, 0.001)
    age.sq.effect[c] ~ dnorm(0, 0.001)

    # c.age.effect[c] ~ dnorm(age.effect, tau.c.age)
    # c.age.sq.effect[c] ~ dnorm(age.sq.effect, tau.c.age.sq)

    for (s in 1:N_sample)
    {
      cage2.effect[c,s] <- intercept + age.effect[c] * age[c] +
        age.sq.effect[c] * pow(age[s], 2)
      logit(cage2.prob[c,s]) <- cage2.effect[c,s]
    }
  }

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  #monitor# full.pd, dic, deviance, clin.effect, age.effect, intercept, age.sq.effect
}
