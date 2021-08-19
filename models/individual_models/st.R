model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability st.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(st.prob[sample_type[h_sample_GUID[p]]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(st.prob[sample_type[gp_sample_GUID[gp]]])
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(st.prob[sample_type[v_sample_GUID[v]]])
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(st.prob[sample_type[o_sample_GUID[o]]])
  }

  # ------------------------

  # Define the priors:
  # Prior distribution for st.effect (log-odds for each sample type). Sample
  # different st.effect from normal distribution for each sample type and
  # convert to a probability). Since there is only one response variable, put
  # intercept here.
  for (st in sampletypes)
  {
    st.effect[st] ~ dnorm(intercept, tau)
    logit(st.prob[st]) <- st.effect[st]
  }

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd <- sqrt(1/tau)

  #monitor# full.pd, dic, deviance, st.prob, intercept, sd
}
