model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability w.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(st.prob[hospital[ward[h_sample_GUID[p]]],
                                sample_type[h_sample_GUID[p]]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(gp.prob)
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(v.prob)
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(o.prob)
  }

  # ------------------------

  # Prior distribution for h.effect (log-odds for each hospital). Sample
  # different h.effect from normal distribution for each hospital and
  # convert to a probability).
  for (h in 1:N_hosp)
  {
    h.effect[h] ~ dnorm(intercept, tau.hosp)
    logit(h.prob[h]) <- h.effect[h]
  }

  # equivalent to h.effect
  nh.effect ~ dnorm(intercept, tau.hosp)
  logit(nh.prob) <- nh.effect

  # ------------------------

  # Define the priors:
  # Prior distribution for st.effect (log-odds for each sample type). Sample
  # different st.effect from normal distribution for each sample type and
  # convert to a probability).
  for (h in 1:N_hosp)
  {
    for (st in sampletypes)
    {
      st.effect[h,st] ~ dnorm(h.effect[h], tau.st)
      logit(st.prob[h,st]) <- st.effect[h,st]
    }
  }

  # equivalent to st.effect
  gp.effect ~ dnorm(nh.effect, tau.st)
  logit(gp.prob) <- gp.effect

  v.effect ~ dnorm(nh.effect, tau.st)
  logit(v.prob) <- v.effect

  o.effect ~ dnorm(nh.effect, tau.st)
  logit(o.prob) <- o.effect

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau.hosp ~ dgamma(0.001, 0.001)
  tau.st ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.hosp <- sqrt(1/tau.hosp)
  sd.st <- sqrt(1/tau.st)

  #monitor# full.pd, dic, deviance, intercept, h.prob, nh.prob, sd.hosp, sd.st
}
