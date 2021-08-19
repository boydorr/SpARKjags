model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability wc.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(cwt.prob[ward[h_sample_GUID[p]],
                                clinical[sample_type[h_sample_GUID[p]]]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(gpc.prob)
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(vc.prob)
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(oc.prob)
  }

  # ------------------------

  # Define the priors:
  clin.effect[ncarr] ~ dnorm(0, 0.001)
  clin.effect[nclin] <- -clin.effect[ncarr]

  # Prior distribution for wt.effect (log-odds for each ward type). Sample
  # different wt.effect from normal distribution for each ward type and
  # convert to a probability). Put intercept here.
  for (wt in hosp_wardtypes)
  {
    wt.effect[wt] ~ dnorm(intercept, tau.wt)
    logit(wt.prob[wt]) <- wt.effect[wt]
  }

  # equivalent to wt.effect
  nh.effect ~ dnorm(intercept, tau.wt)
  logit(nh.prob) <- nh.effect

  # ------------------------

  # Prior distribution for wtw.effect (log-odds for each hospital ward). Sample
  # different wtw.effect from normal distribution for each hospital ward and
  # convert to a probability).
  for (w in hosp_wards)
  {
    wtw.effect[w] ~ dnorm(wt.effect[ward_type[w]], tau.w)
  }

  # equivalent to wtw.effect
  gp.effect ~ dnorm(nh.effect, tau.w)
  v.effect ~ dnorm(nh.effect, tau.w)
  o.effect ~ dnorm(nh.effect, tau.w)

  # ------------------------

  for (w in hosp_wards)
  {
    for (c in c(ncarr, nclin))
    {
      cwt.effect[w,c] <- wtw.effect[w] + clin.effect[c]
      logit(cwt.prob[w,c]) <- cwt.effect[w,c]
    }
  }

  # equivalent to cwt.effect
  gpc.effect <- gp.effect + clin.effect[gp_clinical]
  vc.effect <- v.effect + clin.effect[v_clinical]
  oc.effect <- o.effect + clin.effect[o_clinical]

  # convert to probability
  logit(gpc.prob) <- gpc.effect
  logit(vc.prob) <- vc.effect
  logit(oc.prob) <- oc.effect

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau.wt ~ dgamma(0.001, 0.001)
  tau.w ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.wt <- 1/sqrt(tau.wt)
  sd.w <- 1/sqrt(tau.w)

  # Calculate odds
  c.diff <- clin.effect[ncarr] - clin.effect[nclin]
  odds.c <- exp(c.diff)

  #monitor# full.pd, dic, deviance, intercept, clin.effect, gp.prob, v.prob, o.prob, odds.c, c.diff, sd.cwt, sd.wt, sd.w
}
