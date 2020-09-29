model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability c.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients) {
    h_resist[p] ~ dbern(c.prob[clinical[sample_type[h_sample_GUID[p]]]])
    # For WAIC computation
    # h_like[p] <- logdensity.bin(h_resist[p],
    #                             c.prob[clinical[sample_type[h_sample_GUID[p]]]], 1)
  }

  for (gp in 1:N_gp) {
    gp_resist[gp] ~ dbern(gp.prob)
    # gp_like[gp] <- logdensity.bin(gp_resist[gp], gp.prob, 1)
  }

  for (v in 1:N_volunteers) {
    v_resist[v] ~ dbern(v.prob)
    # v_like[v] <- logdensity.bin(v_resist[v], v.prob, 1)
  }

  for (o in 1:N_outpatients) {
    o_resist[o] ~ dbern(o.prob)
    # o_like[o] <- logdensity.bin(o_resist[o], o.prob, 1)
  }

  # ------------------------

  # Define the priors:
  # Prior distribution for c.effect (log-odds for each clinical class). Sample
  # different c.effect from normal distribution for each clinical class and
  # convert to a probability). Since there is only one response variable, with
  # only two levels, set intercept to 0 (here we're estimating two parameters
  # -- clinical and carriage -- so we can't use more than two things to explain
  # it; infact only use tau.clin).
  for (c in c(ncarr, nclin)) # 1, 2!
  {
    c.effect[c] ~ dnorm(0, tau.clin)
    logit(c.prob[c]) <- c.effect[c]
  }

  # Prior belief about c.prob, gp.effect, v.effect, and o.effect vary
  # with c.effect
  gp.effect <- c.effect[gp_clinical]
  logit(gp.prob) <- gp.effect

  o.effect  <- c.effect[o_clinical]
  logit(o.prob) <- o.effect

  v.effect  <- c.effect[v_clinical]
  logit(v.prob) <- v.effect

  # ------------------------

  # Prior values for precision
  tau.clin ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.clin <- sqrt(1/tau.clin)

  # Calculate difference in log odds of carriage from clinical
  c.diff <- c.effect[ncarr] - c.effect[nclin]
  # Calculate odds
  odds.c <- exp(c.diff)

  #monitor# full.pd, dic, deviance, c.effect, c.prob, odds.c, sd.clin
}