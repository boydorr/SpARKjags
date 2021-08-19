model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability bc.prob
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(bc.prob[h_bacteria[p],
                                clinical[sample_type[h_sample_GUID[p]]]])
    # For WAIC computation
    # h_like[p] <- logdensity.bin(h_resist[p],
    #                             bc.prob[h_bacteria[p],
    #                                     clinical[sample_type[h_sample_GUID[p]]]], 1)
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(bc.prob[gp_bacteria[gp], gp_clinical])
    # For WAIC computation
    # gp_like[gp] <- logdensity.bin(gp_resist[gp],
    #                             bc.prob[gp_bacteria[gp], gp_clinical], 1)

  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(bc.prob[v_bacteria[v], v_clinical])
    # For WAIC computation
    # gp_like[v] <- logdensity.bin(v_resist[v],
    #                              bc.prob[v_bacteria[v], v_clinical], 1)
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(bc.prob[o_bacteria[o], o_clinical])
    # For WAIC computation
    # o_like[o] <- logdensity.bin(o_resist[o],
    #                             bc.prob[o_bacteria[o], o_clinical], 1)
  }

  # ------------------------

  # Prior distribution for mu.clin
  for (c in c(ncarr, nclin))
  {
    mu.clin[c] ~ dnorm(0, tau.clin)
  }

  # Prior distribution for b.effect
  for (b in bact_species)
  {
    b.effect[b] ~ dnorm(intercept, tau)
  }

  for (c in c(ncarr, nclin))
  {
    for (b in bact_species)
    {
      bc.effect[b,c] <- b.effect[b] + mu.clin[c]
      logit(bc.prob[b,c]) <- bc.effect[b,c]
    }
  }

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau ~ dgamma(0.001, 0.001)
  tau.clin ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd <- sqrt(1/tau)
  sd.clin <- sqrt(1/tau.clin)

  # Calculate odds
  c.diff <- mu.clin[ncarr] - mu.clin[nclin]
  odds.c <- exp(c.diff)

  #monitor# full.pd, dic, deviance, intercept, mu.clin, b.effect, odds.c, sd, sd.clin
}
