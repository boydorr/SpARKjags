model {
  # Define likelihood model for data:
  for (p in 1:N_patients)
  {
    for(a in 1:antibiotic_classes)
    {
      response[h_GUID[p],a] ~ dbern(h.prob)
    }
  }

  for (gp in 1:N_gp)
  {
    for(a in 1:antibiotic_classes)
    {
      response[gp_GUID[gp],a] ~ dbern(gp.prob)
    }
  }

  for (v in 1:N_volunteers)
  {
    for(a in 1:antibiotic_classes)
    {
      response[v_GUID[v],a] ~ dbern(v.prob)
    }
  }

  for (o in 1:N_outpatients)
  {
    for(a in 1:antibiotic_classes)
    {
      response[o_GUID[o],a] ~ dbern(o.prob)
    }
  }

  # ------------------------

  # Define the priors:
  logit(h.prob) <- intercept
  logit(o.prob) <- intercept
  logit(gp.prob) <- intercept
  logit(v.prob) <- intercept

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, tau)

  # Prior values for precision
  tau ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd <- sqrt(1/tau)

  #monitor# full.pd, dic, deviance, h.prob, intercept, sd
}
