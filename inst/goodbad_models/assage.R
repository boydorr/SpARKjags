model {
  # Define likelihood model for data:
  for (p in 1:N_patients)
  {
    # probability of belonging to the bad group
    bad.p[p] ~ dbern(prob.of.bad.hosp)
    # index the bad group in a.prob
    index.bad.p[p] <- bad.p[p] + 1

    for(a in 1:antibiotic_classes)
    {
      # Response is different for each antibiotic and depending on
      # which pop it's from
      response[h_GUID[p],a] ~ dbern(aag.prob[a,
                                             index.bad.p[p],
                                             sample_season[h_sample_GUID[p]],
                                             h_sample_GUID[p]])
    }
  }

  for (gp in 1:N_gp)
  {
    # probability of belonging to the bad group
    bad.gp[gp] ~ dbern(prob.of.bad.gp)
    # index the bad group in a.prob
    index.bad.gp[gp] <- bad.gp[gp] + 1

    for(a in 1:antibiotic_classes)
    {
      # Response is different for each antibiotic and depending on
      # which pop it's from
      response[gp_GUID[gp],a] ~ dbern(aag.prob[a,
                                               index.bad.gp[gp],
                                               sample_season[gp_sample_GUID[gp]],
                                               gp_sample_GUID[gp]])
    }
  }

  for (v in 1:N_volunteers)
  {
    # probability of belonging to the bad group
    bad.v[v] ~ dbern(prob.of.bad.vol)
    # index the bad group in a.prob
    index.bad.v[v] <- bad.v[v] + 1

    for(a in 1:antibiotic_classes)
    {
      # Response is different for each antibiotic and depending on
      # which pop it's from
      response[v_GUID[v],a] ~ dbern(aag.prob[a,
                                             index.bad.v[v],
                                             sample_season[v_sample_GUID[v]],
                                             v_sample_GUID[v]])
    }
  }

  for (o in 1:N_outpatients)
  {
    # probability of belonging to the bad group
    bad.o[o] ~ dbern(prob.of.bad.out)
    # index the bad group in a.prob
    index.bad.o[o] <- bad.o[o] + 1

    for(a in 1:antibiotic_classes)
    {
      # Response is different for each antibiotic and depending on
      # which pop it's from
      response[o_GUID[o],a] ~ dbern(aag.prob[a,
                                             index.bad.o[o],
                                             sample_season[o_sample_GUID[o]],
                                             o_sample_GUID[o]])
    }
  }

  # ------------------------

  # Define the priors:
  for(a in 1:antibiotic_classes)
  {
    # probability of being resistant in the good group (less resistances)
    antibiotic.class.effect[a, 1] ~ dnorm(intercept, tau.class)
    logit(a.prob[a,1]) <- antibiotic.class.effect[a, 1]

    # probability of being resistant in the bad group (many resistances) will
    # always be higher than the good group
    antibiotic.class.effect[a, 2] ~ dnorm(intercept.plus, tau.class)
    logit(a.prob[a,2]) <- antibiotic.class.effect[a, 2]
  }

  for(n in 1:N_sample_season)
  {
    sampleseason.effect[n] ~ dnorm(0, tau.sampleseason)
  }

  age.effect ~ dnorm(0, 0.001)

  for (s in 1:N_sample)
  {
    mu.age[s] <- age.effect * age[s]
  }

  for(a in 1:antibiotic_classes)
  {
    for(b in 1:2)
    {
      for(n in 1:N_sample_season)
      {
        for (s in 1:N_sample)
        {
          logit(aag.prob[a,b,n,s]) <- antibiotic.class.effect[a,b] +
            sampleseason.effect[n] + mu.age[s]
        }
      }
    }
  }

  # Prior value for intercept (log-odds of the average resistance in all samples)
  intercept ~ dnorm(0, 0.001)

  # Difference between distribution means of "good" and "bad" groups
  diff ~ dgamma(0.001, 0.001)
  intercept.plus <- intercept + diff

  # Probability of being in the bad group
  prob.of.bad.hosp ~ dbeta(1, 1)
  prob.of.bad.gp ~ dbeta(1, 1)
  prob.of.bad.vol ~ dbeta(1, 1)
  prob.of.bad.out ~ dbeta(1, 1)

  # Prior values for precision
  tau.class ~ dgamma(0.001, 0.001)
  tau.sampleseason ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.class <- sqrt(1/tau.class)
  sd.sampleseason <- sqrt(1/tau.sampleseason)

  #monitor# full.pd, dic, deviance, a.prob, prob.of.bad.hosp, prob.of.bad.gp, prob.of.bad.vol, prob.of.bad.out, bad.p, bad.gp, bad.v, bad.o, intercept, sd.class, sd.sampleseason
}
