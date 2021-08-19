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
      response[h_GUID[p],a] ~ dbern(a.prob[a, index.bad.p[p]])
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
      response[gp_GUID[gp],a] ~ dbern(a.prob[a, index.bad.gp[gp]])
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
      response[v_GUID[v],a] ~ dbern(a.prob[a, index.bad.v[v]])
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
      response[o_GUID[o],a] ~ dbern(a.prob[a, index.bad.o[o]])
    }
  }

  for (x in 1:N_cattle)
  {
    # probability of belonging to the bad group
    bad.cattle[x] ~ dbern(prob.of.bad.cattle)
    # index the bad group in a.prob
    index.bad.cattle[x] <- bad.cattle[x] + 1

    for(a in 1:antibiotic_classes)
    {
      # Response is different for each antibiotic and depending on
      # which pop it's from
      response_animal[cattle_GUID[x],a] ~ dbern(a.prob[a, index.bad.cattle[x]])
    }
  }

  for (y in 1:N_pig)
  {
    # probability of belonging to the bad group
    bad.pig[y] ~ dbern(prob.of.bad.pig)
    # index the bad group in a.prob
    index.bad.pig[y] <- bad.pig[y] + 1

    for(a in 1:antibiotic_classes)
    {
      # Response is different for each antibiotic and depending on
      # which pop it's from
      response_animal[pig_GUID[y],a] ~ dbern(a.prob[a, index.bad.pig[y]])
    }
  }

  for (z in 1:N_chicken)
  {
    # probability of belonging to the bad group
    bad.chicken[z] ~ dbern(prob.of.bad.chicken)
    # index the bad group in a.prob
    index.bad.chicken[z] <- bad.chicken[z] + 1

    for(a in 1:antibiotic_classes)
    {
      # Response is different for each antibiotic and depending on
      # which pop it's from
      response_animal[chicken_GUID[z],a] ~ dbern(a.prob[a, index.bad.chicken[z]])
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

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  # Difference between distribution means of "good" and "bad" groups
  diff ~ dgamma(0.001, 0.001)
  intercept.plus <- intercept + diff

  # Probability of being in the bad group
  prob.of.bad.hosp ~ dbeta(1, 1)
  prob.of.bad.gp ~ dbeta(1, 1)
  prob.of.bad.vol ~ dbeta(1, 1)
  prob.of.bad.out ~ dbeta(1, 1)
  prob.of.bad.cattle ~ dbeta(1, 1)
  prob.of.bad.pig ~ dbeta(1, 1)
  prob.of.bad.chicken ~ dbeta(1, 1)

  # Prior values for precision
  tau.class ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.class <- sqrt(1/tau.class)

  #monitor# full.pd, dic, deviance, a.prob, prob.of.bad.hosp, prob.of.bad.gp, prob.of.bad.vol, prob.of.bad.out, prob.of.bad.cattle, prob.of.bad.pig, prob.of.bad.chicken, bad.p, bad.gp, bad.v, bad.o, bad.cattle, bad.pig, bad.chicken, intercept, sd.class
}
