model {
  # Define likelihood model for data:
  for (p in 1:N_patients)
  {
    for(a in 1:antibiotic_classes)
    {
      response[h_GUID[p],a] ~ dbern(ag.prob[a,
                                            sample_month[h_sample_GUID[p]],
                                            gender[h_sample_GUID[p]]])
    }
  }

  for (gp in 1:N_gp)
  {
    for(a in 1:antibiotic_classes)
    {
      response[gp_GUID[gp],a] ~  dbern(ag.prob[a,
                                               sample_month[gp_sample_GUID[gp]],
                                               gender[gp_sample_GUID[gp]]])
    }
  }

  for (v in 1:N_volunteers)
  {
    for(a in 1:antibiotic_classes)
    {
      response[v_GUID[v],a] ~dbern(ag.prob[a,
                                           sample_month[v_sample_GUID[v]],
                                           gender[v_sample_GUID[v]]])
    }
  }

  for (o in 1:N_outpatients)
  {
    for(a in 1:antibiotic_classes)
    {
      response[o_GUID[o],a] ~ dbern(ag.prob[a,
                                            sample_month[o_sample_GUID[o]],
                                            gender[o_sample_GUID[o]]])
    }
  }

  # ------------------------

  # Define the priors:
  for(a in 1:antibiotic_classes)
  {
    antibiotic.class.effect[a] ~ dnorm(intercept, tau.class)
    logit(a.prob[a]) <- antibiotic.class.effect[a]
  }

  for(s in 1:N_sample_month)
  {
    samplemonth.effect[s] ~ dnorm(0, tau.samplemonth)
  }

  gender.effect[female] ~ dnorm(0, 0.001)
  gender.effect[male] <- -gender.effect[female]


  for(a in 1:antibiotic_classes)
  {
    for(s in 1:N_sample_month)
    {
      for (g in genders)
      {
        logit(ag.prob[a,s,g]) <- antibiotic.class.effect[a] +
          samplemonth.effect[s] + gender.effect[g]
      }
    }
  }

  # Prior value for intercept (log-odds of the average resistance in all samples)
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau.class ~ dgamma(0.001, 0.001)
  tau.samplemonth ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.class <- sqrt(1/tau.class)
  sd.samplemonth <- sqrt(1/tau.samplemonth)

  #monitor# full.pd, dic, deviance, a.prob, intercept, sd.class, sd.samplemonth
}
