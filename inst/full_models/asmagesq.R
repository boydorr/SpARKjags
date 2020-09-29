model {
  # Define likelihood model for data:
  for (p in 1:N_patients)
  {
    for(a in 1:antibiotic_classes)
    {
      response[h_GUID[p],a] ~ dbern(aag.prob[a,
                                             sample_month[h_sample_GUID[p]],
                                             h_sample_GUID[p]])
    }
  }

  for (gp in 1:N_gp)
  {
    for(a in 1:antibiotic_classes)
    {
      response[gp_GUID[gp],a] ~ dbern(aag.prob[a,
                                               sample_month[gp_sample_GUID[gp]],
                                               gp_sample_GUID[gp]])
    }
  }

  for (v in 1:N_volunteers)
  {
    for(a in 1:antibiotic_classes)
    {
      response[v_GUID[v],a] ~ dbern(aag.prob[a,
                                             sample_month[v_sample_GUID[v]],
                                             v_sample_GUID[v]])
    }
  }

  for (o in 1:N_outpatients)
  {
    for(a in 1:antibiotic_classes)
    {
      response[o_GUID[o],a] ~ dbern(aag.prob[a,
                                             sample_month[o_sample_GUID[o]],
                                             o_sample_GUID[o]])
    }
  }

  # ------------------------

  # Define the priors:
  for(a in 1:antibiotic_classes)
  {
    antibiotic.class.effect[a] ~ dnorm(intercept, tau.class)
    logit(a.prob[a]) <- antibiotic.class.effect[a]
  }

  for(m in 1:N_sample_month)
  {
    samplemonth.effect[m] ~ dnorm(0, tau.samplemonth)
  }

  age.effect ~ dnorm(0, 0.001)
  age2.effect ~ dnorm(0, 0.001)

  for (s in 1:N_sample)
  {
    mu.age[s] <- age.effect * age[s] + age2.effect * pow(age[s], 2)
  }

  for(a in 1:antibiotic_classes)
  {
    for(m in 1:N_sample_month)
    {
      for (s in 1:N_sample)
      {
        logit(aag.prob[a,m,s]) <- antibiotic.class.effect[a] +
          samplemonth.effect[m] + mu.age[s]
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
