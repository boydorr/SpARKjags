model {
  # Define likelihood model for data:
  for (p in 1:N_patients)
  {
    # probability of belonging to the bad group
    bad.p[p] ~ dbern(prob.of.bad.hosp)
    # index the bad group in ac.prob
    index.bad.p[p] <- bad.p[p] + 1

    for(a in 1:antibiotic_classes)
    {
      # Response is different for each antibiotic and depending on
      # which pop it's from
      response[h_GUID[p],a] ~
        dbern(acsm.prob[a,
                        index.bad.p[p],
                        clinical[sample_type[h_sample_GUID[p]]],
                        sample_season[h_sample_GUID[p]]])
    }
  }


  for (gp in 1:N_gp)
  {
    # probability of belonging to the bad group
    bad.gp[gp] ~ dbern(prob.of.bad.gp)
    # index the bad group in ac.prob
    index.bad.gp[gp] <- bad.gp[gp] + 1

    for(a in 1:antibiotic_classes)
    {
      # Response is different for each antibiotic and depending on
      # which pop it's from
      response[gp_GUID[gp],a] ~ dbern(asm.gp.prob[a,
                                                  index.bad.gp[gp],
                                                  sample_season[gp_sample_GUID[gp]]])
    }
  }

  for (v in 1:N_volunteers)
  {
    # probability of belonging to the bad group
    bad.v[v] ~ dbern(prob.of.bad.vol)
    # index the bad group in ac.prob
    index.bad.v[v] <- bad.v[v] + 1

    for(a in 1:antibiotic_classes)
    {
      # Response is different for each antibiotic and depending on
      # which pop it's from
      response[v_GUID[v],a] ~ dbern(asm.v.prob[a,
                                               index.bad.v[v],
                                               sample_season[v_sample_GUID[v]]])
    }
  }

  for (o in 1:N_outpatients)
  {
    # probability of belonging to the bad group
    bad.o[o] ~ dbern(prob.of.bad.out)
    # index the bad group in ac.prob
    index.bad.o[o] <- bad.o[o] + 1

    for(a in 1:antibiotic_classes)
    {
      # Response is different for each antibiotic and depending on
      # which pop it's from
      response[o_GUID[o],a] ~ dbern(asm.o.prob[a,
                                               index.bad.o[o],
                                               sample_season[o_sample_GUID[o]]])
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

    for(b in 1:2)
    { # good bad
      for (c in c(ncarr, nclin)) # 1, 2!
      {
        ac.effect[a,b,c] ~ dnorm(antibiotic.class.effect[a,b], tau.clin)

        for(s in 1:N_sample_season)
        {
          logit(acsm.prob[a,b,c,s]) <- ac.effect[a,b,c] + sampleseason.effect[s]
        }
      }

      a.gp.effect[a,b] <- ac.effect[a,b,gp_clinical]
      a.v.effect[a,b] <- ac.effect[a,b,v_clinical]
      a.o.effect[a,b] <- ac.effect[a,b,o_clinical]

      for(s in 1:N_sample_season)
      {
        logit(asm.gp.prob[a,b,s]) <- a.gp.effect[a,b] + sampleseason.effect[s]
        logit(asm.v.prob[a,b,s]) <- a.v.effect[a,b] + sampleseason.effect[s]
        logit(asm.o.prob[a,b,s]) <- a.o.effect[a,b] + sampleseason.effect[s]

      }
    }
  }

  for(s in 1:N_sample_season)
  {
    sampleseason.effect[s] ~ dnorm(0, tau.sampleseason)
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

  # Prior values for precision
  tau.class ~ dgamma(0.001, 0.001)
  tau.clin ~ dgamma(0.001, 0.001)
  tau.sampleseason ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.class <- sqrt(1/tau.class)
  sd.clin <- sqrt(1/tau.clin)
  sd.sampleseason <- sqrt(1/tau.sampleseason)

  #monitor# full.pd, dic, deviance, a.prob, ac.prob, prob.of.bad.hosp, prob.of.bad.gp, prob.of.bad.vol, prob.of.bad.out, bad.p, bad.gp, bad.v, bad.o, intercept, sd.class, sd.clin, sd.sampleseason
}
