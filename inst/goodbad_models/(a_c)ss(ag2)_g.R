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
                        sample_season[h_sample_GUID[p]],
                        age_group2[h_sample_GUID[p]],
                        gender[h_sample_GUID[p]]])
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
      response[gp_GUID[gp],a] ~
        dbern(asm.gp.prob[a,
                          index.bad.gp[gp],
                          sample_season[gp_sample_GUID[gp]],
                          age_group2[gp_sample_GUID[gp]],
                          gender[gp_sample_GUID[gp]]])
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
                                               sample_season[v_sample_GUID[v]],
                                               age_group2[v_sample_GUID[v]],
                                               gender[v_sample_GUID[v]]])
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
                                               sample_season[o_sample_GUID[o]],
                                               age_group2[o_sample_GUID[o]],
                                               gender[o_sample_GUID[o]]])
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

  for(s in 1:N_sample_season)
  {
    sampleseason.effect[s] ~ dnorm(0, tau.sampleseason)
  }

  for (g in 1:N_age_group)
  {
    agegroup.effect[g] ~ dnorm(0, tau.agegroup)

    for(d in 1:N_gender)
    {
      age.gender.effect[g,d] ~ dnorm(agegroup.effect[g], tau.gender)
      logit(age.gender.prob[g,d]) <- age.gender.effect[g,d]
    }
  }

  for(a in 1:antibiotic_classes)
  {
    for(b in 1:2) # good bad
    {
      a.gp.effect[a,b] <- ac.effect[a,b,gp_clinical]
      a.v.effect[a,b] <- ac.effect[a,b,v_clinical]
      a.o.effect[a,b] <- ac.effect[a,b,o_clinical]
      logit(a.gp.prob[a,b]) <- a.gp.effect[a,b]
      logit(a.v.prob[a,b]) <- a.v.effect[a,b]
      logit(a.o.prob[a,b]) <- a.o.effect[a,b]

      for (c in c(ncarr, nclin)) # 1, 2!
      {
        ac.effect[a,b,c] ~ dnorm(antibiotic.class.effect[a,b], tau.clin)
        logit(ac.prob[a,b,c]) <- ac.effect[a,b,c]
      }
    }
  }

  for(s in 1:N_sample_season)
  {
    for (g in 1:N_age_group)
    {
      for (d in genders)
      {
        ss.age.gender.effect[s,g,d] <-
          sampleseason.effect[s] + age.gender.effect[g,d]

        for(a in 1:antibiotic_classes)
        {
          for(b in 1:2) # good bad
          {

            for (c in c(ncarr, nclin)) # 1, 2!
            {
              logit(acsm.prob[a,b,c,s,g,d]) <-
                ac.effect[a,b,c] + ss.age.gender.effect[s,g,d]
            }

            logit(asm.gp.prob[a,b,s,g,d]) <-
              a.gp.effect[a,b] + ss.age.gender.effect[s,g,d]
            logit(asm.v.prob[a,b,s,g,d]) <-
              a.v.effect[a,b] + ss.age.gender.effect[s,g,d]
            logit(asm.o.prob[a,b,s,g,d]) <-
              a.o.effect[a,b] + ss.age.gender.effect[s,g,d]
          }
        }
      }
    }
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
  tau.agegroup ~ dgamma(0.001, 0.001)
  tau.gender ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.class <- sqrt(1/tau.class)
  sd.clin <- sqrt(1/tau.clin)
  sd.sampleseason <- sqrt(1/tau.sampleseason)
  sd.agegroup <- sqrt(1/tau.agegroup)
  sd.gender <- sqrt(1/tau.gender)

  #monitor# full.pd, dic, deviance, a.prob, ac.prob, a.gp.prob, a.v.prob, a.o.prob, prob.of.bad.hosp, prob.of.bad.gp, prob.of.bad.vol, prob.of.bad.out, bad.p, bad.gp, bad.v, bad.o, intercept, sd.class, sd.clin, sd.sampleseason, sd.agegroup, sd.gender
}
