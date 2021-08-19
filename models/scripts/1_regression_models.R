#' ---
#' title: Regression models
#' author: Sonia Mitchell
#' date: 14 Jan 2020
#' output:
#'   html_document:
#'     theme: paper
#'     format: html_clean
#'     code_folding: show
#'     highlight: pygments
#' ---
#'
#' Last updated: `r Sys.time()`

#+ r setup, include=FALSE

library(SpARK)
library(SpARKjags)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)
library(lme4)
library(lubridate)

bin_ages <- function(x, N) {
  df <- x %>% dplyr::mutate(AGE = as.numeric(interval(ymd(DATE_OF_BIRTH),
                                                      ymd(Sys.Date())),
                                             unit = "years"))
  max_age <- ceiling(max(df$AGE))
  jump <- floor(max_age / N)
  tmp <- seq(1, jump*N, jump)
  df %>%
    dplyr::mutate(AGE_GROUP = cut(AGE, c(tmp, max_age))) %>%
    dplyr::mutate(AGE_GROUP = as.factor(AGE_GROUP))
}

get_data <- function(categ, antclass) {
  id <- atb %>%
    filter(used_MIC == "yes",
           Classification == antclass,
           Interpretation == "R") %$%
    UNIQUE_SPARK_ID %>%
    unique()
  df <- meta %>%
    filter(Category == categ) %>%
    mutate(response = case_when(GUID %in% id ~ 1,
                                T ~ 0))
  if(categ == "animal") {
    df %<>%
      select(ASSOCIATED_GROUP, Livestock, Companion_animal,
             Wild_animal, PATHOGEN, response)

  }else if(categ == "human") {
    df %<>%
      dplyr::filter(DATE_OF_BIRTH != "XXXX",
                    DATE_OF_BIRTH != "1986-02-97") %>%
      bin_ages(10) %>%
      select(Clinical, SPECIFIC_GROUP, HOSPITAL_WARD, PATHOGEN,
             SEX, AGE, AGE_GROUP, response)

  }else if(categ == "other") {
    df %<>%
      select(Plant_pre_harvest, Plant_post_harvest, Water,
             Soil, Other, Hospital, ORIGIN, response)
  }else if(categ == "all") {
    df %<>%
      select(Category, PATHOGEN, response)
  }

  df
}

plotROC <- function(dataset, glm, title) {
  prob <- predict(glm, type = c("response"))
  dataset$prob <- prob
  tmp <- pROC::roc(response ~ prob, dataset, quiet = T)
  auc <- paste0("AUC: ", round(as.numeric(tmp$auc), 4))
  aic <- paste0("AIC: ", round(AIC(glm), 1))
  tmp <- cbind.data.frame(x = rev(1 - tmp$specificities),
                          y = rev(tmp$sensitivities),
                          facet = title)

  ggplot() + geom_line(aes(x, y), tmp) +
    theme_bw() + labs(x = "1-specificity", y = "Sensitivity") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    facet_wrap(~facet) +
    annotate("text", x = .8, y = .2, label = auc) +
    annotate("text", x = .8, y = .1, label = aic)
}

#+

# Metadata ----------------------------------------------------------------

meta <- METAdata %>%
  mutate(sample_GUID = gsub("_C[0-9]$", "", GUID))

# Number of GUIDs with metadata entries
length(unique(meta$GUID))
# Number of sample_GUIDs with metadata entries
length(unique(meta$sample_GUID))

# Number of sample_GUIDs with metadata entries in each category and
# Proportion of sample_GUIDs with metadata entries in each category
meta %>%
  select(sample_GUID, Category) %>%
  unique() %>%
  group_by(Category) %>%
  summarise(count = n(),
            proportion = n() / nrow(.))





# Antibiograms ------------------------------------------------------------

atb <- ATBdata %>%
  mutate(sample_GUID = gsub("_C[0-9]$", "", UNIQUE_SPARK_ID))

# Number of GUIDs with antibiogram profiles
length(unique(atb$UNIQUE_SPARK_ID))
# Number of sample_GUIDs with antibiogram profiles
length(unique(atb$sample_GUID))

# Number of sample_GUIDs with antibiograms in each category and
# Proportion of sample_GUIDs with antibiogram profiles in each category
meta %>%
  filter(sample_GUID %in% atb$sample_GUID) %>%
  select(sample_GUID, Category) %>%
  unique() %>%
  group_by(Category) %>%
  summarise(count = n(),
            proportion = n() / nrow(.))


# Number of sample_GUIDs with antibiogram profiles in each category and
# Proportion of sample_GUIDs with antibiogram profiles in each category
# relative to the total number of sample_GUIDs with metadata entries
meta %>%
  filter(sample_GUID %in% atb$sample_GUID) %>%
  select(sample_GUID, Category) %>%
  unique() %>%
  group_by(Category) %>%
  summarise(n_atb = n()) %>%
  merge(meta %>%
          select(sample_GUID, Category) %>%
          unique() %>%
          group_by(Category) %>%
          summarise(n_meta = n())) %>%
  mutate(percentage = n_atb / n_meta)


# Distribution of numbers of antibiotics tested (that use MIC values)
atb %>%
  filter(used_MIC == "yes",
         Interpretation != "ND") %>%
  select(sample_GUID, Antibiotic_name) %>%
  unique() %>%
  group_by(sample_GUID) %>%
  summarise(count = n()) %$%
  count %>%
  summary()


# Find most common resistance profiles (by GUIDs)
tmp <- atb %>%
  filter(used_MIC == "yes",
         Interpretation != "ND")
antibiotics <- as.character(unique(tmp$Antibiotic_name))

profiles <- tmp %>%
  reshape2::dcast(UNIQUE_SPARK_ID + sample_GUID ~ Antibiotic_name,
                  value.var = "Interpretation") %>%
  mutate_at(vars(one_of(antibiotics)), ~ case_when(.=="R" ~ 1, T ~ 0)) %>%
  mutate(profile = apply(.[antibiotics], 1, function(x)
    paste0(x, collapse = "")))

top_profiles <- profiles %>%
  group_by(profile) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
top_profiles

# Most common antibiotic profiles (note that the 3 most common profile is a
# resistance to nothing)
lapply(1:10, function(x) {
  profiles %>%
    filter(profile == top_profiles$profile[x]) %>%
    .[antibiotics] %>%
    colSums() %>%
    .[. != 0] %>%
    names()
})




# Sequences ---------------------------------------------------------------

kleb <- KLEBdata %>%
  mutate(sample_GUID = gsub("_C[0-9]$", "", GUID))

# Number of GUIDs with seqences
length(unique(kleb$GUID))
# Number of sample_GUIDs with seqences
length(unique(kleb$sample_GUID))

# Number of sample_GUIDs with sequences in each category and
# Proportion of sample_GUIDs with sequences in each category
meta %>%
  select(sample_GUID, Category) %>%
  unique() %>%
  filter(sample_GUID %in% kleb$sample_GUID) %>%
  group_by(Category) %>%
  summarise(count = n(),
            proportion = n() / nrow(.))

# Number of sample_GUIDs with sequences from each pathogen and
# Proportion of sample_GUIDs with sequences from each pathogen
kleb %>%
  select(sample_GUID, species) %>%
  unique() %>%
  group_by(species) %>%
  summarise(count = n(),
            proportion = n() / nrow(.)) %>%
  arrange(desc(proportion))



# Start of modelling section ----------------------------------------------

antibiotic_class <- atb %>%
  filter(used_MIC == "yes") %$%
  Classification %>%
  unique()
antibiotic_class

resistances <- lapply(antibiotic_class, function(x) {
  id <- atb %>%
    filter(Classification == x,
           Interpretation == "R") %$%
    UNIQUE_SPARK_ID %>%
    unique()

  data <- meta %>%
    mutate(response = case_when(GUID %in% id ~ "resistant",
                                T ~ "not resistant"))

  data %>%
    group_by(Category, response) %>%
    summarise(Count = n()) %>%
    complete(Category = unique(Category),
             response = c("resistant", "not resistant"),
             fill = list(Count = 0)) %>%
    ungroup() %>%
    mutate(class = x)
}) %>%
  do.call(rbind.data.frame, .)

# Put plots in this order
sort_order <- resistances %>%
  filter(response == "resistant") %>%
  group_by(class) %>%
  summarise(total = sum(Count)) %>%
  arrange(desc(total)) %$%
  class
resistances %<>%
  mutate(class = factor(class, levels = sort_order))

# Plot the number of GUIDs resistant to (red) and not resistant to (grey) each
# antibiotic
resistances %>%
  ggplot() + theme_minimal() + facet_wrap(~class) +
  geom_bar(aes(x = Category, y = Count, fill = response),
           colour = "black", stat = "identity") +
  scale_fill_manual(values = c("grey", "red"))

# Plot the resistance to each antibiotic as a percentage (of GUIDs)
resistances %>%
  group_by(Category, class) %>%
  summarise(percentage = Count[response == "resistant"] / sum(Count)) %>%
  ggplot() + theme_minimal() + facet_wrap(~class) +
  geom_bar(aes(x = Category, y = percentage),
           colour = "black", fill = "red", stat = "identity") +
  labs(title = "Percentage of resistant samples")


# List values used in the above plot (percentage resistant) categorised
# as human, animal, or other
resistances %>%
  group_by(Category, class) %>%
  summarise(percentage = Count[response == "resistant"] / sum(Count)) %>%
  data.frame()


# List percentage resistant (combining all human, animal, and other samples)
resistances %>%
  group_by(response, class) %>%
  summarise(total = sum(Count)) %>%
  group_by(class) %>%
  summarise(percentage = total[response == "resistant"] / sum(total)) %>%
  data.frame()



# Penicillin --------------------------------------------------------------
# --- Animal
df <- get_data("animal", "Penicillin")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Penicillin", deparse(best_model$formula))
apen <- plotROC(df, best_model, title)
apen

# --- Human
df <- get_data("human", "Penicillin")

full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
best_model <- MASS::stepAIC(full_model)
best_model$formula
# Adding ward to this as a random effect results in an error

title <- gsub("response", "Penicillin", deparse(best_model$formula))
hpen <- plotROC(df, best_model, title)
hpen

# --- Other
df <- get_data("other", "Penicillin")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
best_model <- MASS::stepAIC(full_model)

title <- paste(deparse(best_model$formula), collapse = "") %>%
  gsub("response", "Penicillin", .) %>%
  gsub("    ", "", .)
open <- plotROC(df, best_model, title)
open


# Penicillin Combination --------------------------------------------------
# --- Animal
df <- get_data("animal", "Penicillin Combination")
full_model <- glm(response ~ ., binomial(link = "logit"),  df)

tmp <- alias(full_model)$Complete
colinear_variables <- tmp[, which(colSums(tmp) != 0), drop = F]
df %<>% select(-Wild_animal)
full_model <- glm(response ~ ., binomial(link = "logit"), df)
df %<>% select(-PATHOGEN)
full_model <- glm(response ~ ., binomial(link = "logit"), df)
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "PenC", deparse(best_model$formula))
apenc <- plotROC(df, best_model, title)
apenc

# --- Human
df <- get_data("human", "Penicillin Combination")
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "PenC", deparse(best_model$formula))
hpenc <- plotROC(df, best_model, title)
hpenc

# --- Other
df <- get_data("other", "Penicillin Combination")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "PenC", deparse(best_model$formula))
openc <- plotROC(df, best_model, title)
openc



# Cephalosporin -----------------------------------------------------------
# --- Animal
df <- get_data("animal", "Cephalosporin")
full_model <- glm(response ~ ., binomial(link = "logit"), df)

tmp <- alias(full_model)$Complete
colinear_variables <- tmp[, which(colSums(tmp) != 0), drop = F]
df %<>% select(-PATHOGEN, -Wild_animal)
full_model <- glm(response ~ ., binomial(link = "logit"), df)
df %<>% filter(ASSOCIATED_GROUP != "cheese")
full_model <- glm(response ~ ., binomial(link = "logit"), df)

model1 <- glm(response ~ ASSOCIATED_GROUP, binomial(link = "logit"), df)
model2 <- glm(response ~ Livestock + Companion_animal,
              binomial(link = "logit"), df)
AIC(model1, model2) # model1 is better

title <- gsub("response", "Ceph", deparse(best_model$formula))
acep <- plotROC(df, model1, title)
acep

# --- Human
df <- get_data("human", "Cephalosporin")
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Ceph", deparse(best_model$formula))
hcep <- plotROC(df, best_model, title)
hcep

# --- Other
df <- get_data("other", "Cephalosporin")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
best_model <- MASS::stepAIC(full_model)

title <- paste(deparse(best_model$formula), collapse = "") %>%
  gsub("response", "Ceph", .) %>%
  gsub("    ", "", .)
ocep <- plotROC(df, best_model, title)
ocep




# Fosfomycin --------------------------------------------------------------
# --- Animal
df <- get_data("animal", "Fosfomycin")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Fos", deparse(best_model$formula))
afos <- plotROC(df, best_model, title)
afos

# --- Human
df <- get_data("human", "Fosfomycin")
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Fos", deparse(best_model$formula))
hfos <- plotROC(df, best_model, title)
hfos

# --- Other
df <- get_data("other", "Fosfomycin")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
best_model <- MASS::stepAIC(full_model)

title <- paste(deparse(best_model$formula), collapse = "") %>%
  gsub("response", "Fos", .) %>%
  gsub("    ", "", .)
ofos <- plotROC(df, best_model, title)
ofos




# Penicillin (Penams) -----------------------------------------------------
# --- Animal
df <- get_data("animal", "Penicillin (Penams)")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Pen (Pen)", deparse(best_model$formula))
app <- plotROC(df, best_model, title)
app

# --- Human
df <- get_data("human", "Penicillin (Penams)")
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Pen (Pen)", deparse(best_model$formula))
hpp <- plotROC(df, best_model, title)
hpp

# --- Other
df <- get_data("other", "Penicillin (Penams)")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
best_model <- MASS::stepAIC(full_model)

title <- paste(deparse(best_model$formula), collapse = "") %>%
  gsub("response", "Pen (Pen)", .) %>%
  gsub("    ", "", .)
opp <- plotROC(df, best_model, title)
opp





# Fluoroquinolone ---------------------------------------------------------
# --- Animal
df <- get_data("animal", "Fluoroquinolone")
full_model <- glm(response ~ ., binomial(link = "logit"), df)

tmp <- alias(full_model)$Complete
colinear_variables <- tmp[, which(colSums(tmp) != 0), drop = F]
df %<>% select(-PATHOGEN, -Wild_animal)
full_model <- glm(response ~ ., binomial(link = "logit"), df)
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Flu", deparse(best_model$formula))
aflu <- plotROC(df, best_model, title)
aflu

# --- Human
df <- get_data("human", "Fluoroquinolone")
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Flu", deparse(best_model$formula))
hflu <- plotROC(df, best_model, title)
hflu

# --- Other
df <- get_data("other", "Fluoroquinolone")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
table(df$response)




# Trimethoprim/Sulfamethoxazole -------------------------------------------
# --- Animal
df <- get_data("animal", "Trimethoprim/Sulfamethoxazole")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
full_model <- glm(response ~ ., binomial(link = "logit"), df, maxit = 100)

tmp <- alias(full_model)$Complete
colinear_variables <- tmp[, which(colSums(tmp) != 0), drop = F]
df %<>% select(-PATHOGEN, -Wild_animal)
full_model <- glm(response ~ ., binomial(link = "logit"), df)
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Tri/Sul", deparse(best_model$formula))
ats <- plotROC(df, best_model, title)
ats

# --- Human
df <- get_data("human", "Trimethoprim/Sulfamethoxazole")
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))

tmp <- alias(full_model)$Complete
df %<>% select(-PATHOGEN)
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Tri/Sul", deparse(best_model$formula))
hts <- plotROC(df, best_model, title)
hts

# --- Other
df <- get_data("other", "Fluoroquinolone")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
table(df$response)







# Aminoglycoside ----------------------------------------------------------
# --- Animal
df <- get_data("animal", "Aminoglycoside")

full_model <- glm(response ~ ., binomial(link = "logit"), df)
full_model <- glm(response ~ ., binomial(link = "logit"), df, maxit = 100)

tmp <- alias(full_model)$Complete
colinear_variables <- tmp[, which(colSums(tmp) != 0), drop = F]
df %<>% select(-PATHOGEN, -Wild_animal)
full_model <- glm(response ~ ., binomial(link = "logit"), df)

model1 <- glm(response ~ ASSOCIATED_GROUP, binomial(link = "logit"), df)
model2 <- glm(response ~ Livestock + Companion_animal,
              binomial(link = "logit"), df)
AIC(model1, model2) # model2 is better

title <- gsub("response", "Amin", deparse(best_model$formula))
aam <- plotROC(df, model2, title)
aam

# --- Human
df <- get_data("human", "Aminoglycoside")
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))

tmp <- alias(full_model)$Complete
df %<>% select(-PATHOGEN)
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Amin", deparse(best_model$formula))
ham <- plotROC(df, best_model, title)
ham

# --- Other
df <- get_data("other", "Aminoglycoside")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
table(df$response)



# Nitrofurantoin ----------------------------------------------------------
# --- Animal
df <- get_data("animal", "Nitrofurantoin")
full_model <- glm(response ~ ., binomial(link = "logit"), df)

tmp <- alias(full_model)$Complete
colinear_variables <- tmp[, which(colSums(tmp) != 0), drop = F]
df %<>% select(-PATHOGEN, -Wild_animal)
full_model <- glm(response ~ ., binomial(link = "logit"), df)
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Nitr", deparse(best_model$formula))
anit <- plotROC(df, best_model, title)
anit

# --- Human
df <- get_data("human", "Nitrofurantoin")
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))

tmp <- alias(full_model)$Complete
df %<>% select(-PATHOGEN)
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Nitr", deparse(best_model$formula))
hnit <- plotROC(df, best_model, title)
hnit

# --- Other
df <- get_data("other", "Nitrofurantoin")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Nitr", deparse(best_model$formula))
onit <- plotROC(df, best_model, title)
onit



# Carbapenem --------------------------------------------------------------
# --- Animal
df <- get_data("animal", "Carbapenem")
table(df$response)

# --- Human
df <- get_data("human", "Carbapenem")
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))

tmp <- alias(full_model)$Complete
df %<>% select(-PATHOGEN)
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
df %<>% select(-AGE_GROUP)
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
best_model <- MASS::stepAIC(full_model)
deparse(best_model$formula)

model2 <-  glmer(response ~ Clinical + SPECIFIC_GROUP + SEX + AGE + (1|HOSPITAL_WARD),
                 df, binomial(link = "logit"))
# summary(allFit(model2))

title <- gsub("response", "Carb", deparse(best_model$formula))
hcar <- plotROC(df, best_model, title)
hcar

# --- Other
df <- get_data("other", "Carbapenem")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
table(df$response)





# Trimethoprim ------------------------------------------------------------
# --- Animal
df <- get_data("animal", "Trimethoprim")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
full_model <- glm(response ~ ., binomial(link = "logit"), df, maxit = 100)

tmp <- alias(full_model)$Complete
colinear_variables <- tmp[, which(colSums(tmp) != 0), drop = F]
df %<>% select(-PATHOGEN, -Wild_animal)
full_model <- glm(response ~ ., binomial(link = "logit"), df)
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Tri", deparse(best_model$formula))
atri <- plotROC(df, best_model, title)
atri

# --- Human
df <- get_data("human", "Trimethoprim")
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Tri", deparse(best_model$formula))
htri <- plotROC(df, best_model, title)
htri

# --- Other
df <- get_data("other", "Trimethoprim")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
table(df$response)




# Tetracycline ------------------------------------------------------------
# --- Animal
df <- get_data("animal", "Tetracycline")
full_model <- glm(response ~ ., binomial(link = "logit"), df)

tmp <- alias(full_model)$Complete
colinear_variables <- tmp[, which(colSums(tmp) != 0), drop = F]
df %<>% select(-PATHOGEN, -Wild_animal)
full_model <- glm(response ~ ., binomial(link = "logit"), df)

model1 <- glm(response ~ ASSOCIATED_GROUP, binomial(link = "logit"), df)
model2 <- glm(response ~ Livestock + Companion_animal,
              binomial(link = "logit"), df)
AIC(model1, model2) # model2 is better

title <- gsub("response", "Tetr", deparse(best_model$formula))
atet <- plotROC(df, model2, title)
atet

# --- Human
df <- get_data("human", "Tetracycline")
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Tetr", deparse(best_model$formula))
htet <- plotROC(df, best_model, title)
htet

# --- Other
df <- get_data("other", "Tetracycline")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
table(df$response)




# Monobactam --------------------------------------------------------------
# --- Animal
df <- get_data("animal", "Monobactam")
full_model <- glm(response ~ ., binomial(link = "logit"), df)

tmp <- alias(full_model)$Complete
colinear_variables <- tmp[, which(colSums(tmp) != 0), drop = F]
df %<>% select(-PATHOGEN, -Wild_animal)
full_model <- glm(response ~ ., binomial(link = "logit"), df)
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Mono", deparse(best_model$formula))
amon <- plotROC(df, best_model, title)
amon

# --- Human
df <- get_data("human", "Monobactam")
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))

tmp <- alias(full_model)$Complete
df %<>% select(-PATHOGEN)
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Mono", deparse(best_model$formula))
hmon <- plotROC(df, best_model, title)
hmon

# --- Other
df <- get_data("other", "Monobactam")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
table(df$response)



# Colistin ----------------------------------------------------------------
# --- Animal
df <- get_data("animal", "Colistin")
full_model <- glm(response ~ ., binomial(link = "logit"), df)

tmp <- alias(full_model)$Complete
colinear_variables <- tmp[, which(colSums(tmp) != 0), drop = F]
df %<>% select(-PATHOGEN, -Wild_animal)
full_model <- glm(response ~ ., binomial(link = "logit"), df)

model1 <- glm(response ~ ASSOCIATED_GROUP, binomial(link = "logit"), df)
model2 <- glm(response ~ Livestock + Companion_animal,
              binomial(link = "logit"), df)
AIC(model1, model2) # model2 is better

title <- gsub("response", "Col", deparse(best_model$formula))
acol <- plotROC(df, model2, title)
acol

# --- Human
df <- get_data("human", "Colistin")
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))

tmp <- alias(full_model)$Complete
df %<>% select(-PATHOGEN)
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
df %<>% select(-AGE_GROUP)
full_model <- glm(response ~ ., binomial(link = "logit"), df %>%
                    select(-HOSPITAL_WARD))
best_model <- MASS::stepAIC(full_model)

title <- gsub("response", "Col", deparse(best_model$formula))
hcol <- plotROC(df, best_model, title)
hcol

# --- Other
df <- get_data("other", "Colistin")
full_model <- glm(response ~ ., binomial(link = "logit"), df)
table(df$response)





# animals <- cowplot::plot_grid(aam, acep, acol, aflu, afos, amon, anit, apen,
#                               apenc, app, atet, atri, ats, align = "vh",
#                               axis = "tb", ncol = 2, labels = letters[1:13])
# ggsave("animals.pdf", animals, width = 12, height = 24)
#
# humans <- cowplot::plot_grid(ham, hcar, hcep, hcol, hflu, hfos, hmon, hnit,
#                              hpen, hpenc, hpp, htet, htri, hts, align = "vh",
#                              axis = "tb", ncol = 2, labels = letters[1:14])
# ggsave("humans.pdf", humans, width = 12, height = 24)
#
# others <- cowplot::plot_grid(ocep, ofos, onit, open, openc, opp, align = "vh",
#                              axis = "tb", ncol = 2, labels = letters[1:6])
# ggsave("others.pdf", others, width = 12, height = 10)



