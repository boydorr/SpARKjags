library(SpARKjags)
library(runjags)
library(tidybayes)
library(dplyr)
library(ggplot2)
library(forcats)

data <- jags_data(classification = "all",
                  categories = "human",
                  pathogen = "Klebsiella pneumoniae",
                  removeQuinPen = T)

path <- "/Users/Soniam/Desktop/git/SpARK/SpARKjags/results/full_models/a.rds"
m <- get_model(path)

m %>% spread_draws(a.prob[class])

# Plotting posterior medians and 66% and 95% intervals
m %>%
  spread_draws(a.prob[class]) %>%
  mutate(class = factor(class)) %>%
  ggplot(aes(y = fct_rev(class), x = a.prob)) +
  stat_pointinterval()

# Combinations of posterior intervals and densities, drawn as violin plots
m %>%
  spread_draws(a.prob[class]) %>%
  mutate(class = factor(class)) %>%
  ggplot(aes(y = fct_rev(class), x = a.prob)) +
  stat_eye()

# Intervals (90% and 50%) with posterior densities
m %>%
spread_draws(a.prob[class]) %>%
  mutate(class = factor(class)) %>%
  ggplot(aes(y = fct_rev(class), x = a.prob)) +
  stat_halfeye(.width = c(.90, .5))



m %>%
  spread_draws(a.prob[class]) %>%
  mutate(class = factor(class)) %>%
  ggplot(aes(y = fct_rev(class), x = a.prob, fill = stat(abs(x) < .3))) +
  stat_halfeye() +
  geom_vline(xintercept = c(.2, .3), linetype = "dashed") +
  scale_fill_manual(values = c("gray80", "skyblue"))







