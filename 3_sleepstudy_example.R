# This example is adapted from:
# Pavel Chernyavskiy’s lecture series on Applied Bayesian Statistics (Fall 2020).
# https://www.youtube.com/@applied_bayesian_stats6070

################## Sleepstudy Example #####################


library(lme4)
library("ggplot2")
library("dplyr")


############## DATA  ############## 

# Daten laden
data("sleepstudy")
str(sleepstudy)

############## ICC  ############## 

# Fit random intercept model
icc_model <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)

# Calculate ICC
variance_components <- as.data.frame(VarCorr(icc_model))
icc_value <- variance_components$vcov[1] / sum(variance_components$vcov)

# Calculate global mean
global_mean <- mean(sleepstudy$Reaction)

# Calculate subject means
subject_means <- sleepstudy %>%
  group_by(Subject) %>%
  summarise(mean_reaction = mean(Reaction))

# Create ICC plot
ggplot(sleepstudy, aes(x = as.factor(Subject), y = Reaction)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = global_mean, color = "blue", linetype = "dashed", size = 0.6) +
  geom_point(data = subject_means, aes(x = as.factor(Subject), y = mean_reaction),
             color = "red", size = 2) +
  theme_minimal() +
  labs(title = "",
       subtitle = paste("ICC =", round(icc_value, 1), "| Overall Mean =", round(global_mean, 1), "ms"),
       x = "Subject", y = "Reaction Time (ms)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############## SUBJECT INDIVIDUELL ############## 

# Calculate daily means
daily_means <- sleepstudy %>%
  group_by(Days) %>%
  summarise(mean_reaction = mean(Reaction))

# Create time trend plot
ggplot() +
  geom_line(data = sleepstudy, aes(x = Days, y = Reaction, group = Subject, color = factor(Subject)),
            alpha = 0.5, size = 1) +
  geom_line(data = daily_means, aes(x = Days, y = mean_reaction),
            color = "red", linetype = "longdash", size = 1.5) +
  scale_x_continuous(breaks = seq(0, max(sleepstudy$Days), by = 2)) +
  theme_minimal() +
  labs(title = "Reaction Time Over Days",
       subtitle = "Individual subject trends | Mean reaction time in red (dashed)",
       x = "Days", y = "Reaction Time (ms)", color = "Subject") +
  theme(legend.position = "none")



############## MODEL 1 ############## 


# Model 1: random subject intercept + days 
set.seed(05)

my_priors <- c(
  set_prior ("normal(0, 50)", class = "sd"),
  set_prior ("normal (250, 50)", class = "Intercept"),
  set_prior ("normal (10, 10)", class ="b") ,
  set_prior ("normal (0, 50)", class = "sigma"))

set.seed(05)
library("brms")
fit_sleep_m1 <- brm(formula = Reaction ~ (1|Subject) + Days,
                 data = sleepstudy,
                 family = gaussian(),
                 prior = my_priors,
                 seed = 10,
                 control = list(adapt_delta = 0.8, max_treedepth=10),
                 warmup=1000, iter=3000, thin=1, chains=3, cores=3)

summary(fit_sleep_m1)

coef(fit_sleep_m1)

fixef(fit_sleep_m1)
ranef(fit_sleep_m1)


############## POSTERIOR PREDICITVE CHECK MODEL 1 ############## 
pp_check(fit_sleep_m1, ndraws=100)
pp_check(fit_sleep_m1, ndraws = 100, type = "stat_grouped", group = "Subject", stat="mean") + theme_bw()



############## FITTED LINES IN ONE PLOT MODEL 1 ############## 
# Add predicted values to data
sleepstudy$predicted <- fitted(fit_sleep_m1)[, "Estimate"]

# Compute population-level fit (mean prediction per day)
pop_fit <- sleepstudy %>%
  group_by(Days) %>%
  summarise(pop_mean = mean(predicted))

# Plot all subject-level lines and population average
ggplot() +
  geom_line(data = sleepstudy, aes(x = Days, y = predicted, group = Subject, color = Subject),
            alpha = 0.6, size = 1) +
  geom_line(data = pop_fit, aes(x = Days, y = pop_mean),
            color = "red", linetype = "longdash", size = 1.5) +
  scale_x_continuous(breaks = seq(0, max(sleepstudy$Days), by = 2)) +
  theme_minimal() +
  labs(x = "Days", y = "Reaction Time (ms)", color = "Subject") +
  theme(legend.position = "none")

############## SUBJECT SPECIFIC FITTED LINES MODEL 1 ############## 

# Add predicted values to data
sleepstudy$predicted <- fitted(fit_sleep_m1)[, "Estimate"]

# Define custom colors per subject
subject_colors <- c(
  "308" = "#1f77b4", "309" = "#ff7f0e", "310" = "#2ca02c", "330" = "#d62728",
  "331" = "#9467bd", "332" = "#8c564b", "333" = "#e377c2", "334" = "#7f7f7f",
  "335" = "#bcbd22", "337" = "#17becf", "349" = "#aec7e8", "350" = "#ffbb78",
  "351" = "#98df8a", "352" = "#ff9896", "369" = "#c5b0d5", "370" = "#c49c94",
  "371" = "#f7b6d2", "372" = "#c7c7c7"
)

# Faceted plot: actual and predicted values per subject
ggplot(sleepstudy, aes(x = Days, y = Reaction)) +
  geom_point(alpha = 0.8) +
  geom_line(aes(y = predicted, color = factor(Subject)), size = 1) +
  scale_color_manual(values = subject_colors) +
  scale_x_continuous(breaks = seq(0, max(sleepstudy$Days), by = 2)) +
  theme_minimal() +
  labs(x = "Days", y = "Estimated Reaction Time") +
  facet_wrap(~ Subject) +
  theme(legend.position = "none")





############## MODEL 2 ############## 
# Model 2: random subject intercept + random subject slope w/ respect to day 
my_priors <- c(
  set_prior ("normal(0, 50)", class = "sd", coef = "Intercept", group = "Subject"),
  set_prior ("normal(0, 20)", class = "sd", coef = "Days", group = "Subject"),
  set_prior ("normal (250, 50)", class = "Intercept"),
  set_prior ("normal (10, 10)", class ="b") ,
  set_prior ("normal (0, 50)", class = "sigma"))


library("brms")
fit_sleep_m2 <- brm(formula = Reaction ~ (Days||Subject) + Days,
                 data = sleepstudy,
                 family = gaussian(),
                 prior = my_priors,
                 seed = 10,
                 control = list(adapt_delta = 0.8, max_treedepth=10),
                 warmup=1000, iter=3000, thin=1, chains=3, cores=3)
summary(fit_sleep_m2)

coef(fit_sleep_m2)



############## FITTED LINES IN ONE PLOT MODEL 2 ############## 

# Extract predicted values from model
sleepstudy$predicted_m2 <- fitted(fit_sleep_m2)[, "Estimate"]

# Compute population-level mean prediction per day
pop_fit_m2 <- sleepstudy %>%
  group_by(Days) %>%
  summarise(pop_mean = mean(predicted_m2))

# Plot 1: Overall view with fitted lines per subject and population average
ggplot() +
  geom_line(data = sleepstudy, aes(x = Days, y = predicted_m2, group = Subject, color = Subject),
            alpha = 0.6, size = 1) +
  geom_line(data = pop_fit_m2, aes(x = Days, y = pop_mean),
            color = "red", linetype = "longdash", size = 1.5) +
  scale_x_continuous(breaks = seq(0, max(sleepstudy$Days), by = 2)) +
  theme_minimal() +
  labs(x = "Days", y = "Reaction Time (ms)", color = "Subject") +
  theme(legend.position = "none")

############## SUBJECT SPECIFIC FITTED LINES MODEL 2 ############## 

# Define subject-specific colors
subject_colors <- c(
  "308" = "#1f77b4", "309" = "#ff7f0e", "310" = "#2ca02c", "330" = "#d62728",
  "331" = "#9467bd", "332" = "#8c564b", "333" = "#e377c2", "334" = "#7f7f7f",
  "335" = "#bcbd22", "337" = "#17becf", "349" = "#aec7e8", "350" = "#ffbb78",
  "351" = "#98df8a", "352" = "#ff9896", "369" = "#c5b0d5", "370" = "#c49c94",
  "371" = "#f7b6d2", "372" = "#c7c7c7"
)

# Plot 2: Faceted view showing individual fitted lines and observations
ggplot(sleepstudy, aes(x = Days, y = Reaction)) +
  geom_point(alpha = 0.8) +
  geom_line(aes(y = predicted_m2, color = factor(Subject)), size = 1) +
  scale_color_manual(values = subject_colors) +
  scale_x_continuous(breaks = seq(0, max(sleepstudy$Days), by = 2)) +
  theme_minimal() +
  labs(x = "Days", y = "Estimated Reaction Time") +
  facet_wrap(~ Subject) +
  theme(legend.position = "none")

############## POSTERIOR PREDICTIVE CHECK MODEL 2 ############## 

pp_check(fit_sleep_m2, ndraws = 100) + xlab('Rection Time') + theme_bw()
pp_check(fit_sleep_m2, ndraws = 100, type = "ribbon_grouped", x="Days", group = "Subject") + theme_bw()



############## WAIC MODEL 1 & 2 ############## 

library(brms)
waic(fit_sleep_m1)
waic(fit_sleep_m2)


############## MODEL 3 ############## 

my_priors <- c(
  set_prior ("lkj(3)", class = "cor"),
  set_prior ("normal(0, 50)", class = "sd", coef = "Intercept", group = "Subject"),
  set_prior ("normal(0, 20)", class = "sd", coef = "Days", group = "Subject"),
  set_prior ("normal (250, 50)", class = "Intercept"),
  set_prior ("normal (10, 10)", class ="b"),
  set_prior ("normal (0, 50)", class = "sigma"))

library("brms")
fit_sleep_m3 <- brm(formula = Reaction ~ (Days|Subject) + Days,
                 data = sleepstudy,
                 family = gaussian(),
                 prior = my_priors,
                 seed = 10,
                 control = list(adapt_delta = 0.8, max_treedepth=10),
                 warmup=1000, iter=3000, thin=1, chains=3, cores=3)

summary(fit_sleep_m3)
waic(fit_sleep_m3)
plot(fit_sleep_m3)


plot(fit_sleep_m1)
