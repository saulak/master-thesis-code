# This example is adapted from the ANOVA lecture and book material by Dr. Lukas Meier (ETH Zürich, 2024).
# It was originally implemented using lme4 and here extended with a Bayesian version using brms.
#https://doi.org/10.1201/9781003146216


################## Machines Example #####################

install.packages("nlme")
install.packages("lme4")
library(lme4)
library(dplyr)
library(ggplot2)

################## Dataset #####################

data("Machines", package = "nlme")
str(Machines)
head(Machines)

## technical details for shorter output:
class(Machines) <- "data.frame"
Machines[, "Worker"] <- factor(Machines[, "Worker"], levels = 1:6, ordered = FALSE)
str(Machines, give.attr = FALSE) ## give.attr to shorten output


# Adjust structure for cleaner output
class(Machines) <- "data.frame"
Machines[, "Worker"] <- factor(Machines[, "Worker"], levels = 1:6, ordered = FALSE)
str(Machines, give.attr = FALSE)

# Compute mean scores per Worker-Machine combination
mean_data <- Machines %>%
  group_by(Worker, Machine) %>%
  summarise(mean_score = mean(score), .groups = "drop")

# Plot scores and subject-specific means
ggplot(Machines, aes(x = Machine, y = score, group = Worker, 
                     color = as.factor(Worker), shape = as.factor(Worker))) +
  geom_point(size = 3, alpha = 0.6) +
  geom_line(data = mean_data, aes(x = Machine, y = mean_score, group = Worker), linewidth = 0.5) +
  geom_point(data = mean_data, aes(x = Machine, y = mean_score, group = Worker), size = 4) +
  theme_bw() +
  labs(color = "Worker", shape = "Worker") +
  scale_color_manual(values = c("yellow1", "deepskyblue", "lightgreen", "purple", "darkorange", "deeppink"))



############## MODEL LME ############## 

options(contrasts = c("contr.treatment", "contr.poly"))

library(lme4)
install.packages("lmerTest")
library(lmerTest)
fit_machines_lme4 <- lmer(score ~ Machine + (1 | Worker) +
                       (1 | Worker:Machine), data = Machines)
summary(fit_machines_lme4)


coef(fit_machines_lme4)
fixef(fit_machines_lme4)
ranef(fit_machines_lme4)


confint(fit_machines_lme4, oldNames = FALSE)


############## Diagnostics ############## 

plot(fit_machines_lme4)

par(mfrow = c(1, 3))
qqnorm(ranef(fit_machines_lme4)$Worker[,1], main = "Random effects of Worker")
qqnorm(ranef(fit_machines_lme4)$'Worker:Machine'[,1], main = "Random Interaction")
qqnorm(resid(fit_machines_lme4), main = "Residuals")

############## MODEL BRM ############## 

library(brms)

fit_machines_brms <- brm(
  score ~ Machine + (1 | Worker) + (1 | Worker:Machine), 
  data = Machines,
  family = gaussian(),
  prior = c(
    prior(normal(0, 10), class = "b"),                               # Prior für Fixed Effect Machine
    prior(cauchy(0, 2), class = "sd", group = "Worker"),             # Prior für Random Effect Worker
    prior(cauchy(0, 2), class = "sd", group = "Worker:Machine")      # Prior für Random Effect Interaction Worker:Machine
  ),
  control = list(adapt_delta = 0.95, max_treedepth=10),
  warmup=1000, iter=4000, thin=1, chains=3, cores=3)


summary(fit_machines_brms)

coef(fit_machines_brms)
plot(fit_machines_brms)

################## Checks #####################
pp_check(fit_machines_brms, ndraws = 100)  + theme_bw()
pp_check(fit_machines_brms, ndraws = 100) +
  theme_bw() +
  labs(x = "Productivity Score", y = "Density")


