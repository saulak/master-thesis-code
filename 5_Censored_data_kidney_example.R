# This example is adapted from:
# Bürkner, P.-C. (2017). brms: An R Package for Bayesian Multilevel Models Using Stan.
# Journal of Statistical Software, 80(1), 1–28. https://doi.org/10.18637/jss.v080.i01


################## Censored data Example #####################


################## Dataset #####################
install.packages("brms", dependencies = TRUE)
install.packages("rstan", dependencies = TRUE)
install.packages("StanHeaders", dependencies = TRUE)

library("brms")
data("kidney", package = "brms")
head(kidney)
str(kidney)
help(kidney)


################## Model 1 #####################

fit_kidney_m1 <- brm(formula = time | cens(censored) ~ age * sex + disease + (1 + age|patient),
data = kidney, family = lognormal(),
prior = c(set_prior("normal(0,5)", class = "b"),
          set_prior("cauchy(0,2)", class = "sd"),
          set_prior("lkj(2)", class = "cor")), warmup = 1500,
iter = 3000, chains = 4, control = list(adapt_delta = 0.95))

summary(fit_kidney_m1)



get_prior(time | cens(censored) ~ age * sex + disease + (1 + age | patient), data = kidney, family = lognormal())

get_prior(
  time | cens(censored) ~ age * sex + disease + (1 + age | patient),
  data = kidney,
  family = lognormal()
)


################## Interpreting Model Results #####################

stancode(fit_kidney_m1)
standata(fit_kidney_m1)

################## Effects #####################


marginal_effects(fit_kidney_m1)
conditional_effects(fit_kidney_m1)

conditional_effects(fit_kidney_m1, effects = "age")
conditional_effects(fit_kidney_m1, effects = "age", conditions = list(sex = "female", disease = "other"))
conditional_effects(fit_kidney_m1, effects = "age", conditions = list(sex = "male", disease = "other"))

conditional_effects(fit_kidney_m1, effects = "sex")
conditional_effects(fit_kidney_m1, effects = "disease")
conditional_effects(fit_kidney_m1, effects = "disease", conditions = list(sex = "male"))
conditional_effects(fit_kidney_m1, effects = "disease", conditions = list(sex = "female"))


conditional_effects(fit_kidney_m1, effects = "sex", conditions = list(disease = "other"))
conditional_effects(fit_kidney_m1, effects = "disease", conditions = list(sex = "male"))


conditional_effects(
  fit_kidney_m1,
  effects = "age:sex",
  conditions = list(disease = "other")
)


conditional_effects(
  fit_kidney_m1,
  effects = "age:sex",
  conditions = list(disease = "PKD")
)

marginal_effects(fit_kidney_m1)

################## Model 2 #####################
fit_kidney_m2 <- update(fit_kidney_m1, formula. = ~ . - (1 + age | patient) + (1 | patient))


################## LOO #####################
loo(fit_kidney_m1)
loo(fit_kidney_m2)


################## plots #####################
pdf("fit_kidney_m1_diagnostics.pdf", width = 9.5, height = 14)
plot(fit_kidney_m1, nvariables = 11)  
dev.off()
getwd()


################## Check #####################

library(bayesplot)

posterior <- as_draws_df(fit_kidney_m1)  # oder as.matrix(), oder as_draws_array()
mcmc_dens_overlay(
  posterior,
  pars = c("b_age", "b_sexfemale", "b_diseasePKD")
)

