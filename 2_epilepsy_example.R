# This example is adapted from the brms reference manual:
# Bürkner, P.-C. (2024). brms: Bayesian Regression Models using 'Stan'.
# R package version 2.22.0. https://CRAN.R-project.org/package=brms


################## Epilepsy Example #####################

library(brms)
data("epilepsy", package="brms")
str(epilepsy)
head(epilepsy, n=10)

fit1 <- brm(formula = count ~ zBase * Trt + (1 | patient),
            data = epilepsy, 
            family = poisson(),
            prior = c(set_prior("normal(0,10)", class = "b"),
                      set_prior("cauchy(0,2)", class = "sd")),
            warmup = 1000, iter = 2000, chains = 4, 
            control = list(adapt_delta = 0.95)
)

summary(fit1)
plot(fit1)


################## Checks #####################
head(predict(fit1))
fixef(fit1)
ranef(fit1)

pp_check(fit1)

plot(conditional_effects(fit1), ask = FALSE)






