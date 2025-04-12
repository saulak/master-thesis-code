# This example is adapted and extended from:
# Pavel Chernyavskiy’s lecture series on Applied Bayesian Statistics (Fall 2020).
# https://www.youtube.com/@applied_bayesian_stats6070


################## 3D Posterior with 3 Parameters #####################

################## Dataset #####################
library(VGAM)
set.seed(01)
par(mfrow = c(1, 1))

location <- 1.664  # Lageparameter ξ
scale <- 1.143     # Skalierungsparameter ω
shape <- -3.849    # Formparameter α

Y_skew <- rskewnorm(n = 100, location = location, scale = scale, shape = shape)

par(mfrow = c(1, 2))
# Histogram
hist(Y_skew, breaks = 8, col = "lightblue", border = "black",
     main = "Histogram",
     xlab = "log(Secchi)", ylab = "Frequency")

#Q-Q-Plot
qqnorm(Y_skew)
qqline(Y_skew, col = "blue", lwd = 2)

par(mfrow = c(1, 1))

##################### Predictive Prior Simualtion #####################
xi_samp <- rnorm(1000, mean = 0, sd = 2)   
omega_samp <- rexp(1000, rate = 1)          
alpha_samp <- rnorm(1000, mean = -2, sd = 2)            

# Feed into skew-normal likelihood (log-priors)
fake_data <- rskewnorm(1000, location = xi_samp, scale = omega_samp, shape = alpha_samp)

par(mfrow = c(1, 3))  

# Histogram
hist(fake_data, xlab = "Prior log(secchi) values", main = "")

# Dichtefunktion 
plot(density(fake_data), xlab = "Prior log(secchi) values", main = "")

# QQ-Plot
#car::qqPlot(fake_data)
qqnorm(fake_data)
qqline(fake_data, col = "blue", lwd = 2)


##################### GRID APPROXIMATION #####################

# Set up grid for the location, scale, shape parameter 
xi_grid <- seq(from = -5, to = 5, by = 0.1)  
omega_grid <- seq(from = 0.1, to = 10.1, by = 0.1)  
shape_grid <- seq(from = -8, to = 2, by = 0.1)  

# Make a 3-D grid from 1-D vectors 
full_grid <- expand.grid(xi = xi_grid, omega = omega_grid, shape = shape_grid)

# Dimension of full grid
dim(full_grid)  
head(full_grid)
tail(full_grid)

##################### Evaluate priors like professor #####################
#xi_prior <- dnorm(full_grid$xi, 0, 2, log=TRUE)
#omega_prior <- dexp(full_grid$omega, 1, log=TRUE)
#shape_prior <- dnorm(full_grid$shape, -2, 2, log=TRUE)

##################### Plot priors #####################
#plot(xi_prior ~ full_grid$xi, type = 'p', col = 'red', pch = 16,
#     xlab = expression(paste('Values of parameter ', xi)),
#     ylab = "Relative plausibility")

#plot(omega_prior ~ full_grid$omega, type = 'p', col = 'green', pch = 16,
#     xlab = expression(paste('Values of parameter ', omega)),
#     ylab = "Relative plausibility")

#plot(shape_prior ~ full_grid$shape, type = 'p', col = 'blue', pch = 16,
#     xlab = expression(paste('Values of parameter ', alpha)),
#     ylab = "Relative plausibility")

#gc()  


##################### Evaluate priors, SIMPLER WAY  #####################
xi_prior <- dnorm(xi_grid, 0, 2, log=TRUE)
omega_prior <- dexp(omega_grid, 1, log=TRUE)
shape_prior <- dnorm(shape_grid, -2, 2, log=TRUE)

# garbage clean-up after plotting 1 million points
gc()  


par(mfrow = c(1, 3))
str(full_grid$xi)
length(unique(full_grid$xi))

help("plot")

##################### Plot priors, SIMPLER WAY #####################
plot(xi_prior ~ xi_grid, type = 'l', col = 'red', lwd = 2.5,
     xlab = expression(paste('Values of parameter ', xi)),
     ylab = "Relative plausibility")

plot(omega_prior ~ omega_grid, type = 'l', col = 'green', lwd = 2.5,
     xlab = expression(paste('Values of parameter ', omega)),
     ylab = "Relative plausibility")

plot(shape_prior ~ shape_grid, type = 'l', col = 'blue', lwd = 2.5,
     xlab = expression(paste('Values of parameter ', alpha)),
     ylab = "Relative plausibility")

gc()  

##################### Compute the log-likelihood #####################
#Example one combination 
dskewnorm(Y_skew,
          location = full_grid$xi[1],
          scale = full_grid$omega[1000],
          shape = full_grid$shape[1000000], log=TRUE)

sum(dskewnorm(Y_skew,
          location = full_grid$xi[1],
          scale = full_grid$omega[1000],
          shape = full_grid$shape[1000000], log=TRUE))

library(foreach)
install.packages("doParallel")
library(doParallel)
cl <- makeCluster(6)  
registerDoParallel(cl)  
log_lik <- foreach(i = 1:nrow(full_grid), .packages = c('VGAM'),
                   .combine = 'c') %dopar% {
                     sum(dskewnorm(Y_skew,
                                   location = full_grid$xi[i],
                                   scale = full_grid$omega[i],
                                   shape = full_grid$shape[i],
                                   log = TRUE))
                   }


# Stopping the parallel cluster
stopCluster(cl)  

# Saving the log likelihood for later use
saveRDS(log_lik, "log_lik.rds")  

# Loading the pre-computed log-likelihood
log_lik <- readRDS("log_lik.rds")  
summary(log_lik)  #check for possible NAs
length(log_lik)
dim(full_grid)

##################### Compute  log-posterior and posterior #####################
log_post_num <- log_lik + xi_prior + omega_prior + shape_prior

# convert from log-scale using exp() and normalize to make real distribution
post <- exp(log_post_num) / sum(exp(log_post_num))

# prove that sums to 1
sum(post)


##################### Sampling from posterior #####################
set.seed(10111985)
samp <- full_grid[sample(1:nrow(full_grid), 
                         prob = post,         #probability of inclusion = posterior
                         size = 50000, 
                         replace = TRUE), ]
summary(samp)
summary(post)


##################### Plot the approximate 4-D posterior ##################### 
install.packages("misc3d")
install.packages("rgl")
install.packages("viridis")
library(misc3d)
library(rgl)
library(viridis)

# less likely grid points in the posterior don't make it into our sample
dens <- kde3d(x=samp$xi, y=samp$omega, z=samp$shape, 
              h=c(0.3, 0.3, 0.3), n=50)
gc()  
array <- array(data=dens$d,
               dim=c(length(dens$x), length(dens$y), length(dens$z)))
image3d(array, col=magma(12))  #
box3d(col='gray', add=TRUE)
axes3d(c('x', 'y', 'z'), labels = NULL)

# Achsen beschriften
mtext3d(expression(xi), edge = "x--", line = 3, cex = 2, font = 2)     
mtext3d(expression(omega), edge = "y--", line = 3, cex = 2, font = 2)  
mtext3d(expression(alpha), edge = "z--", line = 3, cex = 2, font = 2) 

# other  colors to try: inferno(12), viridis(12), heat.colors(12)

par(mfrow=c(1,1))


##################### histogram und Beispiele ##################### 
#histogram of xi parameter averaged over the other two parameter

#Estimates from our posterior distribution
summary(samp)

#Histograms
par(mfrow = c(1, 3))
hist(samp$xi, breaks = 15, main = "", xlab = expression(xi), col = "lightyellow")
quantile_xi <- quantile(samp$xi, 0.95)  
abline(v = quantile_xi, col = "red", lwd = 2, lty = 2)  
text(quantile_xi, y = 5, labels = "95% quantile", col = "red",  adj = c(1.1, 1.5))


hist(samp$omega, breaks = 15, main = "", xlab = expression(omega), col = "lightblue")
quantile_omega <- quantile(samp$omega, 0.75)  
abline(v = quantile_omega, col = "red", lwd = 2, lty = 2) 
text(quantile_omega, y = 5, labels = "75% quantile", col = "red",  adj = c(2.5, 1.5)) 


hist(samp$shape, breaks = 10, main = "", xlab = expression(alpha), col = "honeydew")
quantile_alpha <- quantile(samp$shape, 0.5)  
abline(v = quantile_alpha, col = "red", lwd = 2, lty = 2)  
text(quantile_alpha, y = 5, labels = "50% quantile", col = "red",  adj = c(2.2, 1.5)) 





##################### 2D Plots  ##################### 
library(MASS)

# Calculate 2D density for xi and omega
dens_x_o <- kde2d(x = samp$xi, y = samp$omega, n = 50, h = c(0.3, 0.3))

# Open new graphics window 
quartz()  

#Contour plot
filled.contour(
  dens_x_o,
  color.palette = viridis::viridis,  # Farbpalette
  xlab = expression(xi),
  ylab = expression(omega),
  cex.lab = 1.6
)


# Berechnung der 2D-Dichte (xi und shape)
dens_x_s <- kde2d(x = samp$xi, y = samp$shape, n = 50, h = c(0.3, 0.3))

quartz()  

# Wiederhole den Plot
filled.contour(
  dens_x_s,
  color.palette = viridis::viridis,  
  xlab = expression(xi),
  ylab = expression(alpha),
  cex.lab = 1.6
)


# Berechnung der 2D-Dichte (omega und shape)
dens_o_s <- kde2d(x = samp$omega, y = samp$shape, n = 50, h = c(0.3, 0.3))

quartz()  # Für macOS

# Wiederhole den Plot
filled.contour(
  dens_o_s,
  color.palette = viridis::viridis,  # Farbpalette
  xlab = expression(omega),
  ylab = expression(alpha),
  cex.lab = 1.6
)


# wieder zurück Plot our priors 
dev.off()
options(device = "RStudioGD")
par(mfrow=c(1,1))
plot(1:10, main = "Test Plot in RStudio")



##################### TRACEPLOTS ##################### 
par(mfrow=c(3,1))
plot(samp$xi, type='l', ylab='Parameter value', xlab='Sample number', main=expression(xi))
plot(samp$omega, type='l', ylab='Parameter value', xlab='Sample number', main=expression(omega))
plot(samp$shape, type='l', ylab='Parameter value', xlab='Sample number', main=expression(alpha))

#not every combination of our parameter is equally likely
summary(post)

#Bayesian estimates 
summary(samp)

#Probability that alpha is less then -1.4 is 95% --> P(alpha < -1.4)
quantile(samp$shape, 0.95)

##################### derive mean and SD of the data #####################
mean <- samp$xi + (sqrt(2/pi) * samp$omega * (samp$shape/sqrt(1+samp$shape^2)))
summary(mean) 
quantile(mean, c(0.025, 0.5, 0.975))  # 95% Credible interval of MEAN log(secchi)

SD <- sqrt(samp$omega^2 * (1 - (2/pi) * (samp$shape/sqrt(1+samp$shape^2))^2))
summary(SD)  
quantile(SD, c(0.025, 0.5, 0.975)) # 95% Credible interval of SD log(secchi)

summary(samp$shape)  
quantile(samp$shape, c(0.025, 0.5, 0.975))  # 95% Credible interval of alpha

##################### Trace plots for the mean, SD, shape ##################### 
par(mfrow=c(3,1))
plot(mean, type="l", ylab="Parameter value", xlab="Sample number", main="Mean of log(secchi)")
plot(SD, type="l", ylab="Parameter value", xlab="Sample number", main="SD of log(secchi)")
plot(samp$shape, type='l', ylab="Parameter value", xlab="Sample number", main=expression(alpha))


par(mfrow=c(1,3))
hist(mean, breaks = 10, main = "", xlab = "Expected value", col = "lightblue")
hist(SD, breaks = 10, main = "", xlab = "Standard deviation", col = "lightblue")
hist(samp$shape, breaks = 10, main = "", xlab = expression(alpha), col = "lightblue")



##################### BRMS ##################### 
library(brms)
fb <- bf(data ~ 1) 

# set priors 
my_priors <- c(
  set_prior("normal(0.8, 2)", class = "Intercept"),   
  set_prior("exponential(1)", class = "sigma"),       
  set_prior("normal(0, 2)", class = "alpha")          
)

fitb <- brm(fb, data = Y_skew, family = skew_normal(),
            prior = my_priors,
            control = list(adapt_delta = 0.9, max_treedepth = 15), seed = 10111985,
            cores = 3, iter = 4000, warmup = 1000, chains = 3, thin = 1, init_r = 1)


summary(fitb)
plot(fitb)

##################### CHECK von vorher  ##################### 
mean <- samp$xi + (sqrt(2/pi) * samp$omega * (samp$shape/sqrt(1+samp$shape^2)))
SD <- sqrt(samp$omega^2 * (1 - ((2/pi) * (samp$shape/sqrt(1+samp$shape^2))^2)))
quantile(mean, c(0.025, 0.5, 0.975))
quantile(SD, c(0.025, 0.5, 0.975))
quantile(samp$shape, c(0.025, 0.5, 0.975))






















