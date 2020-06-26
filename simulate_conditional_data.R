##############################################
#
# A simulated example of the occupancy/image analysis
#
#  Code written by M. Fidino
#
##############################################

# load in the necessary packages
library(runjags)
library(mcmcplots)

# The number of sites for this example
nsite <- 300

# The number of visits to each site
j <- 4

set.seed(-565)

# An environmental covaraite
x1 <- rnorm(nsite) 

# Occupancy  logit-linear predictor
z_det <- 0.5 - 1 * x1 

# Probability of occupancy
z_prob <- plogis(
  z_det
)

# True occupancy status of species
z <- rbinom(
  nsite,
  1,
  z_prob
)

####################
# the observed data
####################

# A different covariate for detection
x2 <- rnorm(
  nsite
)

# Logit-linear predictor
y_det <- -0.5 + 0.5 * x2

# The probability of detection given presence
y_prob <- plogis(
  y_det
)

# The observed data 
#  Detection probability is 0 if species is not present
y <- rbinom(
  nsite,
  j,
  y_prob * z
)

# Mange logit-linear predictor.
#  We will use the same covariate as we did with the latent
#  occupancy state.
w_det <- -0.5 + 0.5 * x1 

# the probability a coyote is mangy at the site
w_prob <- plogis(
  w_det
) 

# The mange status of coyote across the sites
w <- rbinom(
  nsite,
  1,
  w_prob * z
)

# Simulate mange detection per image

# Step 1. figure out where we collected photos
sites_with_photos <- which(
  y>0
)

# Step 2. simulate number of photos per site
photos_per_site <- sample(
  1:30,
  length(sites_with_photos),
  replace = TRUE
)

# This is the total number of images
n_photos <- sum(
  photos_per_site
)

# This is a data.frame where each row has info on the sites
#  with images and the number of images at that site.
my_sites <- data.frame(
  sites = sites_with_photos,
  count = photos_per_site
)

# This is the site index. It is of length n_photos. Each
#  element represents the site an image belongs to.
site_idx <- rep(
  my_sites$sites,
  my_sites$count
)

# Covariates, one binary, one continuous
x3 <- rnorm(
  n_photos
)

x4 <- rbinom(
  n_photos,
  1,
  0.3
)

# Logit-linear predictor for detecting mange given presence
g_det <- -0.5 + 0.25 * x3 + 0.7 * x4

# The probability of detecting mange in an image given presence
g_prob <- plogis(
  g_det
)

# A binary vector that determines which sites we have observed
#  the species. We use this and the 'w' vector to simulate
#  our observed data. Basically we need a photo at the site
#  and mangy coyote must also be at the site.
y_observed <- as.numeric(
  y>0
)

g <- rbinom(
  n_photos,
  1,
  g_prob * as.numeric(
    y_observed[site_idx] * w[site_idx]
  )
)

# This is the initial values for whether mange is at a site. It
#  equals 1 if we observed it and is NA otherwise.
x_guess <- rep(
  NA, 
  nsite
) 
x_guess[unique(site_idx[g == 1])] <- 1

# The number of parameters we are estimating
ncov_psi <- 2
ncov_rho <- 2
ncov_omega <- 2
ncov_gamma <- 3

# put together the data list that we need for this analysis
data_list <- list(
  y = y, 
  q = g, 
  psi_cov = cbind(1, x),
  rho_cov = cbind(1, x2),
  omega_cov = cbind(1, x),
  gamma_cov = cbind(1, x3, x4),
  nphoto = n_photos,
  nsite = nsite,
  site_vec = site_idx,
  ncov_psi = ncov_psi,
  ncov_rho = ncov_rho,
  ncov_omega = ncov_omega,
  ncov_gamma = ncov_gamma,
  J = rep(j,nsite)
)


# generate initial values for the analysis
inits_simulated <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = as.numeric(y_seen),
      x = x_guess,
      psi = rnorm(ncov_psi),
      rho = rnorm(ncov_rho),
      omega = rnorm(ncov_omega),
      gamma = rnorm(ncov_gamma),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,           
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}

# fit the conditional model
m1 <- run.jags(
  "./jags_script/conditional_model_single_season.R",
  monitor = c("psi", "omega","rho","gamma", "n_coyote", "n_mange"),
  data = data_list,
  n.chains = 4,
  inits = inits_simulated,
  burnin = 10000,
  adapt = 10000,
  sample = 20000,
  modules = 'glm',
  method = 'parallel'
)

# put together a plot of the parameter estimates
#  and compare them to the true values
mcmcplots::caterplot(m1, regex = "psi|rho|gam|ome", reorder = FALSE)
my_pars <- c(0.5, -1, -0.5, 0.5, -0.5, 0.5, -0.5, 0.25, 0.7)
points( rev(1:9) ~ my_pars, pch = 19)

# convert output to a matrix for plotting purposes
mm <- as.matrix(as.mcmc.list(m1), chains = TRUE)

# calculate coyote without mange
newx <- cbind(1,seq(-3,3, 0.01))

psi <- mm[,2:3]
omega <- mm[,6:7]

psi_pred <- plogis(psi %*% t(newx))
omega_pred <- plogis(omega %*% t(newx))

c_nomange <- psi_pred * (1 - omega_pred)
c_mange <- psi_pred * omega_pred

c_nomange <- t(apply(c_nomange, 2, quantile, probs = c(0.025,0.5,0.975)))
c_mange <- t(apply(c_mange, 2, quantile, probs = c(0.025,0.5,0.975)))

plot(c_nomange[,2] ~ newx[,2], xlab = "Environmental covariate",
     ylab = "Probability of...", ylim = c(0,1), type = 'l', bty = 'l',
     lwd = 2, las = 1)
lines(c_nomange[,1] ~ newx[,2], lty = 2)
lines(c_nomange[,3] ~ newx[,2], lty = 2)

lines(c_mange[,2] ~ newx[,2], col = "red", lwd = 2)
lines(c_mange[,1] ~ newx[,2], col = "red", lty = 2)
lines(c_mange[,3] ~ newx[,2], col = "red", lty = 2)

legend('topright', legend = c('Species without mange',
                              'Species with mange'),
       col = c('black', 'red'), lwd = 2, bty = "n")

