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

# the number of sites for this example
nsite <- 300
# species occupancy linear predictor
set.seed(-565)

# an environmental covaraite
x <- rnorm(nsite) 

# occupancy  logit-linear predictor
z_det <- 0.5 - 1 * x 
z_prob <- plogis(z_det)

# true location of species
z <- rbinom(
  nsite,
  1,
  z_prob
)

# the observed data

# species detection logit-linear predictor
x2 <- rnorm(nsite)

y_det <- -0.5 + 0.5 * x2
y_prob <- plogis(y_det) * z
# assuming 4 weeks of sampling
y <- rbinom(nsite, 4, y_prob)

# mange linear predictor assuming it's the same covariate as in occupancy
w_det <- -0.5 + 0.5 * x 
w_prob <- plogis(w_det) * z

w <- rbinom(nsite, 1, w_prob)

# mange, by photo, linear predictor

# step 1. Figure out where we got photos
sites_with_photos <- which(y>0)

# step 2. simulate number of photos per site
photos_per_site <- sample(1:30, length(sites_with_photos), replace = TRUE)
n_photos <- sum(photos_per_site)

# covaraites, one binary, one continuous
x3 <- rnorm(n_photos)
x4 <- rbinom(n_photos, 1, 0.3)


my_sites <- data.frame(sites = sites_with_photos,
                       count = photos_per_site)
true_sv <- rep(my_sites$sites, my_sites$count)
site_photos <- vector('list', length = sum(y>0))

k_det <- -0.5 + 0.25 * x3 + 0.7 * x4

y_seen <- as.numeric(y>0)
k_prob <- plogis(k_det) * as.numeric(y_seen[true_sv] * w[true_sv])

k <- rbinom(n_photos, 1, k_prob)

hm <- glm(k ~ x3 + x4, family = 'binomial')

# put together the data list that we need for this analysis

data_list <- list(y = y, 
                  q = k, 
                  psi_cov = cbind(1,x),
                  rho_cov = cbind(1, x2),
                  omega_cov = cbind(1, x),
                  gamma_cov = cbind(1, x3, x4),
                  nphoto = n_photos,
                  nsite = nsite,
                  site_vec = true_sv,
                  ncov_psi = 2,
                  ncov_rho = 2,
                  ncov_omega = 2,
                  ncov_gamma = 3,
                  J = rep(4, nsite))


x_guess <- rep(NA, nsite) 
x_guess[unique(true_sv[k == 1])] <- 1

ncov_gamma <- 3
ncov_rho <- 2
ncov_psi <- 2
ncov_omega <- 2
nseason <- 1

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
m1 <- run.jags("./jags_script/conditional_model_single_season.R",
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

