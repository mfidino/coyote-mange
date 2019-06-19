
nsite <- 300
# coyote occupancy linear predictor
set.seed(-565)
x <- rnorm(nsite) 

z_det <- 0.5 - 1 * x 
z_prob <- plogis(z_det)

z <- rbinom(nsite, 1, z_prob)

# the observed data

# coyote detection linear predictor
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



data_list <- list(y = y, 
                  q = k, 
                  x1 = x, 
                  x2 = x2,
                  x3 = x3, 
                  x4 = x4,
                  nphoto = n_photos,
                  nsite = nsite,
                  site_vec = true_sv)


x_guess <- rep(NA, nsite) 
x_guess[unique(true_sv[k == 1])] <- 1

inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      psi0 = rnorm(1),
      psi1 = rnorm(1),
      rho0 = rnorm(1),
      rho1 = rnorm(1),
      ome0 = rnorm(1),
      ome1 = rnorm(1),
      gam0 = rnorm(1, -1, 0.1),
      gam1 = rnorm(1),
      gam2 = rnorm(1),
      z = y_seen,
      x = x_guess,
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



library(runjags)

m1 <- run.jags("conditional_model.R",
               monitor = c("psi0", "psi1", "rho0", "rho1", "ome0", "ome1",
                           "gam0", "gam1", "gam2", "n_coyote", "n_mange"),
               data = data_list,
               n.chains = 4,
               inits = inits,
               burnin = 10000,
               adapt = 10000,
               sample = 20000,
               modules = 'glm',
               method = 'parallel'
              )

caterplot(m1, regex = "psi|rho|gam|ome", reorder = FALSE)
my_pars <- c(0.5, -1, -0.5, 0.5, -0.5, 0.5, -0.5, 0.25, 0.7)
points( rev(1:9) ~ my_pars, pch = 19)

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

plot(c_nomange[,2] ~ newx[,2], xlab = "environmental covariate",
     ylab = "probability of", ylim = c(0,1), type = 'l', bty = 'l',
     lwd = 2)
lines(c_nomange[,1] ~ newx[,2], lty = 2)
lines(c_nomange[,3] ~ newx[,2], lty = 2)

lines(c_mange[,2] ~ newx[,2], col = "red", lwd = 2)
lines(c_mange[,1] ~ newx[,2], col = "red", lty = 2)
lines(c_mange[,3] ~ newx[,2], col = "red", lty = 2)

legend('topright', legend = c('Coyote without mange',
                              'Coyote with mange'),
       col = c('black', 'red'), lwd = 2)
