###########################
#
# Evaluate fit of mange model
#
# Written by M. Fidino
#
###########################
library(coda)
library(dplyr)
library(tidyr)

# read in the data
source("./R/prep_data.R")

# read in the mange model output
pars <- readRDS(
  "./results/coyote_mcmc_autologistic.RDS"
)

# make it into a matrix
pars <- as.matrix(
  as.mcmc.list(pars)
)

# remove sites with no data
psi_cov <- data_list$psi_cov#[-which(data_list$J == 0),]
rho_cov <- data_list$rho_cov#[-which(data_list$J == 0),]
ome_cov <- data_list$omega_cov
gam_cov <- data_list$gamma_cov#[data_list$site_vec,]
y <- data_list$y#[-which(data_list$J == 0)]
q <- data_list$q
sample_vec <- data_list$sample_vec#[-which(data_list$J == 0)]
last_sample_vec <- data_list$last_sample_vec#[-which(data_list$J == 0)]
j <- data_list$J[-which(data_list$J == 0)]
site_vec <- data_list$site_vec

theta_psi <- as.numeric(y>0)
theta_psi[1:sum(sample_vec == 1)] <- 0

# we can only evaluate on data that we observed, I use this object a 
#  bunch to drop out some of the rows from the data in the preceeding
#  rows.
to_go <- which(data_list$J == 0)

# get the random effects and the sd term
pars_ranef <- pars[,grep("ranef", colnames(pars))]
pars_sd <- pars[,grep("sd", colnames(pars))]

# get just the parameters
pars <- pars[,-grep("ranef|coyote|mange|sd", colnames(pars))]

# this function will grab parameters from the specific linear
#  predictor we are interested in
state <- function(x,y){
  x[,grep(y, colnames(x))]
}
# inverse logit!
ilogit <- function(x) exp(x) / (1 + exp(x))

set.seed(567)

nsamps <- 5000
psims <- sample(1:nrow(pars), nsamps)

# we need to generate these by season because of the auto-logistic term in the
#  model. So this list stores the log odds of occupancy.
psi_lpred <- vector("list", data_list$nseason)
my_z <- vector("list", data_list$nseason)
 
# do first season, which does not depend on the previous state. This tmp
#  vector will index the first season of data.
tmp <- which(data_list$sample_vec == 1)
 psi_lpred[[1]] <-    state(pars, "psi")[psims,1:4] %*% t(psi_cov[tmp,]) +
    state(pars_ranef, "psi")[psims,][,sample_vec[tmp]]
 # convert to probability
psi_fs <- ilogit(psi_lpred[[1]])
 # and generate a posterior estimate of occupancy
my_z[[1]] <- matrix(
  rbinom(prod(dim(psi_lpred[[1]])), 1, psi_fs),
  ncol = data_list$nfirst_season,
  nrow = nsamps
)
# but if we actually detected the species we know it's there
my_z[[1]][,which(y[tmp] >0)] <- 1
# follow suit for the rest of the seasons
for(i in 2:data_list$nseason){
  tmp <- which(sample_vec == i)
  psi_lpred[[i]] <-    state(pars, "psi")[psims,1:4] %*% t(psi_cov[tmp,]) +
    state(pars_ranef, "psi")[psims,][,sample_vec[tmp]]
  # this adds the auto-logistic term to the linear predictor if
  #  a species was present during the previous timestep
  lz <- sweep(my_z[[i-1]], 1, state(pars, "theta")[psims,1], "*")
  psi_lpred[[i]] <- psi_lpred[[i]] + lz
  psi_fs <- ilogit(psi_lpred[[i]])
  my_z[[i]] <- matrix(
    rbinom(prod(dim(psi_lpred[[i]])), 1, psi_fs),
    ncol = data_list$nfirst_season, nrow = nsamps)
  # if it's there it is there
  my_z[[i]][,which(y[tmp] >0)] <- 1
}

# now store ALL of it in one big matrix
logit_psi <- matrix(NA, ncol = data_list$nsite, nrow = nsamps)
for(i in 1:data_list$nseason){
  logit_psi[,sample_vec == i] <- psi_lpred[[i]]
}

# drop data we did not observe
logit_psi <- logit_psi[,-to_go]

# convert to a probability
psi_pred <- ilogit(logit_psi)
rm(psi_lpred)
# Second linear predictor. rho. Much easier to calculate
#  given there is no autologistic term
rho_lpred <- state(pars, "rho")[psims,] %*% t(rho_cov) +
  state(pars_ranef, "rho")[psims,][sample_vec]

# same thing, drop data we did not observe
rho_lpred <- rho_lpred[,-to_go]

# convert to a probability
rho_pred <- ilogit(rho_lpred)
rm(rho_lpred)

# generate z given psi, following wright et al. 2019

# This is the proability a species was not detected over the
#  entire survey. 
p_notdetected <- sweep(1 - rho_pred, 2, j, FUN = "^")

# the probability that z = 1 given the species is not detected
z_prob <- (psi_pred * p_notdetected) / 
  ((1 - psi_pred) + (psi_pred * p_notdetected)) 

# estimate whether species is there
z_pred <- matrix(
  rbinom(
    prod(
      dim(psi_pred)
    ),
    1,
    z_prob
  ),
  ncol = ncol(psi_pred),
  nrow = nrow(psi_pred)
)

# clean up y now so it's only for the data we have observed.
y <- y[-to_go]
# Make it a 1 if we detected it, because we know coyote are there.
z_pred[,which(y>0)] <- 1

# This is a modified version of the function from Wright et al.'s 2019
#  Ecology paper. Changed it around so that it can be applied
#  to the entire posterior sample instead of one iteration.
#  Note: You'll want to play around with prob.seq in this function
#  to determine the correct binning size to have roughly equal
#  groups. We went with eight (i.e., 1/8).
occ_binned_df <- function(z_mat, psi_mat, cov_vec){
  # raw residuals
  resid_raw <- z_mat - psi_mat
  # the sequence of values to split on.
  prob.seq <- seq(
    0, 1, by = 1/8
  )
  # figure out where the breaks are in our covariate based on quantiles
  breaks.temp <- as.vector(
    quantile(
      cov_vec,
      prob.seq
    )
  )
  # going from those breaks to factors.
  grp.fact <- factor(
    .bincode(
      cov_vec,
      breaks.temp,
      include.lowest = TRUE
    )
  )
  # Splitting the raw residuals into groups based on grp.fact
  resid_split <- split(
    data.frame(
      t(resid_raw)
    ),
    grp.fact
  )
  # Calculate the mean residual from each group and each iteration.
  resid_split <- lapply(
    resid_split,
    colMeans
  )
  # The xvalue for the binned residual.
  xvals <- unlist(
    lapply(
      split(
        cov_vec,
        grp.fact
      ),
      mean
    )
  )
  # The mean binned residual for each iteration.
  resid_df <- data.frame(matrix(
    unlist(resid_split),
    nrow = length(resid_split[[1]]),
    ncol = length(resid_split)
  ),
  sim = 1:length(resid_split[[1]])
)
  # go from wide to long format.
  resid_long <- tidyr::pivot_longer(
    resid_df,
    cols = starts_with('X'),
    names_to = "xvals",
    values_to = "residual"
  )
  # make the xvalues numeric
  resid_long$xvals <- xvals[
    as.numeric(gsub("^X","", resid_long$xvals))
  ]
  # order by xvalue first
  resid_long <- resid_long[order(resid_long$xvals, resid_long$sim),]
 
  return(resid_long)
}

# plot out the diagnostics for psi
psi_urb1 <- occ_binned_df(
  z_pred,
  psi_pred,
  psi_cov[-to_go,2]
)

jpeg("./plots/model_evaluation/psi_urb1_residual.jpg")
plot(1 ~ 1, pch = 19, col = scales::alpha("black", 0.005),
     ylim = c(-0.35, 0.35), xlim = range(psi_urb1$xvals), bty = 'l',
     xlab = "URB1", ylab = "z residual", type = 'n')
for(i in 1:max(psi_urb1$sim)){
  hm <- psi_urb1[psi_urb1$sim == i,]
  lines(hm$residual ~ hm$xvals, col = scales::alpha("black", 0.01))
}

  abline(h = 0, lwd = 2, lty = 2, col = "gray80")

dev.off()

psi_urb2 <- occ_binned_df(z_pred, psi_pred, psi_cov[-to_go,3])

jpeg("./plots/model_evaluation/psi_urb2_residual.jpg")
plot(1 ~ 1, pch = 19, col = scales::alpha("black", 0.005),
     ylim = c(-0.35, 0.35), xlim = range(psi_urb2$xvals), bty = 'l',
     xlab = "URB2", ylab = "z residual", type = 'n')
for(i in 1:max(psi_urb2$sim)){
  hm <- psi_urb2[psi_urb2$sim == i,]
  lines(hm$residual ~ hm$xvals, col = scales::alpha("black", 0.01))
}
abline(h = 0, lwd = 2, lty = 2, col = "gray80")
#}
dev.off()

psi_urb3 <- occ_binned_df(z_pred, psi_pred, psi_cov[-to_go,4])
jpeg("./plots/model_evaluation/psi_urb1x2_residual.jpg")
plot(1 ~ 1, pch = 19, col = scales::alpha("black", 0.005),
     ylim = c(-0.35, 0.35), xlim = range(psi_urb3$xvals), bty = 'l',
     xlab = "URB1x2", ylab = "z residual", type = 'n')
for(i in 1:max(psi_urb3$sim)){
  hm <- psi_urb3[psi_urb3$sim == i,]
  lines(hm$residual ~ hm$xvals, col = scales::alpha("black", 0.01))
}
abline(h = 0, lwd = 2, lty = 2, col = "gray80")
#}
dev.off()

# Calculate detection residuals

# Another modified function from Wright et al's (2019) ecology paper.
#  Stuck with using the iteration due to the subsetting of the z_mat that
#  happens each step.
det_binned_df <- function(y_vec, p_mat, z_mat, cov_vec, iter = 1, j){
  # y is conditional on z, so only assess y if z is 1.
  keep <- which(z_mat[iter, ] == 1)
  # raw residual of y. 
  resid_raw <- apply(
    cbind(y_vec[keep], p_mat[iter,keep], j[keep]),
    1, 
    function(x){
      if(is.na(x[1])){
        return(NA)
      } else {
        ndets <- rep(0, x[3])
        if(x[1]>0){
          ndets[1:x[1]] <- 1
        }
        ndets - rep(x[2], x[3])
      }
    }
  )
  # What quantiles are we going to bin y on
  prob.seq <- seq(
    0,
    1,
    by = 1/15
  )
  # what are the quantiles of the covaraite
  breaks.temp <- as.vector(
    quantile(
      cov_vec[keep],
      prob.seq
    )
  )
  # Make it a factor
  grp.fact <- factor(
    .bincode(
      cov_vec[keep],
      breaks.temp,
      include.lowest = TRUE
    )
  )
  # Calculate mean binned residual
  yvals <- unlist(
    lapply(
      lapply(
        split(
          resid_raw,
          grp.fact
        ),
        unlist
      ),
      mean
    )
  )
  # x value for each binned residual
  xvals <- unlist(
    lapply(
      split(
        cov_vec[keep],
        grp.fact
      ),
      mean
    )
  )
  # return the residual.
  df <- data.frame(
    sim = iter,
    xvals = xvals,
    residual = yvals
  )
  return(df)
}

# Store each calculation in a list
rho_temp <- vector(
  "list",
  length = nrow(rho_pred)
)

for(i in 1:length(rho_temp)){
  rho_temp[[i]] <- det_binned_df(
    y,
    rho_pred,
    z_pred,
    rho_cov[-to_go,2],
    iter = i,
    j = j
  )
}

# Bind rows so it's on big data.frame like the z residual.
rho_temp <- dplyr::bind_rows(rho_temp)

# Plot it out
jpeg("./plots/model_evaluation/rho_temp_residual.jpg")
plot(1 ~ 1, pch = 19, col = scales::alpha("black", 0.005),
     ylim = c(-0.3, 0.3), xlim = range(rho_temp$xvals), bty = 'l',
     xlab = "temperature", ylab = "y residual", type = 'n')
for(i in 1:max(rho_temp$sim)){
  hm <- rho_temp[rho_temp$sim == i,]
  lines(hm$residual ~ hm$xvals, col = scales::alpha("black", 0.01))
}
abline(h = 0, lwd = 2, lty = 2, col = "gray80")
#}
dev.off()


# Coyote with mange calculation

# We have to condition on z for this, but also need to assess instances
#  where we missed coyote but there may be a mangy one present (and
#  therefore have no images of them). Additionally, our indexing of
#  q is going to be easiest if we do this for all of the data 
#  and then subset at the end, which means we need to recalculate
#  all of the predictions from before.

z_pred <- matrix(NA, nrow = nsamps, ncol = data_list$nsite)

for(i in 1:length(my_z)){
  tmp <- which(sample_vec == i)
  z_pred[,tmp] <- my_z[[i]]
}

# just recreating this object as a reminder of what needs to be removed
#  (i.e., times when we did not sample)
to_go <- which(data_list$J == 0)



# simulate a dataset! We need my_z for this, as detecting mange is 
#  conditional on coyote being present.



ome_lpred <- vector("list", length = data_list$nseason)
my_x <- vector("list", length = data_list$nseason)

tmp <- which(data_list$sample_vec == 1)
ome_lpred[[1]] <-    state(pars, "ome")[psims,1:4] %*% t(ome_cov[tmp,]) +
  state(pars_ranef, "ome")[psims,][,sample_vec[tmp]]
# convert to probability
ome_fs <- ilogit(ome_lpred[[1]])
# multiply by presence of coyote
ome_fs <- ome_fs * my_z[[1]]
# and generate a posterior estimate of occupancy
my_x[[1]] <- matrix(
  rbinom(prod(dim(ome_lpred[[1]])), 1, ome_fs),
  ncol = data_list$nfirst_season,
  nrow = nsamps
)
# but if we actually detected the species we know it's there

my_x[[1]][,which(x_guess[tmp] == 1)] <- 1
# follow suit for the rest of the seasons
for(i in 2:data_list$nseason){
  tmp <- which(sample_vec == i)
  ome_lpred[[i]] <-    state(pars, "ome")[psims,1:4] %*% t(ome_cov[tmp,]) +
    state(pars_ranef, "ome")[psims,][,sample_vec[tmp]]
  # this adds the auto-logistic term to the linear predictor if
  #  a species was present during the previous timestep
  lx <- sweep(my_x[[i-1]], 1, state(pars, "theta")[psims,2], "*")
  ome_lpred[[i]] <- ome_lpred[[i]] + lx
  ome_fs <- ilogit(ome_lpred[[i]])
  ome_lpred[[i]] <- ome_lpred[[i]] * my_z[[i]]
  my_x[[i]] <- matrix(
    rbinom(prod(dim(ome_lpred[[i]])), 1, psi_fs),
    ncol = data_list$nfirst_season, nrow = nsamps)
  my_x[[i]][,which(x_guess[tmp] == i)] <- 1
}

# now store ALL of it in one big matrix
logit_ome <- matrix(NA, ncol = data_list$nsite, nrow = nsamps)
for(i in 1:data_list$nseason){
  logit_ome[,sample_vec == i] <- ome_lpred[[i]]
}


# convert it to a probability
ome_pred <- ilogit(logit_ome)
#rm(ome_lpred)
# Fourth linear predictor. rho. Again, easier without the auto-logistic term.
gam_lpred <- state(pars, "gam")[psims,] %*% t(data_list$gamma_cov)

gam_pred <- ilogit(gam_lpred)
rm(gam_lpred)


# generate z given psi, following wright et al. 2019



# generate x given z, following wright et al. 2019. This is a bit trickier
#  as there are actually two seperate probabilities that we need to
#  calculate.

# 1. If we observed mange then give x_pred a 1.
# 2. If we have images but did not observe mange, then the use the
#      probability that we did not detect mange in the images.


# calculate 2. 

# we need to get the probability that mange was not detected at each
#  site given the unequal number of images per site.
gam_list <- split(
  data.frame(
    t(gam_pred)
  ),
  factor(data_list$site_vec) 
)

# we then iterate through this list. Calculate the compliment, and take
#  the row product, which represents the probability that mange was
#  not detected (given it's presence).
for(i in 1:length(gam_list)){
  gam_list[[i]] <- apply(
    1 - t(gam_list[[i]]),
    1,
    prod
  )
}

# we can then put these values in a detection matrix
# make a matrix for gam stuff
gam_mat <- matrix(NA, ncol = ncol(ome_pred), nrow = nrow(ome_pred))
for(i in 1:length(gam_list)){
  gam_mat[, as.numeric(names(gam_list))[i]] <- as.numeric(gam_list[[i]])
}

# we can then use gam_mat for calculating 2 & 3 above.

x_prob <- (ome_pred * gam_mat) / 
  ((1 - ome_pred) + (ome_pred * gam_mat)) 

# make x_prediction.

x_pred <- matrix(
  rbinom(
    prod(
      dim(x_prob)
    ),
    1,
    x_prob
  ),
  ncol = ncol(x_prob),
  nrow = nrow(x_prob)
)

# make it a 1 if we observed mange. This is the vector x_guess
x_pred[,which(x_guess == 1)] <- 1

# this is going to be conditional on z, so we'll use a similar function
#  to rho.


omega_binned_df <- function(x_mat, omega_mat, z_mat, cov_vec, iter = 1){
  # z_mat already
  #  has NA values for no sampling so it will remove the values we should
  #  not access.
  keep <- which(z_mat[iter,] == 1 & !is.na(x_mat[iter,]))
  # raw residual of y.
  resid_raw <- x_mat[iter,keep] - omega_mat[iter,keep]
  
  # What quantiles are we going to bin y on
  prob.seq <- seq(
    0,
    1,
    by = 1/8
  )
  # what are the quantiles of the covaraite
  breaks.temp <- as.vector(
    quantile(
      cov_vec[keep],
      prob.seq
    )
  )
  # Make it a factor
  grp.fact <- factor(
    .bincode(
      cov_vec[keep],
      breaks.temp,
      include.lowest = TRUE
    )
  )
  # Calculate mean binned residual
  yvals <- unlist(
    lapply(
      lapply(
        split(
          resid_raw,
          grp.fact
        ),
        unlist
      ),
      mean
    )
  )
  # x value for each binned residual
  xvals <- unlist(
    lapply(
      split(
        cov_vec[keep],
        grp.fact
      ),
      mean
    )
  )
  # return the residual.
  df <- data.frame(
    sim = iter,
    xvals = xvals,
    residual = yvals
  )
  return(df)
}


omega_1 <- vector("list", length = nrow(x_pred))

pb <- txtProgressBar(min = 1, max = length(omega_1))
for(i in 1:length(omega_1)){
  setTxtProgressBar(pb, i)
  omega_1[[i]] <- omega_binned_df(x_pred, ome_pred, z_pred, data_list$omega_cov[,2], iter = i)
}

omega_1 <- dplyr::bind_rows(omega_1)


# Plot it out
jpeg("./plots/model_evaluation/omega_urb1_given_coyote_detection.jpg")
plot(1 ~ 1, pch = 19, col = scales::alpha("black", 0.005),
     ylim = c(-0.3, 0.3), xlim = range(omega_1$xvals), bty = 'l',
     xlab = "URB 1", ylab = "x residual", type = 'n')
for(i in 1:max(omega_1$sim)){
  hm <- omega_1[omega_1$sim == i,]
  lines(hm$residual ~ hm$xvals, col = scales::alpha("black", 0.01))
}
abline(h = 0, lwd = 2, lty = 2, col = "gray80")
#}
dev.off()


omega_2 <- vector("list", length = nrow(x_pred))

pb <- txtProgressBar(min = 1, max = length(omega_2))
for(i in 1:length(omega_2)){
  setTxtProgressBar(pb, i)
  omega_2[[i]] <- omega_binned_df(x_pred, ome_pred, z_pred, data_list$omega_cov[,3], iter = i)
}

omega_2 <- dplyr::bind_rows(omega_2)


# Plot it out
jpeg("./plots/model_evaluation/omega_urb2_given_coyote_detection.jpg")
plot(1 ~ 1, pch = 19, col = scales::alpha("black", 0.005),
     ylim = c(-0.3, 0.3), xlim = range(omega_2$xvals), bty = 'l',
     xlab = "URB 2", ylab = "x residual", type = 'n')
for(i in 1:max(omega_2$sim)){
  hm <- omega_2[omega_2$sim == i,]
  lines(hm$residual ~ hm$xvals, col = scales::alpha("black", 0.01))
}
abline(h = 0, lwd = 2, lty = 2, col = "gray80")
#}
dev.off()


omega_3 <- vector("list", length = nrow(x_pred))

pb <- txtProgressBar(min = 1, max = length(omega_3))
for(i in 1:length(omega_3)){
  setTxtProgressBar(pb, i)
  omega_3[[i]] <- omega_binned_df(x_pred, ome_pred, z_pred, data_list$omega_cov[,4], iter = i)
}

omega_3 <- dplyr::bind_rows(omega_3)


# Plot it out
jpeg("./plots/model_evaluation/omega_urb1x2_given_coyote_detection.jpg")
plot(1 ~ 1, pch = 19, col = scales::alpha("black", 0.005),
     ylim = c(-0.3, 0.3), xlim = range(omega_3$xvals), bty = 'l',
     xlab = "URB 1x2", ylab = "x residual", type = 'n')
for(i in 1:max(omega_3$sim)){
  hm <- omega_3[omega_3$sim == i,]
  lines(hm$residual ~ hm$xvals, col = scales::alpha("black", 0.01))
}
abline(h = 0, lwd = 2, lty = 2, col = "gray80")
#}
dev.off()


# calculate the same as above but condition on observing coyote

# Now do the by-image mange evaluation, which will be conditional on 
#  mange presence.

gamma_binned <- function(q_vec, g_mat, x_mat, cov_vec, iter = 1, svec){
  # y is conditional on z, so only assess y if z is 1.
  tmp_x <- x_mat[iter,]
  
  tmp_x <- tmp_x[svec]
  
  keep <- which(tmp_x == 1)
  # raw residual of y. 
  resid_raw <- q_vec[keep] - gam_pred[iter,keep]

  # What quantiles are we going to bin q on
  if(all(cov_vec %in% c(0,1))){
    prob.seq <- c(0,1)
  } else {
  prob.seq <- seq(
    0,
    1,
    by = 1/4
  )
  }
  # what are the quantiles of the covaraite
  breaks.temp <- as.vector(
    quantile(
      cov_vec[keep],
      prob.seq
    )
  )
  # Make it a factor
  if(all(cov_vec %in% c(0,1))){
    grp.fact <- factor(cov_vec[keep])
  } else {
    grp.fact <- factor(
      .bincode(
        cov_vec[keep],
        breaks.temp,
        include.lowest = TRUE
      )
    )
  }
  

  # Calculate mean binned residual
  yvals <- unlist(
    lapply(
      lapply(
        split(
          resid_raw,
          grp.fact
        ),
        unlist
      ),
      mean
    )
  )
  # x value for each binned residual
  xvals <- unlist(
    lapply(
      split(
        cov_vec[keep],
        grp.fact
      ),
      mean
    )
  )
  # return the residual.
  df <- data.frame(
    sim = iter,
    xvals = xvals,
    residual = yvals
  )
  return(df)
}

gamma_1 <- vector("list", length = nrow(x_pred))

pb <- txtProgressBar(min = 1, max = length(gamma_1))
for(i in 1:length(gamma_1)){
  setTxtProgressBar(pb, i)
  gamma_1[[i]] <- gamma_binned(data_list$q, gamma_pred, x_pred, data_list$gamma_cov[,2], iter = i, data_list$site_vec)
}

gamma_1 <- dplyr::bind_rows(gamma_1)


# Plot it out
jpeg("./plots/model_evaluation/gamma_clarity.jpg")
plot(1 ~ 1, pch = 19, col = scales::alpha("black", 0.005),
     ylim = c(-0.3, 0.3), xlim = range(gamma_1$xvals), bty = 'l',
     xlab = "Clarity (1 = clear image)", ylab = "q residual", type = 'n')
for(i in 1:max(gamma_1$sim)){
  hm <- gamma_1[gamma_1$sim == i,]
  lines(hm$residual ~ hm$xvals, col = scales::alpha("black", 0.01))
}
abline(h = 0, lwd = 2, lty = 2, col = "gray80")
#}
dev.off()


gamma_2 <- vector("list", length = nrow(x_pred))

pb <- txtProgressBar(min = 1, max = length(gamma_2))
for(i in 1:length(gamma_2)){
  setTxtProgressBar(pb, i)
  gamma_2[[i]] <- gamma_binned(data_list$q, gamma_pred, x_pred, data_list$gamma_cov[,3], iter = i, data_list$site_vec)
}

gamma_2 <- dplyr::bind_rows(gamma_2)


# Plot it out
jpeg("./plots/model_evaluation/gamma_in_color.jpg")
plot(1 ~ 1, pch = 19, col = scales::alpha("black", 0.005),
     ylim = c(-0.3, 0.3), xlim = range(gamma_2$xvals), bty = 'l',
     xlab = "In color (1 = in color)", ylab = "q residual", type = 'n')
for(i in 1:max(gamma_2$sim)){
  hm <- gamma_2[gamma_2$sim == i,]
  lines(hm$residual ~ hm$xvals, col = scales::alpha("black", 0.01))
}
abline(h = 0, lwd = 2, lty = 2, col = "gray80")
#}
dev.off()


gamma_3 <- vector("list", length = nrow(x_pred))

pb <- txtProgressBar(min = 1, max = length(gamma_3))
for(i in 1:length(gamma_3)){
  setTxtProgressBar(pb, i)
  gamma_3[[i]] <- gamma_binned(data_list$q, gamma_pred, x_pred, data_list$gamma_cov[,4], iter = i, data_list$site_vec)
}

gamma_3 <- dplyr::bind_rows(gamma_3)

# there are too few values > 1.4 and so the grouped residuals
#  are not creating sufficiently equal groups. Removing them.
gamma_3 <- gamma_3[gamma_3$xvals < 1.4,]

# Plot it out
jpeg("./plots/model_evaluation/gamma_body_visible.jpg")
plot(1 ~ 1, pch = 19, col = scales::alpha("black", 0.005),
     ylim = c(-0.3, 0.3), xlim = range(gamma_3$xvals), bty = 'l',
     xlab = "Proportion body visible", ylab = "q residual", type = 'n')
for(i in 1:max(gamma_3$sim)){
  hm <- gamma_3[gamma_3$sim == i,]
  lines(hm$residual ~ hm$xvals, col = scales::alpha("black", 0.01))
}
abline(h = 0, lwd = 2, lty = 2, col = "gray80")
#}
dev.off()

