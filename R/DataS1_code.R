####Overview####
##Updating the real data analyses for Montana bats
##silver-haired bats, Lasionycteris noctivagans (LANO) and
##little brown bats, Myotis lucifugus (MYLU).
##Ecology manuscript - Identifying occupancy model inadequacies: can residuals
##separately assess detection and presence?

##08/31/2018
##Author: Wilson Wright

##R version 3.4.1 (2017-06-30)

####Setup####
##set working directory
setwd()

##Load packages
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

library(dplyr)
library(Rcpp)
library(raster)
library(rgeos)
library(rgdal)
library(compiler)
library(Kmisc)

##Load data - includes species detections, covariates, and Montana shapefile
load('DataS1_data.RData')

##Compile three Stan models
#first is for a basic model without spatial correlation
#second model has spatial correlation among det. probs - RSR approach
basic_model <- stan_model('basic_model.stan')
spat_det2 <- stan_model('spat_det_model2.stan')

##Additional functions used for doing some of the posterior predictive checks

#extract function
#extract posterior draws for different parameters from a fitted model
fit_extract <- function(fit){
  post_logit.psi <- rstan::extract(fit, 'logit_psi')$logit_psi
  post_logit.p <- rstan::extract(fit, 'logit_p')$logit_p
  
  post_z <- rstan::extract(fit, 'z')$z
  post_z.rep <- rstan::extract(fit, 'z_rep')$z_rep
  post_y.rep <- rstan::extract(fit, 'y_rep')$y_rep
  
  post_psi <- plogis(post_logit.psi)
  post_p <- plogis(post_logit.p)
  
  return(list(z = post_z,
              z.rep = post_z.rep,
              y.rep = post_y.rep,
              psi = post_psi,
              p = post_p))
}

#need a functions to create binned residual plots
occ_binned <- function(z_mat, psi_mat, cov_vec, iter = 1){
  resid_raw <- z_mat[iter, ] - psi_mat[iter, ]
  prob.seq <- seq(0, 1, by = 1/40)
  breaks.temp <- as.vector(quantile(cov_vec, prob.seq))
  grp.fact <- factor(.bincode(cov_vec, breaks.temp, include.lowest = TRUE))
  
  yvals <- unlist(lapply(split(resid_raw, grp.fact), mean))
  xvals <- unlist(lapply(split(cov_vec, grp.fact), mean))
  
  df <- data.frame(x = xvals, y = yvals)
  out <- ggplot(df, aes(x = x, y = y)) +
    geom_point(colour = 'blue') +
    geom_hline(yintercept = 0, lty = 2) +
    theme_bw()
  return(out)
}

#note that the z_mat should be expanded across the visits so each 
#z is replicated for its corresponding visits
det_binned <- function(y_vec, p_mat, z_mat, cov_vec, iter = 1){
  keep <- which(z_mat[iter, ] == 1)
  resid_raw <- y_vec[keep] - p_mat[iter, keep]
  prob.seq <- seq(0, 1, by = 1/100)
  breaks.temp <- as.vector(quantile(cov_vec[keep], prob.seq))
  grp.fact <- factor(.bincode(cov_vec[keep], breaks.temp, include.lowest = TRUE))
  
  yvals <- unlist(lapply(split(resid_raw, grp.fact), mean))
  xvals <- unlist(lapply(split(cov_vec[keep], grp.fact), mean))
  
  df <- data.frame(x = xvals, y = yvals)
  out <- ggplot(df, aes(x = x, y = y)) +
    geom_point(colour = 'blue') +
    geom_hline(yintercept = 0, lty = 2) +
    theme_bw()
  return(out)
}

#just save a data frame
occ_binned_df <- function(z_mat, psi_mat, cov_vec, iter = 1){
  resid_raw <- z_mat[iter, ] - psi_mat[iter, ]
  prob.seq <- seq(0, 1, by = 1/40)
  breaks.temp <- as.vector(quantile(cov_vec, prob.seq))
  grp.fact <- factor(.bincode(cov_vec, breaks.temp, include.lowest = TRUE))
  
  yvals <- unlist(lapply(split(resid_raw, grp.fact), mean))
  xvals <- unlist(lapply(split(cov_vec, grp.fact), mean))
  
  df <- data.frame(x = xvals, y = yvals)
  out <- df
  return(out)
}

#note that the z_mat should be expanded across the visits so each 
#z is replicated for its corresponding visits
det_binned_df <- function(y_vec, p_mat, z_mat, cov_vec, iter = 1){
  keep <- which(z_mat[iter, ] == 1)
  resid_raw <- y_vec[keep] - p_mat[iter, keep]
  prob.seq <- seq(0, 1, by = 1/100)
  breaks.temp <- as.vector(quantile(cov_vec[keep], prob.seq))
  grp.fact <- factor(.bincode(cov_vec[keep], breaks.temp, include.lowest = TRUE))
  
  yvals <- unlist(lapply(split(resid_raw, grp.fact), mean))
  xvals <- unlist(lapply(split(cov_vec[keep], grp.fact), mean))
  
  df <- data.frame(x = xvals, y = yvals)
  out <- df
  return(out)
}

#function to calculate Moran's I correlogram
correl_fun <- function(dist.mat, resid.vec, dist.incr){
  n <- length(resid.vec)
  resid.bar <- mean(resid.vec)
  resid.var <- var(resid.vec)
  resid.scale <- resid.vec - resid.bar
  resid.pairs <- outer(resid.scale, resid.scale)
  resid.pairs <- resid.pairs[lower.tri(resid.pairs)]
  
  dgrp.mat <- ceiling(dist.mat / dist.incr)
  dgrp.mat <- dgrp.mat[lower.tri(dgrp.mat)]
  
  dist.mat <- dist.mat[lower.tri(dist.mat)]
  
  moran <- tapply_(resid.pairs, dgrp.mat, mean, na.rm = TRUE) * n / (n - 1) / resid.var
  dist.means <- tapply(dist.mat, dgrp.mat, mean, na.rm = TRUE)
  return(list(dist = dist.means, moran = moran))
}

##Combine appropriate data frames and standardize some of the variables
site_all <- bind_rows(site_r1, site_r2, site_r3)
site_all <- dplyr::select(site_all, -X)

visit_all <- bind_rows(visit_r1, visit_r2, visit_r3)
visit_all <- dplyr::select(visit_all, -X, -X.1)
visit_all$jul2 <- as.vector(scale(visit_all$jul) / 2)
visit_all$tmin1.std <- as.vector(scale(visit_all$tmin1) / 2)

visit_all$prcp.ind <- ifelse(visit_all$prcp1 > 0, "prcp", "none")
visit_all$prcp.ind <- factor(visit_all$prcp.ind)

visit_all2 <- left_join(visit_all, site_all, by = 'Cell')

##Setup the spatial data and determine the neighboring sites
#aggregate the cells to the 10 km level
mt.grid10 <- aggregate(mt.grid, by = "CONUS_10KM")

#identify neighbors matrix of the entire MT grid
neighbors.mat <- gTouches(mt.grid10, byid = TRUE)
neighbors.mat[which(neighbors.mat == FALSE)] <- 0
neighbors.mat[which(neighbors.mat == TRUE)] <- 1

#figure out which of these are sampled and
#expand appropriately for the detections as well
#identify which rows in the neighbors.mat correspond to the observed sites
cells.vis <- rep(NA, nrow(site_all))
for(i in 1:nrow(site_all)){
  cells.vis[i] <- which(covs_all$Cell == site_all$Cell[i])
}

cells.vis2 <- rep(NA, nrow(visit_all))
temp <- 1
for(i in 1:nrow(site_all)){
  cells.vis2[temp:(temp+site_all$visits[i]-1)] <- cells.vis[i]
  temp <- temp + site_all$visits[i]
}

neighbors1 <- neighbors.mat[cells.vis, cells.vis]
neighbors2 <- neighbors.mat[cells.vis2, cells.vis2]

for(i in 1:length(cells.vis2)){
  for(j in 1:length(cells.vis2)){
    if(cells.vis2[i] == cells.vis2[j]){
      neighbors2[i, j] <- 1
    }
  }
}
diag(neighbors2) <- 0

#save this information in terms of edges and nodes
#first for sites
N.edges1 = sum(neighbors1) / 2

node1a <- node1b <- rep(NA, N.edges1)

pos1 <- 1
for(i in 1:(nrow(site_all) - 1)){
  for(j in i:nrow(site_all)){
    if(neighbors1[i, j] == 1){
      node1a[pos1] <- i
      node1b[pos1] <- j
      pos1 <- pos1 + 1
    }
  }
}

#next for detections
N.edges2 <- sum(neighbors2) / 2
node2a <- node2b <- rep(NA, N.edges2)

pos2 <- 1
for(i in 1:(nrow(visit_all) - 1)){
  for(j in i:nrow(visit_all)){
    if(neighbors2[i, j] == 1){
      node2a[pos2] <- i
      node2b[pos2] <- j
      pos2 <- pos2 + 1
    }
  }
}

####Basic model fit####
#setup data for Stan and fit models for each species
N <- nrow(site_all)
J <- site_all$visits
NJ <- sum(J)

occ_mm <- model.matrix(~log.fpc2 + rugg + elev + ddays, site_all)
det_mm <- model.matrix(~method + jul2 + I(jul2^2) + tmin1.std + prcp.ind,
                       visit_all2)

K <- ncol(occ_mm)
L <- ncol(det_mm)

dets.lano <- visit_all$Lano
dets.mylu <- visit_all$Mylu

dets.lano[which(dets.lano > 1)] <- 1
dets.mylu[which(dets.mylu > 1)] <- 1

site_ind <- rep.int(1:N, times = J)

data_lano <- list('N' = N,
                  'J' = J,
                  'NJ' = NJ,
                  'occ_mat' = occ_mm,
                  'det_mat' = det_mm,
                  'K' = K,
                  'L' = L,
                  'dets' = dets.lano,
                  'site_ind' = site_ind,
                  'N_edges1' = N.edges1,
                  'node1a' = node1a,
                  'node1b' = node1b,
                  'N_edges2' = N.edges2,
                  'node2a' = node2a,
                  'node2b' = node2b)

fit_lano1 <- sampling(basic_model, data_lano)

data_mylu <- list('N' = N,
                  'J' = J,
                  'NJ' = NJ,
                  'occ_mat' = occ_mm,
                  'det_mat' = det_mm,
                  'K' = K,
                  'L' = L,
                  'dets' = dets.mylu,
                  'site_ind' = site_ind,
                  'N_edges1' = N.edges1,
                  'node1a' = node1a,
                  'node1b' = node1b,
                  'N_edges2' = N.edges2,
                  'node2a' = node2a,
                  'node2b' = node2b)

fit_mylu1 <- sampling(basic_model, data_mylu)

#summarize the results
summary(fit_lano1, pars = c('alphas', 'betas'), use_cache = FALSE)$summary
summary(fit_mylu1, pars = c('alphas', 'betas'), use_cache = FALSE)$summary

extract_lano1 <- fit_extract(fit_lano1)
extract_mylu1 <- fit_extract(fit_mylu1)

#creating some graphical assessments with these residuals
#setting up the distance based information needed
mt.grid10_proj <- spTransform(mt.grid10,
                              '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
names(mt.grid10_proj) <- 'Cell'

grid_keep1 <- rep(NA, nrow(site_all))
for(i in 1:length(grid_keep1)){
  grid_keep1[i] <- which(mt.grid10_proj@data$Cell %in% site_all$Cell[i])
}

grid_keep2 <- rep(NA, nrow(visit_all2))
for(i in 1:length(grid_keep2)){
  grid_keep2[i] <- which(mt.grid10_proj@data$Cell %in% visit_all2$Cell[i])
}

grid_coords1 <- coordinates(mt.grid10_proj)[grid_keep1, ]
grid_coords2 <- coordinates(mt.grid10_proj)[grid_keep2, ]

distance.mat1 <- matrix(0, ncol = 406, nrow = 406)
distance.mat2 <- matrix(0, ncol = 1455, nrow = 1455)

for(i in 1:405){
  for(j in (i + 1):406){
    distance.mat1[j, i] <- ncf::gcdist(grid_coords1[i, 1], grid_coords1[i, 2],
                                       grid_coords1[j, 1], grid_coords1[j, 2])
    distance.mat1[i, j] <- distance.mat1[j, i]
  }
}

for(i in 1:1454){
  for(j in (i + 1):1455){
    distance.mat2[j, i] <- ncf::gcdist(grid_coords2[i, 1], grid_coords2[i, 2],
                                       grid_coords2[j, 1], grid_coords2[j, 2])
    distance.mat2[i, j] <- distance.mat2[j, i]
  }
}

#occupancy
lano_corr_occ1 <- mylu_corr_occ1 <- array(NA, c(10, 2, 100))
for(i in 1:100){
  lano_temp <- correl_fun(distance.mat1, (extract_lano1$z[i, ] - extract_lano1$psi[i, ]), 15)
  lano_corr_occ1[, 1, i] <- lano_temp$dist[1:10]
  lano_corr_occ1[, 2, i] <- lano_temp$moran[1:10]
  
  mylu_temp <- correl_fun(distance.mat1, (extract_mylu1$z[i, ] - extract_mylu1$psi[i, ]), 15)
  mylu_corr_occ1[, 1, i] <- mylu_temp$dist[1:10]
  mylu_corr_occ1[, 2, i] <- mylu_temp$moran[1:10]
}

lano_corrm1 <- apply(lano_corr_occ1, 2, rbind)
lano_corrm1 <- as.data.frame(lano_corrm1)
lano_corrm1$iter <- rep(1:100, each = 10)
names(lano_corrm1) <- c('x', 'y', 'iter')

ggplot(lano_corrm1, aes(x = x, y = y, group = iter)) +
  geom_line()

mylu_corrm1 <- apply(mylu_corr_occ1, 2, rbind)
mylu_corrm1 <- as.data.frame(mylu_corrm1)
mylu_corrm1$iter <- rep(1:100, each = 10)
names(mylu_corrm1) <- c('x', 'y', 'iter')

ggplot(mylu_corrm1, aes(x = x, y = y, group = iter)) +
  geom_line()

#detection
lano_corr_det1 <- mylu_corr_det1 <- array(NA, c(11, 2, 100))
for(i in 1:100){
  lano_kt <- which(extract_lano1$z[i, site_ind] == 1)
  lano_temp <- correl_fun(distance.mat2[lano_kt, lano_kt], (dets.lano[lano_kt] - extract_lano1$p[i, lano_kt]), 15)
  lano_corr_det1[, 1, i] <- lano_temp$dist[1:11]
  lano_corr_det1[, 2, i] <- lano_temp$moran[1:11]
  
  mylu_kt <- which(extract_mylu1$z[i, site_ind] == 1)
  mylu_temp <- correl_fun(distance.mat2[mylu_kt, mylu_kt], (dets.mylu[mylu_kt] - extract_mylu1$p[i, mylu_kt]), 15)
  mylu_corr_det1[, 1, i] <- mylu_temp$dist[1:11]
  mylu_corr_det1[, 2, i] <- mylu_temp$moran[1:11]
}

lano_corrd1 <- apply(lano_corr_det1, 2, rbind)
lano_corrd1 <- as.data.frame(lano_corrd1)
lano_corrd1$iter <- rep(1:100, each = 11)
names(lano_corrd1) <- c('x', 'y', 'iter')

ggplot(lano_corrd1, aes(x = x, y = y, group = iter)) +
  geom_line()

mylu_corrd1 <- apply(mylu_corr_det1, 2, rbind)
mylu_corrd1 <- as.data.frame(mylu_corrd1)
mylu_corrd1$iter <- rep(1:100, each = 11)
names(mylu_corrd1) <- c('x', 'y', 'iter')

ggplot(mylu_corrd1, aes(x = x, y = y, group = iter)) +
  geom_line()

#some posterior predictive checks
#based on Moran's I
mean(rstan::extract(fit_lano1, 'occ_moran1')$occ_moran1 >
       rstan::extract(fit_lano1, 'occ_moran1_rep')$occ_moran1_rep)
mean(rstan::extract(fit_lano1, 'det_moran1')$det_moran1 >
       rstan::extract(fit_lano1, 'det_moran1_rep')$det_moran1_rep)

mean(rstan::extract(fit_mylu1, 'occ_moran1')$occ_moran1 >
       rstan::extract(fit_mylu1, 'occ_moran1_rep')$occ_moran1_rep)
mean(rstan::extract(fit_mylu1, 'det_moran1')$det_moran1 >
       rstan::extract(fit_mylu1, 'det_moran1_rep')$det_moran1_rep)

####RSR model fit####
##Setup the matrices used for this approach
#first need to remove some cells that have a few missing values
site.drop <- which(is.na(covs_all$ddays) == TRUE)
n.neighbors <- rowSums(neighbors.mat)[-site.drop]
N.total <- length(n.neighbors)
D.mat <- diag(n.neighbors,
              ncol = N.total,
              nrow = N.total)
Q.mat <- D.mat - neighbors.mat[-site.drop, -site.drop]

##find K matrix based on that occupancy model with covariates from all sites
pred_mm <- model.matrix(~log.fpc2 + rugg + elev + ddays, covs_all)

P.t <- diag(1, ncol = N.total, nrow = N.total) - 
  (pred_mm %*% solve(t(pred_mm) %*% pred_mm) %*% t(pred_mm) )

omega.mat <- as.vector(N.total /
                         (matrix(1, ncol = N.total, nrow=1) %*%
                            neighbors.mat[-site.drop, -site.drop] %*%
                            matrix(1, ncol = 1, nrow = N.total))) *
  (P.t %*% neighbors.mat[-site.drop, -site.drop] %*% P.t)

omega.sd <- eigen(omega.mat)

#keep the first 200 eigenvectors for the K matrix
k.keep <- which(omega.sd$values > 0)[1:400]
K.mat <- omega.sd$vectors[, k.keep]

#calculate the sigma.mat beforehand
Sigma.mat <- t(K.mat) %*% Q.mat %*% K.mat

##Setup data for Stan and fit models
data_lano3 <- list('N' = N,
                   'J' = J,
                   'NJ' = NJ,
                   'occ_mat' = occ_mm,
                   'det_mat' = det_mm,
                   'K' = K,
                   'L' = L,
                   'dets' = dets.lano,
                   'site_ind' = site_ind,
                   'N_edges1' = N.edges1,
                   'node1a' = node1a,
                   'node1b' = node1b,
                   'N_edges2' = N.edges2,
                   'node2a' = node2a,
                   'node2b' = node2b,
                   'n_alpha2' = 400,
                   'K_mat' = K.mat[cells.vis, ],
                   'Sigma_mat' = Sigma.mat)

fit_lano3 <- sampling(spat_det2, data_lano3,
                      control = list(adapt_delta = 0.9))

summary(fit_lano3, pars = c('alphas', 'betas', 'sigma_eta'),
        use_cache = FALSE)$summary

traceplot(fit_lano3, pars = c('alphas'))
traceplot(fit_lano3, pars = c('betas'))

data_mylu3 <- list('N' = N,
                   'J' = J,
                   'NJ' = NJ,
                   'occ_mat' = occ_mm,
                   'det_mat' = det_mm,
                   'K' = K,
                   'L' = L,
                   'dets' = dets.mylu,
                   'site_ind' = site_ind,
                   'N_edges1' = N.edges1,
                   'node1a' = node1a,
                   'node1b' = node1b,
                   'N_edges2' = N.edges2,
                   'node2a' = node2a,
                   'node2b' = node2b,
                   'n_alpha2' = 400,
                   'K_mat' = K.mat[cells.vis, ],
                   'Sigma_mat' = Sigma.mat)

fit_mylu3 <- sampling(spat_det2, data_mylu3,
                      control = list(adapt_delta = 0.9))

summary(fit_mylu3, pars = c('alphas', 'betas', 'sigma_eta'),
        use_cache = FALSE)$summary

traceplot(fit_mylu3, pars = c('alphas'))
traceplot(fit_mylu3, pars = c('betas'))

extract_lano3 <- fit_extract(fit_lano3)
extract_mylu3 <- fit_extract(fit_mylu3)

#occupancy
lano_corr_occ3 <- mylu_corr_occ3 <- array(NA, c(10, 2, 100))
for(i in 1:100){
  lano_temp <- correl_fun(distance.mat1, (extract_lano3$z[i, ] - extract_lano3$psi[i, ]), 15)
  lano_corr_occ3[, 1, i] <- lano_temp$dist[1:10]
  lano_corr_occ3[, 2, i] <- lano_temp$moran[1:10]
  
  mylu_temp <- correl_fun(distance.mat1, (extract_mylu3$z[i, ] - extract_mylu3$psi[i, ]), 15)
  mylu_corr_occ3[, 1, i] <- mylu_temp$dist[1:10]
  mylu_corr_occ3[, 2, i] <- mylu_temp$moran[1:10]
}

lano_corrm3 <- apply(lano_corr_occ3, 2, rbind)
lano_corrm3 <- as.data.frame(lano_corrm3)
lano_corrm3$iter <- rep(1:100, each = 10)
names(lano_corrm3) <- c('x', 'y', 'iter')

ggplot(lano_corrm3, aes(x = x, y = y, group = iter)) +
  geom_line()

mylu_corrm3 <- apply(mylu_corr_occ3, 2, rbind)
mylu_corrm3 <- as.data.frame(mylu_corrm3)
mylu_corrm3$iter <- rep(1:100, each = 10)
names(mylu_corrm3) <- c('x', 'y', 'iter')

ggplot(mylu_corrm3, aes(x = x, y = y, group = iter)) +
  geom_line()

#detection
lano_corr_det3 <- mylu_corr_det3 <- array(NA, c(11, 2, 100))
for(i in 1:100){
  lano_kt <- which(extract_lano3$z[i, site_ind] == 1)
  lano_temp <- correl_fun(distance.mat2[lano_kt, lano_kt], (dets.lano[lano_kt] - extract_lano3$p[i, lano_kt]), 15)
  lano_corr_det3[, 1, i] <- lano_temp$dist[1:11]
  lano_corr_det3[, 2, i] <- lano_temp$moran[1:11]
  
  mylu_kt <- which(extract_mylu3$z[i, site_ind] == 1)
  mylu_temp <- correl_fun(distance.mat2[mylu_kt, mylu_kt], (dets.mylu[mylu_kt] - extract_mylu3$p[i, mylu_kt]), 15)
  mylu_corr_det3[, 1, i] <- mylu_temp$dist[1:11]
  mylu_corr_det3[, 2, i] <- mylu_temp$moran[1:11]
}

lano_corrd3 <- apply(lano_corr_det3, 2, rbind)
lano_corrd3 <- as.data.frame(lano_corrd3)
lano_corrd3$iter <- rep(1:100, each = 11)
names(lano_corrd3) <- c('x', 'y', 'iter')

ggplot(lano_corrd3, aes(x = x, y = y, group = iter)) +
  geom_line()

mylu_corrd3 <- apply(mylu_corr_det3, 2, rbind)
mylu_corrd3 <- as.data.frame(mylu_corrd3)
mylu_corrd3$iter <- rep(1:100, each = 11)
names(mylu_corrd3) <- c('x', 'y', 'iter')

ggplot(mylu_corrd3, aes(x = x, y = y, group = iter)) +
  geom_line()

#some posterior predictive checks
#based on Moran's I
mean(rstan::extract(fit_lano3, 'occ_moran1')$occ_moran1 >
       rstan::extract(fit_lano3, 'occ_moran1_rep')$occ_moran1_rep)
mean(rstan::extract(fit_lano3, 'det_moran1')$det_moran1 >
       rstan::extract(fit_lano3, 'det_moran1_rep')$det_moran1_rep)

mean(rstan::extract(fit_mylu3, 'occ_moran1')$occ_moran1 >
       rstan::extract(fit_mylu3, 'occ_moran1_rep')$occ_moran1_rep)
mean(rstan::extract(fit_mylu3, 'det_moran1')$det_moran1 >
       rstan::extract(fit_mylu3, 'det_moran1_rep')$det_moran1_rep)

#some comparison plots
plot(fit_lano1, pars = c('alphas', 'betas')) + xlim(-7.5, 7.5)
plot(fit_lano3, pars = c('alphas', 'betas')) + xlim(-7.5, 7.5)

plot(fit_mylu1, pars = c('alphas', 'betas')) + xlim(-3, 4)
plot(fit_mylu3, pars = c('alphas', 'betas')) + xlim(-3, 4)

#residual plots
with(extract_lano1, occ_binned(z.rep, psi, site_all$log.fpc1, 2))

with(extract_lano1, det_binned(dets.lano, p, z[site_ind, site_ind], visit_all2$rugg, 1))
with(extract_lano1, det_binned(y.rep[1, ], p, z[site_ind, site_ind], visit_all2$rugg, 1))

with(extract_mylu1, det_binned(dets.mylu, p, z[site_ind, site_ind], visit_all2$jul, 2)) + ylim(-0.5, 0.5)
with(extract_mylu1, det_binned(y.rep[2, ], p, z[site_ind, site_ind], visit_all2$jul, 2)) + ylim(-0.5, 0.5)

with(extract_lano3, det_binned(dets.lano, p, z[site_ind, site_ind], visit_all2$jul, 2)) + ylim(-0.4, 0.4)
with(extract_lano3, det_binned(y.rep[2, ], p, z[site_ind, site_ind], visit_all2$jul, 2)) + ylim(-0.4, 0.4)

with(extract_mylu3, det_binned(dets.mylu, p, z[site_ind, site_ind], visit_all2$jul, 3)) + ylim(-0.35, 0.35)
with(extract_mylu3, det_binned(y.rep[3, ], p, z[site_ind, site_ind], visit_all2$jul, 3)) + ylim(-0.35, 0.35)
