model{
  # priors
  for(pc in 1:ncov_psi){
    psi[pc] ~ dlogis(0, 1)
  }
  for(rc in 1:ncov_rho){
    rho[rc] ~ dlogis(0,1)
  }
  for(oc in 1:ncov_omega){
    omega[oc] ~ dlogis(0,1)
  }
  for(gc in 1:ncov_gamma){
    gamma[gc] ~ dlogis(0,1)
  }
  for(ns in 1:nseason){
    psi_ranef[ns] ~ dnorm(0, tau_psi)
    rho_ranef[ns] ~ dnorm(0, tau_rho)
    omega_ranef[ns] ~ dnorm(0, tau_omega)
  }
  tau_psi ~ dgamma(1, 1)
  tau_rho ~ dgamma(1, 1)
  tau_omega ~ dgamma(1, 1)
  sd_psi <- 1 / sqrt(tau_psi)
  sd_rho <- 1 / sqrt(tau_rho)
  sd_omega <- 1 / sqrt(tau_omega)
  # observation level residual variance
  psi_resid_tau ~ dgamma(1, 1)
  rho_resid_tau ~ dgamma(1, 1)
  ome_resid_tau ~ dgamma(1, 1)
  gam_resid_tau ~ dgamma(1, 1)
  sd_psi_resid <- 1/sqrt(psi_resid_tau)
  sd_rho_resid <- 1/sqrt(rho_resid_tau)
  sd_ome_resid <- 1/sqrt(ome_resid_tau)
  sd_gam_resid <- 1/sqrt(gam_resid_tau)
  for(site in 1:nsite){
    psi_mu_logit[site] <- inprod(psi, psi_cov[site,]) + 
                                  psi_ranef[sample_vec[site]] 
    psi_mu_tmp[site] ~ dnorm(psi_mu_logit[site], psi_resid_tau)
    logit(psi_mu[site]) <- psi_mu_tmp[site]
    z[site] ~ dbern(psi_mu[site])
  }
  for(site in 1:nsite){
    rho_mu_logit[site] <- inprod(rho, rho_cov[site,]) + 
                                  rho_ranef[sample_vec[site]]
    rho_mu_tmp[site] ~ dnorm(rho_mu_logit[site], rho_resid_tau)
    logit(rho_mu[site]) <- rho_mu_tmp[site]
    y[site] ~ dbin(rho_mu[site] * z[site], J[site])
  }
  for(site in 1:nsite){
    ome_resid[site] ~ dnorm(0, ome_resid_tau)
    logit(ome_mu_logit[site]) <- inprod(omega, omega_cov[site,]) + 
                                  omega_ranef[sample_vec[site]]
    ome_mu_tmp[site] ~ dnorm(ome_mu_logit[site], ome_resid_tau)
    logit(ome_mu[site]) <- ome_mu_tmp[site]
    x[site] ~ dbern(ome_mu[site] * z[site])
  }
  for(photo in 1:nphoto){
    gam_mu_logit[photo] <- inprod(gamma, gamma_cov[photo,]) 
    gam_mu_tmp[photo] ~ dnorm(gam_mu_logit[photo], gam_resid_tau)
    logit(gam_mu[photo]) <- gam_mu_tmp[photo]
    q[photo] ~ dbern(gam_mu[photo] * x[site_vec[photo]])
  }
  # derived quantities
  for(ns in 1:nseason){
    n_coyote[ns] <- sum(z[c(st[ns,1]:st[ns,2])])
    n_mange[ns] <-  sum(x[c(st[ns,1]:st[ns,2])])
  }
}