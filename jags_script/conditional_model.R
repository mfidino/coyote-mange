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
  
  for(site in 1:nsite){
    logit(psi_mu[site]) <- inprod(psi, psi_cov[site,]) + 
                                  psi_ranef[sample_vec[site]]
    z[site] ~ dbern(psi_mu[site])
  }
  for(site in 1:nsite){
    logit(rho_mu[site]) <- inprod(rho, rho_cov[site,]) + 
                                  rho_ranef[sample_vec[site]]
    y[site] ~ dbin(rho_mu[site] * z[site], J[site])
  }
  for(site in 1:nsite){
    logit(ome_mu[site]) <- inprod(omega, omega_cov[site,]) + 
                                  omega_ranef[sample_vec[site]]
    x[site] ~ dbern(ome_mu[site] * z[site])
  }
  for(photo in 1:nphoto){
    logit(gam_mu[photo]) <- inprod(gamma, gamma_cov[photo,])
    q[photo] ~ dbern(gam_mu[photo] * x[site_vec[photo]])
  }
  # derived quantities
  for(ns in 1:nseason){
    n_coyote[ns] <- sum(z[c(st[ns,1]:st[ns,2])])
    n_mange[ns] <-  sum(x[c(st[ns,1]:st[ns,2])])
  }
}