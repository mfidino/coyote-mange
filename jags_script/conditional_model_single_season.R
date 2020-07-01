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
  for(site in 1:nsite){
    logit(psi_mu[site]) <- inprod(psi, psi_cov[site,])
    z[site] ~ dbern(psi_mu[site])
  }
  for(site in 1:nsite){
    logit(rho_mu[site]) <- inprod(rho, rho_cov[site,])
    y[site] ~ dbin(rho_mu[site] * z[site], J[site])
  }
  for(site in 1:nsite){
    logit(ome_mu[site]) <- inprod(omega, omega_cov[site,])
    x[site] ~ dbern(ome_mu[site] * z[site])
  }
  for(photo in 1:nphoto){
    logit(gam_mu[photo]) <- inprod(gamma, gamma_cov[photo,])
    q[photo] ~ dbern(gam_mu[photo] * x[site_vec[photo]])
  }
  # derived quantities

    n_coyote <- sum(z)
    n_mange <-  sum(x)

}