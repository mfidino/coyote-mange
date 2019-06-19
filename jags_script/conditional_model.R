model{
  psi0 ~ dlogis(0,1)
  psi1 ~ dlogis(0,1)
  rho0 ~ dlogis(0,1)
  rho1 ~ dlogis(0,1)
  ome0 ~ dlogis(0,1)
  ome1 ~ dlogis(0,1)
  gam0 ~ dlogis(0,1)
  gam1 ~ dlogis(0,1)
  gam2 ~ dlogis(0,1)
  
  for(site in 1:nsite){
    logit(psi_mu[site]) <- psi0 + psi1 * x1[site]
    z[site] ~ dbern(psi_mu[site])
  }
  for(site in 1:nsite){
    logit(rho_mu[site]) <- rho0 + rho1 * x2[site]
    y[site] ~ dbin(rho_mu[site] * z[site], 4)
  }
  for(site in 1:nsite){
    logit(ome_mu[site]) <- ome0 + ome1 * x1[site]
    x[site] ~ dbern(ome_mu[site] * z[site])
  }
  for(photo in 1:nphoto){
    logit(gam_mu[photo]) <- gam0 + gam1 * x3[photo] + gam2 * x4[photo]
    q[photo] ~ dbern(gam_mu[photo] * x[site_vec[photo]])
  }
  n_coyote <- sum(z)
  n_mange <- sum(x)
}