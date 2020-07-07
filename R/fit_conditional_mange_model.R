# prepare all the data for analysis
source("./R/prep_data.R")

# fit the model
start <- Sys.time()
mout <- run.jags(model = "./jags_script/conditional_model_residual_variance.R",
                 monitor = c('psi',
                             'rho',
                             'omega',
                             'gamma',
                             'sd_psi',
                             'sd_rho',
                             'sd_omega',
                             'psi_ranef',
                             'rho_ranef',
                             'omega_ranef',
                             "sd_psi_resid",
                             "sd_rho_resid",
                             "sd_ome_resid",
                             "sd_gam_resid",
                             'n_coyote',
                             'n_mange'),
                 data = data_list,
                 n.chains = 6,
                 inits = inits_resid,
                 adapt = 1000,
                 burnin = 25000,
                 sample = ceiling(200000/6),
                 thin = 2,
                 module = 'glm',
                 method = 'parallel')
end <- Sys.time()
end - start


m2 <- as.mcmc.list(mout)
saveRDS(mout, "./results/coyote_mcmc_inxs_update.RDS")
ans <- summary(mout)
str(ans)

