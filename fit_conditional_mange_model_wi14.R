# prepare all the data for analysis
source("prep_data.R")

# fit the model
start <- Sys.time()
mout <- run.jags(model = "./jags_script/conditional_model.R",
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
                             'n_coyote',
                             'n_mange'),
                 data = data_list,
                 n.chains = 6,
                 inits = inits,
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
#ans2 <- summary(mout2)
str(ans)

