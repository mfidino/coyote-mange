model{
#
# latent state - OPEN
# latent state - OPEN
# latent state - OPEN
#
# This associates the parameters from the psi transition matrix to the
# imperfectly observed true occupancy status of coyote.
# 1 = no coyote, 2 = coyote w/o mange, 3 = coyote w/ mange
# When indexed like so psi is a vector of length 3 that contains the 
# unnormalized probability weights for z being in each state. JAGS
# does the normalization for us, which effectively just completes
# the multinomial logit link function (i.e., psi[ site, yr, 1:nstate] are
# three numerators of the multinomial logit link and JAGS divides each element
# by sum(psi[ site, yr, 1:nstate]) to normalize.)
#
  for(site in 1:nsite){
    for(yr in 1:nyear){
      z[ site, yr ] ~ dcat(psi[ site, yr, 1:nstate ])
    }
  }
#
# latent state - CLOSE
# latent state - CLOSE
# latent state - CLOSE
#
# detection - OPEN
# detection - OPEN 
# detection - OPEN
#
# This links the latent state, z, to the actual observed data through the 
# detection matrix in the model. y is the data that we supply and is used
# to estimate the probability of detecting each state GIVEN the true imperfectly
# observed state z. We do this through by indexing the fourth dimension of the
# detection probability matrix lambda by z. Given this indexing lambda is a 
# vector of length 3 associated to the probability of observing the given
# detection.
#
  for(site in 1:nsite){
    for(yr in 1:nyear){
      for(surv in 1:nsurvey){
        y[ site, yr, surv ] ~ dcat(lambda[ site, yr,1:nstate, z[site, yr] ])
      } 
    }    
  } 
#
# detection - CLOSE
# detection - CLOSE 
# detection - CLOSE
#
# psi tpm - OPEN 
# psi tpm - OPEN
# psi tpm - OPEN
#
# This is the linear predictors for each of the 3 possible latent states. We
# have set no coyotes (state 1) as the refereence category for the multinomial
# logit link (i.e., it is set to 1). All parameters for coyote w/o mange,
# state 2, begin with the letter 'a' while those with mange, state 3,
# begin with 'b'. a and b are included in the third state because 'b' parameters
# are the relative change in occupancy given the 'a' parameters. 
# We exponentiate each of these probabilities to ensure that the elements in
# psi represent the numerators of the multinomial logit link.
# We've also incorporated nested indexing here to estimate seasonal and yearly
# effects.
  for(site in 1:nsite){
    for(yr in 1:nyear){
      pr_occ_coyote[site, yr] <- ilogit(a0 + inprod(a, acovs[ site, ]) +  
                                aseason[ sea_vec[yr] ] )
      pr_occ_mcoyote[site, yr] <- ilogit(b0 + inprod(b, bcovs[ site, yr, ]) +  
                                    bseason[ sea_vec[yr] ] )
      # no coyote, set as reference category
      psi[ site, yr, 1 ] <- 1 - pr_occ_coyote[site, yr]
      # coyote without mange occupancy linear predictor
      psi[ site, yr, 2 ] <- pr_occ_coyote[site, yr] * 
                            (1 - pr_occ_mcoyote[site, yr])
      # coyote with mange occupancy linear predictor
      psi[ site, yr, 3 ] <- pr_occ_coyote[site, yr] * pr_occ_mcoyote[site, yr]
    }
  } 
#
# psi tpm - CLOSE 
# psi tpm - CLOSE
# psi tpm - CLOSE
#
# lambda tpm - OPEN
# lambda tpm - OPEN
# lambda tpm - OPEN
#
  for(site in 1:nsite){
    for(yr in 1:nyear){
    # linear predictor for probability of detecting a coyote per week
    pr_coyote[site, yr] <- ilogit( d0 + inprod(d, dcovs[site, ]) + 
                                    dseason[ sea_vec[yr] ] )
    for(week in 1:nweek){
      # TS = no coyote = 1
      lambda[site, yr, week, 1, 1] <- 1 # OS = no coyote
      lambda[site, yr, week, 2, 1] <- 0 # OS = coyote w/o mange (not possible | TS)
      lambda[site, yr, week, 3, 1] <- 0 # OS = coyote w/ mange  (not possible | TS)
      # TS = coyote w/o mange = 2
      lambda[site, yr, week, 1, 2] <- 1 - pr_coyote[site, yr]
      lambda[site, yr, week, 2, 2] <- pr_coyote[site, yr]
      lambda[site, yr, week, 3, 2] <- 0 # OS = coyote w/ mange (not possible | TS)
      # TS = coyote w/ mange = 3
      lambda[site, yr, week, 1, 3] <- 1 - pr_coyote[site, yr]
      lambda[site, yr, week, 2, 3] <- pr_coyote[site, yr] * 
                                        (1 - pr_mangy_coyote[site, yr, week])
      lambda[site, yr, week, 3, 3] <- pr_coyote[site, yr] * 
                                        pr_mangy_coyote[site, yr, week]
    }
    }
  }
# lambda tpm - CLOSE
# lambda tpm - CLOSE
# lambda tpm - CLOSE
#
# Detection sub-model of mange per image at site and season - OPEN
# Detection sub-model of mange per image at site and season - OPEN
# Detection sub-model of mange per image at site and season - OPEN
for(photo in 1:nphoto){
  logit(mange_mu[photo]) <- f0 + inprod(f, mange_covs[photo,])
  mange_signs_present[photo] ~ dbern(mange_mu[photo])
}
# Detection sub-model of mange per image at site and season - CLOSE
# Detection sub-model of mange per image at site and season - CLOSE
# Detection sub-model of mange per image at site and season - CLOSE
#  
#
# priors - OPEN
# priors - OPEN  
# priors - OPEN
#  
# Latent, coyote w/o mange
a0 ~ dlogis(0, 1)
for(a_covs in 1:1){
  a[a_covs] ~ dlogis(0, 1)
}
#
# Latent, coyote w/ mange
b0 ~ dlogis(0, 1)
for(b_covs in 1:2){
  b[b_covs] ~ dlogis(0, 1)
}
#
# Pr(detecting coyotes)
d0 ~ dlogis(0, 1)
for(d_covs in 1:1){
  d[d_covs] ~ dlogis(0, 1)
}
#
# Pr(detecting coyotes with mange)
f0 ~ dlogis(0, 1)
for(f_covs in 1:4){
  f[f_covs] ~ dlogis(0, 1)
}
#
# Seasonal effects
aseason[1] <- 0
bseason[1] <- 0
dseason[1] <- 0
fseason[1] <- 0
for(sea in 1:3){
  aseasonpar[sea] ~ dlogis(0, 1)
  bseasonpar[sea] ~ dlogis(0, 1)
  dseasonpar[sea] ~ dlogis(0, 1)
  fseasonpar[sea] ~ dlogis(0, 1)
}
for(fillsea in 2:4){
  aseason[fillsea] <- aseasonpar[fillsea-1]
  bseason[fillsea] <- bseasonpar[fillsea-1]
  dseason[fillsea] <- dseasonpar[fillsea-1]
  fseason[fillsea] <- fseasonpar[fillsea-1]
}
#

# priors - CLOSE
# priors - CLOSE
# priors - CLOSE
}