model{
  # the mange detection model
  
  for(photo in 1:nphoto){
    logit(mange_mu[photo]) <- f0 + inprod(f, mange_covs[photo,])
    mange_signs_present[photo] ~ dbern(mange_mu[photo])
  }
  
  f0 ~ dlogis(0,1)
  for(i in 1:I){
    f[i] ~ dlogis(0,1)
  }
  #bern_trial ~ dbeta(4,4)
  #
  #for(site in 1:nsite){
  #  for(yr in 1:nyear){
  #    for(week in 1:nweek){
  #      pr_mangy_coyote[week,yr,site] <- hp[week,yr,site] * 
  #        (1 - ( exp(sum(log((1-mange_mu[(mmin[week,yr,site]):(mmax[week,yr,site])])))))) +
  #        (1 - hp[week,yr,site]) * bern_trial
  #      
  #    }}}
  
}