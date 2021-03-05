# prepare all the data for analysis
source("./R/prep_data.R")

# fit the model
start <- Sys.time()
cl <- parallel::makeCluster(4)
mout <- run.jags(
  model = "./jags_script/conditional_model.R",
  monitor = c(
    'psi',
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
    'n_mange'
  ),
  data = data_list,
  n.chains = 4,
  inits = inits,
  adapt = 1000,
  burnin = 25000,
  sample = ceiling(200000/4),
  thin = 2,
  module = 'glm',
  method = 'rjparallel',
  cl = cl
)
end <- Sys.time()
end - start

parallel::stopCluster(cl)

m2 <- as.mcmc.list(mout)
saveRDS(mout, "./results/coyote_mcmc_moredata.RDS")
ans <- summary(mout)
str(ans)


diagMCMC = function( codaObject , parName=varnames(codaObject)[1] ,
                     saveName=NULL , saveType="jpg" ) {
  DBDAplColors = c("skyblue","black","royalblue","steelblue")
  openGraph(height=5,width=7)
  par( mar=0.5+c(3,4,1,0) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
       cex.lab=1.5 )
  layout(matrix(1:4,nrow=2))
  # traceplot and gelman.plot are from CODA package:
  require(coda)
  coda::traceplot( codaObject[,c(parName)] , main="" , ylab="Param. Value" ,
                   col=DBDAplColors ) 
  
  DbdaAcfPlot(codaObject,parName,plColors=DBDAplColors)
  DbdaDensPlot(codaObject,parName,plColors=DBDAplColors)
  mtext( text=parName , outer=TRUE , adj=c(0.5,0.5) , cex=2.0 )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName), type=saveType)
  }
}


