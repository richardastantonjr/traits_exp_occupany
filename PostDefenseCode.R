

## Prepare trait-based analysis of bird occupancy responses to shrub encroachment and land-use intensification

## Packages needed
##-----------------
library(jagsUI)
library(coda)

## Import and process data
##----------------------------------------------------------------------------------
## array of detection-nondetection histories 371 sites*5 visits*209 spp.
## The name is an artifact of the raw data being abundance rather than detected/not
## Also contains covariate data 
load("C:/Users/Owner/Dropbox/Manuscripts/Shrub encroachment and vertebrate communities/community str across gradients/GenerateAbunHists.RData")

## Ancillary data -- spp traits
TraitData <- read.csv("C:/Users/Owner/Dropbox/Manuscripts/Shrub encroachment and vertebrate communities/community str across gradients/data/TraitData.csv")

## subset to remove uncommon species and irrelevant/incomplete trait data
## sort alphabetically to match AbunHists.subset
TraitData <- TraitData[c(which(as.character(TraitData$Species) %in% dimnames(AbunHists.subset)[[3]])),
                       c(1,2,4,6,9:14,23:26)]
TraitData <- plyr::arrange(TraitData,Species)

## Classify species as predatory or not
predator<-rep(NA, nrow(TraitData))
for (i in 1:nrow(TraitData)){
  if((TraitData[i,11] > 0) || (TraitData[i,12] > 0) || (TraitData[i,14] > 0)){
    predator[i] = 1
  } else {
    predator[i] = 0} 
}
TraitData$predator <- predator
TraitData$pseudo_loading <- TraitData$Mass / TraitData$Wing ## proxy for wing loading


##------------------------------------------------------------------------------------
## Modify site and sampling covariates to work with JAGS/BUGS
#site-level covariates CANNOT have NAs
point<-as.numeric(factor(c(1:nrow(siteCovs))))
land<-as.numeric(factor(siteCovs$LandUse))
nland<-length(levels(siteCovs$LandUse))
shrub<-siteCovs$ShrubMean
grid = as.numeric(factor(siteCovs$grid))

#sampling-level covariates CANNOT have NAs
date=scale(SamplingCovs[,c(2:5)])
date[is.na(date)] <- 0      ## fill missing values with mean imputation
time=scale(SamplingCovs[,c(12:15)])
time[is.na(time)] <- 0
wind=SamplingCovs[,c(17:20)]
temp=scale(SamplingCovs[,c(22:25)]) 
cloud=scale(SamplingCovs[,c(27:30)])


##---------------------------------------------------------------------------------------
## fit Occupancy models for each spp.
##-------------------------------------
# Specify model in BUGS language

#for (i in c(1,2,3,4,5,6,7,9,11,12,13,14,15,16,17,20,21,22,24,25,26,27,28,31,32,33,34,35,36,37,39,40,41,42,43, 44, 46,47,48)) { 
#for (i in c(10, 18)){
#for (i in c(19, 23, 29)){
#for (i in c(45)){
#for (i in c(8)){  
#for (i in c(30, 38)){
for (i in c(44)){
  detect.i <- AbunHists.subset[,,i]
  y <- as.matrix(detect.i)
  # Bundle data
  jags.data <- list(y = y, 
                    R = length(point), 
                    J = ncol(y),
                    ngrid = length(unique(siteCovs$grid)),
                    grid = grid,
                    land = land,
                    
                    shrub = shrub,
                    date = date,
                    time = time
  )
  
  # Initial values
  zst <- apply(cbind(y, rep(0,nrow(y))), 1, max, na.rm = TRUE)  # use max of each row as starting value for z (0 / 1)
  
  inits <- function(){list(z = zst,
                           beta.l = rnorm(nland, mean = mean(y, na.rm = T)),
                           ## beta1 = mean(y, na.rm = T)
                           grid.effect = rnorm(length(unique(grid)), mean = mean(y, na.rm = T)), sd.grid = runif(1,0,5)
  )}
  
  # Parameters monitored
  params <- c("beta.l", "beta1")  
  # MCMC settings
  ni <- 190000                    ## 190000   ## 250000 for spp. 10, 18; 450000 for 19, 23, 29; 1200000 for 45; 2m for 8; 4.5m for 30, 38
  #ni <- 250000
  #ni <- 450000
  #ni <- 1200000
  #ni <- 2000000
  #ni <- 4500000
  nt <- 50                       ##     50    ## 5000 for 30, 38
  #nt <- 5000
  nb <- 40000                     ##  40000    ## 2,000,000 for 30, 38
  #nb <- 2000000
  nc <- 4
  
  # Call JAGS from R 
  out.rep <- jags.basic(jags.data, inits, params, "./code/jagsmodel1.R", n.thin = nt, n.chains = nc, n.burnin = nb, 
                  n.iter = ni, parallel = TRUE)
  allPost <- as.matrix(out.rep)
  assign(paste('Output',i,sep= ''), allPost)

  print (dim(allPost))
  save.image("./Occupancy_analysis_Swazi_birds_Bayesian20170807.RData")
} 

save.image("./Occupancy_analysis_Swazi_birds_Bayesian20170807.RData")

