
## libraries
#library(plotrix)
library(tidyverse)
library(reshape2) 
library(gridExtra)

## functions
#return mean and 95% CRI
meanAndCRI <- function(x) {
  c(mean <- mean(x),quantile(x,c(0.025, 0.975))) 
}

#return mean and 90% CRI
meanAnd90CRI <- function(x) {
  c(mean <- mean(x),quantile(x,c(0.05, 0.95))) 
}


## load model output, organized as matrices named "Output1," "Output2," etc. and most covariate data
load("./Occupancy_analysis_Swazi_birds_Bayesian20170807.RData")

## Put output from jagsUI in a list- each corresponds to a single species
Outputs <- list(Output1,Output2,Output3,Output4,Output5,Output6,Output7,Output8,
                                       Output9,Output10,Output11,Output12,
                                       Output13,Output14,Output15,Output16,
                                       Output17,Output18,Output19,Output20,
                                       Output21,Output22,Output23,Output24,
                                       Output25,Output26,Output27,Output28,
                                       Output29,Output30,Output31,Output32,
                                       Output33,Output34,Output35,Output36,
                                       Output37,Output38,Output39,Output40,
                                       Output41,Output42,Output43,Output44,
                                       Output45,Output46,Output47,Output48)

## take a random draw of ~2000 iterations from each posterior to get equal sample sizes
## create vectors of numbers for each "Output"
## subset each output to include samples "%in%" their respective vectors
to_keep <- matrix(nrow = 2000, ncol = 48)
for(i in 1:48){
to_keep[,i] <- sample(1:dim(Outputs[[i]])[1], size = 2000)
Outputs[[i]] <- subset(Outputs[[i]][c(to_keep[,i]),])   ## keep row numbers that are in column i of to_keep
}

## Assemble outputs from each species in a dataframe in long format
Outputs_long <- as.data.frame (do.call(rbind, Outputs))
Outputs_long <- Outputs_long[,1:5]   ## drop the deviance since it is not needed
## create a column for species name for downstream processing
sppName <- vector(mode = "character", length = dim(Outputs[[1]])[1])  
sppName <- rep(dimnames(AbunHists.subset)[[3]],  dim(Outputs[[1]])[1])
sppName <- sort(sppName)
Outputs_long$Species <- sppName 
Outputs_long$X <- rep(1:dim(Outputs[[1]])[1],48)   ## add a row number for dcast to convert long to wide format

## remove extraneous material
rm(list=ls()[! ls() %in% c("Outputs_long","TraitData")])
##----------------------------------------------------------------------------------
## summarize results for the community of common species
pooledBetas <- colMeans(Outputs_long[1:5])  
pooledHomesteadCRI <-  quantile(Outputs_long[,1], probs = c(0.025,0.975)) 
pooledHomesteadCRI
pooledPastureCRI <-  quantile(Outputs_long[,2], probs = c(0.025,0.975)) 
pooledPastureCRI
pooledProtectedCRI <-  quantile(Outputs_long[,3], probs = c(0.025,0.975)) 
pooledProtectedCRI
pooledSugarEstateCRI <-  quantile(Outputs_long[,4], probs = c(0.025,0.975)) 
pooledSugarEstateCRI 
pooledShrubCRI <-  quantile(Outputs_long[,5], probs = c(0.025,0.975)) 
pooledShrubCRI

#create vectors of species by trait
## nesting substrates
nest.shrub <- subset(TraitData, ShrubNest == 1, select=c(Species))
nest.cavity <- subset(TraitData, CavityNest == 1, select=c(Species))
nest.tree <- subset(TraitData, TreeNest == 1, select=c(Species))
nest.grass <- subset(TraitData, GrassNest == 1, select=c(Species))
## diets
predators <- subset(TraitData, predator == 1, select=c(Species))
diet.invert <- subset(TraitData, Diet == "invertebrates", select=c(Species))
diet.fruit <- subset(TraitData, Diet == "fruit", select=c(Species))
diet.seed <- subset(TraitData, Diet == "seeds", select=c(Species))
diet.nectar <- subset(TraitData, Diet == "nectar", select=c(Species))

## merge posteriors and species traits
Outputs_long <-merge(Outputs_long, TraitData)
#subset by effect size or beta, keeping posterior draw, species beta:
shrub.post<-subset(Outputs_long, select=c("X", "Species","beta1"))
homestead.post<-subset(Outputs_long, select=c("X", "Species","beta.l[1]"))
pasture.post<-subset(Outputs_long, select=c("X", "Species","beta.l[2]"))
pasture.post<-subset(Outputs_long, select=c("X", "Species","beta.l[3]"))
sugarEstate.post<-subset(Outputs_long, select=c("X", "Species","beta.l[4]"))

#reformat to make species columns and posterior samples rows with dcast
shrub.post <- dcast(shrub.post, X~Species)
homestead.post <- dcast(homestead.post, X~Species)
pasture.post <- dcast(pasture.post, X~Species)
protected.post <- dcast(protected.post, X~Species)
sugarEstate.post <- dcast(sugarEstate.post, X~Species)

#take rowMeans based on species list
shrub.post$nest.shrub <- rowMeans(shrub.post[,names(shrub.post) %in% nest.shrub$Species])
shrub.post$nest.tree <- rowMeans(shrub.post[,names(shrub.post) %in% nest.tree$Species])
shrub.post$nest.cavity <- rowMeans(shrub.post[,names(shrub.post) %in% nest.cavity$Species])
shrub.post$nest.grass <- rowMeans(shrub.post[,names(shrub.post) %in% nest.grass$Species])

shrub.post$predators <- rowMeans(shrub.post[,names(shrub.post) %in% predators$Species])
shrub.post$diet.invert <- rowMeans(shrub.post[,names(shrub.post) %in% diet.invert$Species])
shrub.post$diet.fruit <- rowMeans(shrub.post[,names(shrub.post) %in% diet.fruit$Species])
shrub.post$diet.seeds <- rowMeans(shrub.post[,names(shrub.post) %in% diet.seed$Species])
shrub.post$diet.nectar <- rowMeans(shrub.post[,names(shrub.post) %in% diet.nectar$Species])

protected.post$nest.shrub <- rowMeans(protected.post[,names(protected.post) %in% nest.shrub$Species])
protected.post$nest.tree <- rowMeans(protected.post[,names(protected.post) %in% nest.tree$Species])
protected.post$nest.cavity <- rowMeans(protected.post[,names(protected.post) %in% nest.cavity$Species])
protected.post$nest.grass <- rowMeans(protected.post[,names(protected.post) %in% nest.grass$Species])

protected.post$predators <- rowMeans(protected.post[,names(protected.post) %in% predators$Species])
protected.post$diet.invert <- rowMeans(protected.post[,names(protected.post) %in% diet.invert$Species])
protected.post$diet.fruit <- rowMeans(protected.post[,names(protected.post) %in% diet.fruit$Species])
protected.post$diet.seeds <- rowMeans(protected.post[,names(protected.post) %in% diet.seed$Species])
protected.post$diet.nectar <- rowMeans(protected.post[,names(protected.post) %in% diet.nectar$Species])

pasture.post$nest.shrub <- rowMeans(pasture.post[,names(pasture.post) %in% nest.shrub$Species])
pasture.post$nest.tree <- rowMeans(pasture.post[,names(pasture.post) %in% nest.tree$Species])
pasture.post$nest.cavity <- rowMeans(pasture.post[,names(pasture.post) %in% nest.cavity$Species])
pasture.post$nest.grass <- rowMeans(pasture.post[,names(pasture.post) %in% nest.grass$Species])

pasture.post$predators <- rowMeans(pasture.post[,names(pasture.post) %in% predators$Species])
pasture.post$diet.invert <- rowMeans(pasture.post[,names(pasture.post) %in% diet.invert$Species])
pasture.post$diet.fruit <- rowMeans(pasture.post[,names(pasture.post) %in% diet.fruit$Species])
pasture.post$diet.seeds <- rowMeans(pasture.post[,names(pasture.post) %in% diet.seed$Species])
pasture.post$diet.nectar <- rowMeans(pasture.post[,names(pasture.post) %in% diet.nectar$Species])

homestead.post$nest.shrub <- rowMeans(homestead.post[,names(homestead.post) %in% nest.shrub$Species])
homestead.post$nest.tree <- rowMeans(homestead.post[,names(homestead.post) %in% nest.tree$Species])
homestead.post$nest.cavity <- rowMeans(homestead.post[,names(homestead.post) %in% nest.cavity$Species])
homestead.post$nest.grass <- rowMeans(homestead.post[,names(homestead.post) %in% nest.grass$Species])

homestead.post$predators <- rowMeans(homestead.post[,names(homestead.post) %in% predators$Species])
homestead.post$diet.invert <- rowMeans(homestead.post[,names(homestead.post) %in% diet.invert$Species])
homestead.post$diet.fruit <- rowMeans(homestead.post[,names(homestead.post) %in% diet.fruit$Species])
homestead.post$diet.seeds <- rowMeans(homestead.post[,names(homestead.post) %in% diet.seed$Species])
homestead.post$diet.nectar <- rowMeans(homestead.post[,names(homestead.post) %in% diet.nectar$Species])

sugarEstate.post$nest.shrub <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% nest.shrub$Species])
sugarEstate.post$nest.tree <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% nest.tree$Species])
sugarEstate.post$nest.cavity <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% nest.cavity$Species])
sugarEstate.post$nest.grass <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% nest.grass$Species])

sugarEstate.post$predators <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% predators$Species])
sugarEstate.post$diet.invert <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% diet.invert$Species])
sugarEstate.post$diet.fruit <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% diet.fruit$Species])
sugarEstate.post$diet.seeds <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% diet.seed$Species])
sugarEstate.post$diet.nectar <- rowMeans(sugarEstate.post[,names(sugarEstate.post) %in% diet.nectar$Species])


## create a data frame of effect sizes and CRIs for each species for plotting
## species names
Species <- levels(as.factor(Outputs_long$Species))

## create a data frame of one effect and CRIs for all species
shrubSppEfx <- apply(shrub.post, 2, meanAndCRI )
shrubSppEfx <- as.data.frame(t(shrubSppEfx[,2:49]))
shrubSppEfx <- cbind(Species, shrubSppEfx)
colnames(shrubSppEfx) <- c("Species", "effect", "LCL", "UCL")
shrubSppEfx$betaName <- "Shrub cover"
rownames(shrubSppEfx)<- NULL

## repeat for the remaining effects
homesteadSppEfx <- apply(homestead.post, 2, meanAndCRI )
homesteadSppEfx <- as.data.frame(t(homesteadSppEfx[,2:49]))
homesteadSppEfx <- cbind(Species, homesteadSppEfx)
colnames(homesteadSppEfx) <- c("Species", "effect", "LCL", "UCL")
homesteadSppEfx$betaName <- "Homestead"
rownames(homesteadSppEfx)<- NULL

pastureSppEfx <- apply(pasture.post, 2, meanAndCRI )
pastureSppEfx <- as.data.frame(t(pastureSppEfx[,2:49]))
pastureSppEfx <- cbind(Species, pastureSppEfx)
colnames(pastureSppEfx) <- c("Species", "effect", "LCL", "UCL")
pastureSppEfx$betaName <- "Pasture"
rownames(pastureSppEfx)<- NULL

protectedSppEfx <- apply(protected.post, 2, meanAndCRI )
protectedSppEfx <- as.data.frame(t(protectedSppEfx[,2:49]))
protectedSppEfx <- cbind(Species, protectedSppEfx)
colnames(protectedSppEfx) <- c("Species", "effect", "LCL", "UCL")
protectedSppEfx$betaName <- "Protected"
rownames(protectedSppEfx)<- NULL

sugarSppEfx <- apply(sugarEstate.post, 2, meanAndCRI )
sugarSppEfx <- as.data.frame(t(sugarSppEfx[,2:49]))
sugarSppEfx <- cbind(Species, sugarSppEfx)
colnames(sugarSppEfx) <- c("Species", "effect", "LCL", "UCL")
sugarSppEfx$betaName <- "Plantation"
rownames(sugarSppEfx)<- NULL

## bind rows to create one frame for plotting
sppEfx <- rbind.data.frame(shrubSppEfx, protectedSppEfx, pastureSppEfx, homesteadSppEfx, sugarSppEfx )

##-------------------------------------
## build data frames with nest and diet summary effects for plotting
## Nest substrates
r1 <- meanAndCRI(shrub.post$nest.shrub)
r2 <- meanAndCRI(shrub.post$nest.tree)
r3 <- meanAndCRI(shrub.post$nest.cavity)
r4 <- meanAndCRI(shrub.post$nest.grass)
shrub.post.nest <- rbind.data.frame(r1,r2,r3,r4)
colnames(shrub.post.nest) <- c("effect", "LCL", "UCL")
shrub.post.nest$betaName <- "Shrub cover"
shrub.post.nest$nest <- c("shrubs", "trees", "cavities","grass")

r5 <- meanAndCRI (protected.post$nest.shrub)
r6 <- meanAndCRI (protected.post$nest.tree)
r7 <- meanAndCRI (protected.post$nest.cavity)
r8 <- meanAndCRI (protected.post$nest.grass)
protected.post.nest <- rbind.data.frame(r5,r6,r7,r8)
colnames(protected.post.nest) <- c("effect", "LCL","UCL")
protected.post.nest$betaName <- "Protected"
protected.post.nest$nest <- c("shrubs", "trees", "cavities","grass")

r9 <- meanAndCRI (pasture.post$nest.shrub)
r10 <- meanAndCRI (pasture.post$nest.tree)
r11 <- meanAndCRI (pasture.post$nest.cavity)
r12<- meanAndCRI (pasture.post$nest.grass)
pasture.post.nest <- rbind.data.frame(r9,r10,r11,r12)
colnames(pasture.post.nest) <- c("effect", "LCL","UCL")
pasture.post.nest$betaName <- "Pasture"
pasture.post.nest$nest <- c("shrubs", "trees", "cavities","grass")

r13 <- meanAndCRI (homestead.post$nest.shrub)
r14 <- meanAndCRI (homestead.post$nest.tree)
r15 <- meanAndCRI (homestead.post$nest.cavity)
r16 <- meanAndCRI (homestead.post$nest.grass)
homestead.post.nest <- rbind.data.frame(r13,r14,r15,r16)
colnames(homestead.post.nest) <- c("effect", "LCL","UCL")
homestead.post.nest$betaName <- "Homestead"
homestead.post.nest$nest <- c("shrubs", "trees", "cavities","grass")

r17 <- meanAndCRI (sugarEstate.post$nest.shrub)
r18 <- meanAndCRI (sugarEstate.post$nest.tree)
r19 <- meanAndCRI (sugarEstate.post$nest.cavity)
r20 <- meanAndCRI (sugarEstate.post$nest.grass)
sugarEstate.post.nest <- rbind.data.frame(r17,r18,r19,r20)
colnames(sugarEstate.post.nest) <- c("effect", "LCL","UCL")
sugarEstate.post.nest$betaName <- "Plantation"
sugarEstate.post.nest$nest <- c("shrubs", "trees", "cavities","grass")


## Diets
r21 <- meanAndCRI(shrub.post$predators)
r22 <- meanAndCRI(shrub.post$diet.invert)
r23 <- meanAndCRI(shrub.post$diet.fruit)
r24 <- meanAndCRI(shrub.post$diet.seeds)
r25 <- meanAndCRI(shrub.post$diet.nectar)
shrub.post.diet <- rbind.data.frame(r21,r22,r23,r24,r25)
colnames(shrub.post.diet) <- c("effect", "LCL","UCL")
shrub.post.diet$betaName <- "Shrub cover"
shrub.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r26 <- meanAndCRI (protected.post$predators)
r27 <- meanAndCRI (protected.post$diet.invert)
r28 <- meanAndCRI (protected.post$diet.fruit)
r29 <- meanAndCRI (protected.post$diet.seeds)
r30 <- meanAndCRI (protected.post$diet.nectar)
protected.post.diet <- rbind.data.frame(r26,r27,r28,r29,r30)
colnames(protected.post.diet) <- c("effect", "LCL","UCL")
protected.post.diet$betaName <- "Protected"
protected.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r31 <- meanAndCRI (pasture.post$predators)
r32 <- meanAndCRI (pasture.post$diet.invert)
r33 <- meanAndCRI (pasture.post$diet.fruit)
r34 <- meanAndCRI (pasture.post$diet.seeds)
r35 <- meanAndCRI (pasture.post$diet.nectar)
pasture.post.diet <- rbind.data.frame(r31,r32,r33,r34,r35)
colnames(pasture.post.diet) <- c("effect", "LCL","UCL")
pasture.post.diet$betaName <- "Pasture"
pasture.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r36 <- meanAndCRI (homestead.post$predators)
r37 <- meanAndCRI (homestead.post$diet.invert)
r38 <- meanAndCRI (homestead.post$diet.fruit)
r39 <- meanAndCRI (homestead.post$diet.seeds)
r40 <- meanAndCRI (homestead.post$diet.nectar)
homestead.post.diet <- rbind.data.frame(r36,r37,r38,r39,r40)
colnames(homestead.post.diet) <- c("effect", "LCL","UCL")
homestead.post.diet$betaName <- "Homestead"
homestead.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

r41 <- meanAndCRI (sugarEstate.post$predators)
r42 <- meanAndCRI (sugarEstate.post$diet.invert)
r43 <- meanAndCRI (sugarEstate.post$diet.fruit)
r44 <- meanAndCRI (sugarEstate.post$diet.seeds)
r45 <- meanAndCRI (sugarEstate.post$diet.nectar)
sugarEstate.post.diet <- rbind.data.frame(r41,r42,r42,r44,r45)
colnames(sugarEstate.post.diet) <- c("effect", "LCL","UCL")
sugarEstate.post.diet$betaName <- "Plantation"
sugarEstate.post.diet$diet <- c("predatory","invertebrates","fruit","seeds","nectar")

nest.summary <- rbind.data.frame(shrub.post.nest, protected.post.nest, pasture.post.nest, homestead.post.nest,
                                 sugarEstate.post.nest, stringsAsFactors = FALSE)
diet.summary <- rbind.data.frame(shrub.post.diet, protected.post.diet, pasture.post.diet ,homestead.post.diet,
                                 sugarEstate.post.diet, stringsAsFactors = FALSE)


#### Assess effects of continuous covariates using a regression approach
#### Fit a linear model to each row, where rows are posterior samples of a covariate's value and each
#### column is a species. Take the slopes, use to get CRIs and for plotting

## Shrub cover responses
## wing chord explains shrub cover occupancy response; linear and quadratic forms
## moderate positive effect of wing chord on shrub cover effect size
wingSq <- TraitData$Wing * TraitData$Wing
wing_chord_slopes <- vector()
wing_chord_quad_slopes <- vector()
pvals_temp <- vector()
for(i in 1:dim(shrub.post)[1]) {
  wing_chord_quad_slopes[i] <- summary(lm( as.numeric(shrub.post[i,2:49]) ~ Wing + wingSq, data = TraitData))$coefficients[3,1]
  wing_chord_slopes[i] <- summary(lm( as.numeric(shrub.post[i,2:49]) ~ TraitData$Wing))$coefficients[2,1]
  pvals_temp[i] <- summary(lm( as.numeric(shrub.post[i,2:49]) ~ Wing, data = TraitData))$coefficients[2,4]
  }
meanAndCRI(wing_chord_slopes)
meanAndCRI(wing_chord_quad_slopes)
meanAndCRI(pvals_temp)

## 90% CRI excludes 0
meanAnd90CRI(wing_chord_slopes)
meanAnd90CRI(wing_chord_quad_slopes)

## mass explains shrub cover occupancy response; linear and quadratic forms
massSq <- TraitData$Mass * TraitData$Mass
mass_slopes <- vector()
mass_quad_slopes <- vector()
for(i in 1:dim(shrub.post)[1]) {
  mass_quad_slopes[i] <- summary(lm( as.numeric(shrub.post[i,2:49]) ~ Mass + massSq, data = TraitData))$coefficients[3,1]
  mass_slopes[i] <- summary(lm( as.numeric(shrub.post[i,2:49]) ~ TraitData$Mass))$coefficients[2,1]
  pvals_temp[i] <- summary(lm( as.numeric(shrub.post[i,2:49]) ~ Mass, data = TraitData))$coefficients[2,4]
  }
meanAndCRI(mass_slopes)
meanAndCRI(mass_quad_slopes)
meanAndCRI(pvals_temp)

## positive effect but effect is essentially 0
meanAnd90CRI(mass_slopes)
meanAnd90CRI(mass_quad_slopes)

## wing loading proxy #(wing chord / mass)# explains shrub cover occupancy response; linear and quadratic forms
pseudo_load_Sq <- TraitData$pseudo_loading * TraitData$pseudo_loading
pseudo_load_slopes <- vector()
pseudo_load_quad_slopes <- vector()
for(i in 1:dim(shrub.post)[1]) {
  pseudo_load_quad_slopes[i] <- summary(lm( as.numeric(shrub.post[i,2:49]) ~ pseudo_loading + pseudo_load_Sq, data = TraitData))$coefficients[3,1]
  pseudo_load_slopes[i] <- summary(lm( as.numeric(shrub.post[i,2:49]) ~ TraitData$pseudo_loading))$coefficients[2,1]
  pvals_temp[i] <- summary(lm( as.numeric(shrub.post[i,2:49]) ~ pseudo_loading, data = TraitData))$coefficients[2,4]
  }
meanAndCRI(pseudo_load_slopes)
meanAndCRI(pseudo_load_quad_slopes)
meanAndCRI(pvals_temp)

meanAnd90CRI(pseudo_load_slopes)
meanAnd90CRI(pseudo_load_quad_slopes)

## Land-use responses
##-----------------------------------------------------------------------------------------
## wing chord explains land use (protection) occupancy response; linear and quadratic forms
## Effect of protection increases with wing chord
wing_chord_slopes_protected <- vector()
wing_chord_quad_slopes_protected <- vector()
for(i in 1:dim(protected.post)[1]) {
  wing_chord_quad_slopes_protected[i] <- summary(lm( as.numeric(protected.post[i,2:49]) ~ Wing + wingSq, data = TraitData))$coefficients[3,1]
  wing_chord_slopes_protected[i] <- summary(lm( as.numeric(protected.post[i,2:49]) ~ TraitData$Wing))$coefficients[2,1]
  pvals_temp[i] <- summary(lm( as.numeric(protected.post[i,2:49]) ~ Wing, data = TraitData))$coefficients[2,4]
}
meanAndCRI(wing_chord_slopes_protected)
meanAndCRI(wing_chord_quad_slopes_protected)
meanAndCRI(pvals_temp)

meanAnd90CRI(wing_chord_slopes_protected)
meanAnd90CRI(wing_chord_quad_slopes_protected)

## mass explains land use (protection) occupancy response; linear and quadratic forms
mass_slopes_protected <- vector()
mass_quad_slopes_protected <- vector()
for(i in 1:dim(protected.post)[1]) {
  mass_quad_slopes_protected[i] <- summary(lm( as.numeric(protected.post[i,2:49]) ~ Mass + massSq, data = TraitData))$coefficients[3,1]
  mass_slopes_protected[i] <- summary(lm( as.numeric(protected.post[i,2:49]) ~ TraitData$Mass))$coefficients[2,1]
  pvals_temp[i] <- summary(lm( as.numeric(protected.post[i,2:49]) ~ Mass, data = TraitData))$coefficients[2,4]
}
meanAndCRI(mass_slopes_protected)
meanAndCRI(mass_quad_slopes_protected)
meanAndCRI(pvals_temp)

## at 90%, mass has a positive effect
meanAnd90CRI(mass_slopes_protected)
meanAnd90CRI(mass_quad_slopes_protected)

## wing loading proxy #(wing chord / mass)# explains land use (protection) occupancy response; linear and quadratic forms
pseudo_load_slopes_protected <- vector()
pseudo_load_quad_slopes_protected <- vector()
for(i in 1:dim(protected.post)[1]) {
  pseudo_load_quad_slopes_protected[i] <- summary(lm( as.numeric(protected.post[i,2:49]) ~ pseudo_loading + pseudo_load_Sq, data = TraitData))$coefficients[3,1]
  pseudo_load_slopes_protected[i] <- summary(lm( as.numeric(protected.post[i,2:49]) ~ TraitData$pseudo_loading))$coefficients[2,1]
  pvals_temp[i] <- summary(lm( as.numeric(protected.post[i,2:49]) ~ pseudo_loading, data = TraitData))$coefficients[2,4]
}
meanAndCRI(pseudo_load_slopes_protected)
meanAndCRI(pseudo_load_quad_slopes_protected)
meanAndCRI(pvals_temp)

meanAnd90CRI(pseudo_load_slopes_protected)
meanAnd90CRI(pseudo_load_quad_slopes_protected)

## wing chord explains land use (pasture) occupancy response; linear and quadratic forms
wing_chord_slopes_pasture <- vector()
wing_chord_quad_slopes_pasture <- vector()
for(i in 1:dim(pasture.post)[1]) {
  wing_chord_quad_slopes_pasture[i] <- summary(lm( as.numeric(pasture.post[i,2:49]) ~ Wing + wingSq, data = TraitData))$coefficients[3,1]
  wing_chord_slopes_pasture[i] <- summary(lm( as.numeric(pasture.post[i,2:49]) ~ TraitData$Wing))$coefficients[2,1]
  pvals_temp[i] <- summary(lm( as.numeric(pasture.post[i,2:49]) ~ Wing, data = TraitData))$coefficients[2,4]
}
meanAndCRI(wing_chord_slopes_pasture)
meanAndCRI(wing_chord_quad_slopes_pasture)
meanAndCRI(pvals_temp)

meanAnd90CRI(wing_chord_slopes_pasture)
meanAnd90CRI(wing_chord_quad_slopes_pasture)

## mass explains land use (pasture) occupancy response; linear and quadratic forms
mass_slopes_pasture <- vector()
mass_quad_slopes_pasture <- vector()
for(i in 1:dim(pasture.post)[1]) {
  mass_quad_slopes_pasture[i] <- summary(lm( as.numeric(pasture.post[i,2:49]) ~ Mass + massSq, data = TraitData))$coefficients[3,1]
  mass_slopes_pasture[i] <- summary(lm( as.numeric(pasture.post[i,2:49]) ~ TraitData$Mass))$coefficients[2,1]
  pvals_temp[i] <- summary(lm( as.numeric(pasture.post[i,2:49]) ~ Mass, data = TraitData))$coefficients[2,4]
}
meanAndCRI(mass_slopes_pasture)
meanAndCRI(mass_quad_slopes_pasture)
meanAndCRI(pvals_temp)

meanAnd90CRI(mass_slopes_pasture)
meanAnd90CRI(mass_quad_slopes_pasture)

## wing loading proxy #(wing chord / mass)# explains land use (pasture) occupancy response; linear and quadratic forms
pseudo_load_slopes_pasture <- vector()
pseudo_load_quad_slopes_pasture <- vector()
for(i in 1:dim(pasture.post)[1]) {
  pseudo_load_quad_slopes_pasture[i] <- summary(lm( as.numeric(pasture.post[i,2:49]) ~ pseudo_loading + pseudo_load_Sq, data = TraitData))$coefficients[3,1]
  pseudo_load_slopes_pasture[i] <- summary(lm( as.numeric(pasture.post[i,2:49]) ~ TraitData$pseudo_loading))$coefficients[2,1]
  pvals_temp[i] <- summary(lm( as.numeric(pasture.post[i,2:49]) ~ pseudo_loading, data = TraitData))$coefficients[2,4]
}
meanAndCRI(pseudo_load_slopes_pasture)
meanAndCRI(pseudo_load_quad_slopes_pasture)
meanAndCRI(pvals_temp)

meanAnd90CRI(pseudo_load_slopes_pasture)
meanAnd90CRI(pseudo_load_quad_slopes_pasture)

## wing chord explains land use (homestead) occupancy response; linear and quadratic forms
wing_chord_slopes_homestead <- vector()
wing_chord_quad_slopes_homestead <- vector()
for(i in 1:dim(homestead.post)[1]) {
  wing_chord_quad_slopes_homestead[i] <- summary(lm( as.numeric(homestead.post[i,2:49]) ~ Wing + wingSq, data = TraitData))$coefficients[3,1]
  wing_chord_slopes_homestead[i] <- summary(lm( as.numeric(homestead.post[i,2:49]) ~ TraitData$Wing))$coefficients[2,1]
  pvals_temp[i] <- summary(lm( as.numeric(homestead.post[i,2:49]) ~ Wing, data = TraitData))$coefficients[2,4]
}
meanAndCRI(wing_chord_slopes_homestead)
meanAndCRI(wing_chord_quad_slopes_homestead)
meanAndCRI(pvals_temp)

meanAnd90CRI(wing_chord_slopes_homestead)
meanAnd90CRI(wing_chord_quad_slopes_homestead)

## mass explains land use (homestead) occupancy response; linear and quadratic forms
mass_slopes_homestead <- vector()
mass_quad_slopes_homestead <- vector()
for(i in 1:dim(homestead.post)[1]) {
  mass_quad_slopes_homestead[i] <- summary(lm( as.numeric(homestead.post[i,2:49]) ~ Mass + massSq, data = TraitData))$coefficients[3,1]
  mass_slopes_homestead[i] <- summary(lm( as.numeric(homestead.post[i,2:49]) ~ TraitData$Mass))$coefficients[2,1]
  pvals_temp[i] <- summary(lm( as.numeric(homestead.post[i,2:49]) ~ Mass, data = TraitData))$coefficients[2,4]
}
meanAndCRI(mass_slopes_homestead)
meanAndCRI(mass_quad_slopes_homestead)
meanAndCRI(pvals_temp)

meanAnd90CRI(mass_slopes_homestead)
meanAnd90CRI(mass_quad_slopes_homestead)

## wing loading proxy #(wing chord / mass)# explains land use (homestead) occupancy response; linear and quadratic forms
pseudo_load_slopes_homestead <- vector()
pseudo_load_quad_slopes_homestead <- vector()
for(i in 1:dim(homestead.post)[1]) {
  pseudo_load_quad_slopes_homestead[i] <- summary(lm( as.numeric(homestead.post[i,2:49]) ~ pseudo_loading + pseudo_load_Sq, data = TraitData))$coefficients[3,1]
  pseudo_load_slopes_homestead[i] <- summary(lm( as.numeric(homestead.post[i,2:49]) ~ TraitData$pseudo_loading))$coefficients[2,1]
  pvals_temp[i] <- summary(lm( as.numeric(homestead.post[i,2:49]) ~ pseudo_loading, data = TraitData))$coefficients[2,4]
}
meanAndCRI(pseudo_load_slopes_homestead)
meanAndCRI(pseudo_load_quad_slopes_homestead)
meanAndCRI(pvals_temp)

meanAnd90CRI(pseudo_load_slopes_homestead)
meanAnd90CRI(pseudo_load_quad_slopes_homestead)

## wing chord explains land use (sugarEstate) occupancy response; linear and quadratic forms
## minuscule quadratic effect of wing_chord on response to sugarEstate
wing_chord_slopes_sugarEstate <- vector()
wing_chord_quad_slopes_sugarEstate <- vector()
for(i in 1:dim(sugarEstate.post)[1]) {
  wing_chord_quad_slopes_sugarEstate[i] <- summary(lm( as.numeric(sugarEstate.post[i,2:49]) ~ Wing + wingSq, data = TraitData))$coefficients[3,1]
  wing_chord_slopes_sugarEstate[i] <- summary(lm( as.numeric(sugarEstate.post[i,2:49]) ~ TraitData$Wing))$coefficients[2,1]
  pvals_temp[i] <- summary(lm( as.numeric(sugarEstate.post[i,2:49]) ~ Wing + wingSq, data = TraitData))$coefficients[3,4]
}
meanAndCRI(wing_chord_slopes_sugarEstate)
meanAndCRI(wing_chord_quad_slopes_sugarEstate)
meanAndCRI(pvals_temp)

meanAnd90CRI(wing_chord_slopes_sugarEstate)
meanAnd90CRI(wing_chord_quad_slopes_sugarEstate)

## mass explains land use (sugarEstate) occupancy response; linear and quadratic forms
## Miniscule quadratic effect of mass on plantation effect size
mass_slopes_sugarEstate <- vector()
mass_quad_slopes_sugarEstate <- vector()
for(i in 1:dim(sugarEstate.post)[1]) {
  mass_quad_slopes_sugarEstate[i] <- summary(lm( as.numeric(sugarEstate.post[i,2:49]) ~ Mass + massSq, data = TraitData))$coefficients[3,1]
  mass_slopes_sugarEstate[i] <- summary(lm( as.numeric(sugarEstate.post[i,2:49]) ~ TraitData$Mass))$coefficients[2,1]
  pvals_temp[i] <- summary(lm( as.numeric(sugarEstate.post[i,2:49]) ~ Mass + massSq, data = TraitData))$coefficients[3,4]
}
meanAndCRI(mass_slopes_sugarEstate)
meanAndCRI(mass_quad_slopes_sugarEstate)
meanAndCRI(pvals_temp)

meanAnd90CRI(mass_slopes_sugarEstate)
meanAnd90CRI(mass_quad_slopes_sugarEstate)

## wing loading proxy #(wing chord / mass)# explains land use (sugarEstate) occupancy response; linear and quadratic forms
## Massive quadratic effect of pseudo wing chord on plantation effect size
pseudo_load_slopes_sugarEstate <- vector()
pseudo_load_quad_slopes_sugarEstate <- vector()
for(i in 1:dim(sugarEstate.post)[1]) {
  pseudo_load_quad_slopes_sugarEstate[i] <- summary(lm( as.numeric(sugarEstate.post[i,2:49]) ~ pseudo_loading + pseudo_load_Sq, data = TraitData))$coefficients[3,1]
  pseudo_load_slopes_sugarEstate[i] <- summary(lm( as.numeric(sugarEstate.post[i,2:49]) ~ TraitData$pseudo_loading))$coefficients[2,1]
  pvals_temp[i] <-summary(lm( as.numeric(sugarEstate.post[i,2:49]) ~ pseudo_loading + pseudo_load_Sq, data = TraitData))$coefficients[3,4]
}
meanAndCRI(pseudo_load_slopes_sugarEstate)
meanAndCRI(pseudo_load_quad_slopes_sugarEstate)
meanAndCRI(pvals_temp)

meanAnd90CRI(pseudo_load_slopes_sugarEstate)
meanAnd90CRI(pseudo_load_quad_slopes_sugarEstate)

####-----------------------------------------------------------------------------------
###------------------------------------------------------------------------------------

###                    FIGURES
##
###-----------------------------------------------------------------------------------

## community wide effects, 95% CRIs
par(mfrow = c(2,3))
for(j in 2:6){           ### corresponds to column numbers for each covariate
  plot(density(Outputs_long[,j])) 
  abline(v = quantile(Outputs_long[,j], probs = 0.025), col = "blue" )
  abline(v = quantile(Outputs_long[,j], probs = 0.975), col= "blue" )
}

## community wide effects, 90% CRIs
for(j in 2:6){           
  plot(density(Outputs_long[,j])) 
  abline(v = quantile(Outputs_long[,j], probs = 0.05), col = "green" )
  abline(v = quantile(Outputs_long[,j], probs = 0.95), col= "green" )
}

#some graph code from Isabel's paper
#Fig 4: species-specific forest plots with multiple effect sizes, categorical variables (land use)
ggplot(data = sppEfx[49:240,], aes(x = Species, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4) +                     
  #geom_hline(aes(yintercept = 0))+                ### keep zero line???
  ylab("Effect size (95% CRI)")+
  scale_x_discrete(limits = rev(levels(sppEfx$Species)))+
  coord_flip()+
  facet_wrap(~betaName, scales = "free_x", ncol = 4)+ 
  theme_bw()+
  theme(axis.title.x = element_text(vjust = 1.5, size = 18,colour = "black"),
        axis.text.x  = element_text(vjust = 0.5,size = 10,colour = "black", margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size=12, colour = "black", margin = unit(c(0.3,0.3,0.3,0.3), "cm")),
        axis.ticks.length = unit(-0.15, "cm"),
        strip.background = element_rect(fill = "gray95"),
        strip.text.x = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

## forest plot of shrub cover effects
ggplot(data = sppEfx[1:48,], aes(x = Species, y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4) +                     
  geom_hline(aes(yintercept = 0))+                ### keep zero line???
  ylab("Effect size (95% CRI)")+
  scale_x_discrete(limits = rev(levels(sppEfx$Species)))+
  coord_flip()+
  facet_wrap(~betaName, scales = "free_x", ncol = 4)+ 
  theme_bw()+
  theme(axis.title.x = element_text(vjust = 1.5, size = 18,colour = "black"),
        axis.text.x  = element_text(vjust = 0.5,size = 10,colour = "black", margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size=12, colour = "black", margin = unit(c(0.3,0.3,0.3,0.3), "cm")),
        axis.ticks.length = unit(-0.15, "cm"),
        strip.background = element_rect(fill = "gray95"),
        strip.text.x = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

#nesting substrate and diet fig
#nest effect
pd <- position_dodge(0.5) # move them .05 to the left and right
nest.summary$betaName <- factor(nest.summary$betaName,levels = c("Protected", "Pasture", "Homestead","Plantation","Shrub cover"))
diet.summary$betaName <- factor(diet.summary$betaName,levels = c("Protected", "Pasture", "Homestead","Plantation","Shrub cover"))

#nest substrates
nest.fig <- ggplot(data = nest.summary[nest.summary$betaName != "Shrub cover", ], aes(x = betaName, y = effect, group = nest)) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 0.3, position = pd) +
  geom_point(colour = "gray20", shape = 21, size = 2, position = pd, aes(fill = factor(nest))) + 
  scale_fill_manual(values = c("white", "yellow", "gray50", "blue"), name = "Nest substrates ")+
  ylab("Effect size (95% CRI)")+
  ylim(c(-15,15))+
  scale_x_discrete(limits = rev(levels(nest.summary$betaName)[1:4]))+
  coord_flip()+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size=9,colour = "black",margin=unit(c(0.2,0.2,0.2,0.2), "cm")),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.text=element_text(size=9),legend.title=element_text(face="bold",size=10))


diet.fig <- ggplot(data = diet.summary[diet.summary$betaName != "Shrub cover", ], aes(x = betaName, y = effect, group = diet)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 0.3, position = pd) +
  geom_point(colour = "gray20", shape = 21, size = 2, position = pd, aes(fill = factor(diet))) + 
  scale_fill_manual(values = c("white", "yellow", "gray50","blue","black"), name = "Diet")+
  ## geom_hline(aes(yintercept = 0))+
  ## annotation_custom(my_grob_b)+
  ylab("Effect size (95% CRI)")+
  ylim(c(-15,15))+
  scale_x_discrete(limits = rev(levels(nest.summary$betaName)[1:4]))+
  coord_flip()+
  theme_bw()+
  theme(axis.title.x = element_text(face="bold", vjust=0.3, size=11,colour = "black"),
        axis.text.x  = element_text(vjust = 0.5,size = 9,colour = "black",margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
        axis.title.y = element_blank(), 
        axis.text.y  = element_text(size = 9, colour = "black", margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
        axis.ticks.length = unit(-0.15, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9),legend.title = element_text(face = "bold", size = 10))

## Why is nectar effect missing?
nest.fig2 <- ggplot(data = nest.summary[nest.summary$betaName == "Shrub cover", ], aes(x = betaName, y = effect, group = nest)) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 0.3, position = pd) +
  geom_point(colour = "gray20", shape = 21, size = 2, position = pd, aes(fill = factor(nest))) + 
  scale_fill_manual(values = c("white", "yellow", "gray50", "blue"), name = "Nest substrates ")+
  geom_hline(aes(yintercept = 0))+
  ylab("Effect size (95% CRI)")+
  ylim(c(-1,1))+
  scale_x_discrete(limits = rev(levels(nest.summary$betaName)[5]))+
  coord_flip()+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size=9,colour = "black",margin=unit(c(0.2,0.2,0.2,0.2), "cm")),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.text=element_text(size=9),legend.title=element_text(face="bold",size=10))

## Why is nectar effect missing?
diet.fig2 <- ggplot(data = diet.summary[diet.summary$betaName == "Shrub cover", ], aes(x = betaName, y = effect, group = diet)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 0.3, position = pd) +
  geom_point(colour = "gray20", shape = 21, size = 2, position = pd, aes(fill = factor(diet))) + 
  scale_fill_manual(values = c("white", "yellow", "gray50","blue","black"), name = "Diet")+
  geom_hline(aes(yintercept = 0))+
  ## annotation_custom(my_grob_b)+
  ylab("Effect size (95% CRI)")+
  ylim(c(-1,1))+
  scale_x_discrete(limits = rev(levels(nest.summary$betaName)[5]))+
  coord_flip()+
  theme_bw()+
  theme(axis.title.x = element_text(face="bold", vjust=0.3, size=11,colour = "black"),
        axis.text.x  = element_text(vjust = 0.5,size = 9,colour = "black",margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
        axis.title.y = element_blank(), 
        axis.text.y  = element_text(size = 9, colour = "black", margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
        axis.ticks.length = unit(-0.15, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9),legend.title = element_text(face = "bold", size = 10))

#line up plots in 1 or 2 figures
grid.arrange(nest.fig, diet.fig, ncol = 1)
grid.arrange(nest.fig2, diet.fig2, ncol = 1)
grid.arrange(nest.fig, nest.fig2,diet.fig,  diet.fig2, ncol = 2)  ## legend repreats needlessly

##--------------------------------------------------------------------
##                Plotting effects by mass / "wing loading"
##--------------------------------------------------------------------
## filter sppEfx by type of beta
shrub_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Shrub cover")
prot_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Protected")
past_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Pasture")
home_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Homestead")
plant_mass <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Plantation")


## Attempt to plot species responses in order of log(mass) and overlay best-fit lines
mass.fig.shrub <- ggplot(data = shrub_mass, aes(x = log(Mass), y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4)+
  geom_abline(aes(intercept=0.07816,slope=0.01418), colour = "red", size = 2)+
  ylab("Effect size (95% CRI)")+
  theme_bw()+
  theme(axis.title.x = element_text(vjust = 1.5, size = 18,colour = "black"),
        axis.text.x  = element_text(vjust = 0.5,size = 10,colour = "black", margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size=12, colour = "black", margin = unit(c(0.3,0.3,0.3,0.3), "cm")),
        axis.ticks.length = unit(-0.15, "cm"),
        strip.background = element_rect(fill = "gray95"),
        strip.text.x = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")


mass.fig.prot <- ggplot(data = prot_mass, aes(x = log(Mass), y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4)+
  geom_abline(aes(intercept= -3.1605, slope = 1.3869), colour = "red", size = 2)+
  ylab("Effect size (95% CRI)")+
  theme_bw()+
  theme(#axis.title.x = element_text(vjust = 1.5, size = 18,colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.y = element_text(vjust = 1.5, size = 18,colour = "black"),
        axis.text.y  = element_text(size=12, colour = "black", margin = unit(c(0.3,0.3,0.3,0.3), "cm")),
        axis.ticks.length = unit(-0.15, "cm"),
        strip.background = element_rect(fill = "gray95"),
        strip.text.x = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.position = "none")

mass.fig.past <- ggplot(data = past_mass, aes(x = log(Mass), y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4)+
  geom_abline(aes(intercept= 1.6228, slope = -0.9917), colour = "red", size = 2)+
  theme_bw()+
  theme(#axis.title.x = element_text(vjust = 1.5, size = 18,colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.length = unit(-0.15, "cm"),
    strip.background = element_rect(fill = "gray95"),
    strip.text.x = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    legend.position = "none")


mass.fig.home <- ggplot(data = home_mass, aes(x = log(Mass), y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4)+
  geom_abline(aes(intercept= -2.6681, slope = 0.4170), colour = "red", size = 2)+
  ylab("Effect size (95% CRI)")+
  theme_bw()+
  theme(axis.title.x = element_text(vjust = 1.5, size = 18,colour = "black"),
    axis.title.y = element_text(vjust = 1.5, size = 18,colour = "black"),
    axis.text.y  = element_text(size=12, colour = "black", margin = unit(c(0.3,0.3,0.3,0.3), "cm")),
    axis.ticks.length = unit(-0.15, "cm"),
    strip.background = element_rect(fill = "gray95"),
    strip.text.x = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    legend.position = "none")


mass.fig.plant <- ggplot(data = plant_mass, aes(x = log(Mass), y = effect)) + 
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0) +
  geom_point(colour="gray20", shape=21, size = 4)+
  geom_abline(aes(intercept= -2.5434, slope = -0.8382), colour = "red", size = 2)+
  theme_bw()+
  theme(axis.title.x = element_text(vjust = 1.5, size = 18,colour = "black"),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.length = unit(-0.15, "cm"),
    strip.background = element_rect(fill = "gray95"),
    strip.text.x = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    legend.position = "none")

grid.arrange(mass.fig.shrub, ncol = 1)
grid.arrange(mass.fig.prot,mass.fig.past,mass.fig.home,mass.fig.plant, ncol = 2)



## Plotting code to delete
## Dirty plots of continuous effects, body size and "wing loading"
##-----------------------------------------------------------------------------------
## view mass/wing chord plotted against mean betas; only the shrub effect is a priori
shrubEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Shrub cover")
protEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Protected")
pastEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Pasture")
homeEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Homestead")
sugarEfxAndTraits <- merge(sppEfx, TraitData) %>% 
  filter(betaName == "Plantation")

## full posterior samples
par(mfrow = c(2,3))
summary(lm(Outputs_long$beta1~Outputs_long$pseudo_loading))
plot(Outputs_long$pseudo_loading,Outputs_long$beta1 )
abline(m1<-lm(Outputs_long$beta1~Outputs_long$pseudo_loading))

summary(lm(Outputs_long$`beta.l[1]`~Outputs_long$pseudo_loading))
plot(Outputs_long$pseudo_loading,Outputs_long$`beta.l[1]` )
abline(m1<-lm(Outputs_long$`beta.l[1]`~Outputs_long$pseudo_loading))

summary(lm(Outputs_long$`beta.l[2]`~Outputs_long$pseudo_loading))
plot(Outputs_long$pseudo_loading,Outputs_long$`beta.l[2]` )
abline(m1<-lm(Outputs_long$`beta.l[2]`~Outputs_long$pseudo_loading))

summary(lm(Outputs_long$`beta.l[3]`~Outputs_long$pseudo_loading))
plot(Outputs_long$pseudo_loading,Outputs_long$`beta.l[3]`)
abline(m1<-lm(Outputs_long$`beta.l[3]`~Outputs_long$pseudo_loading))

summary(lm(Outputs_long$`beta.l[4]`~Outputs_long$pseudo_loading))
plot(Outputs_long$pseudo_loading,Outputs_long$`beta.l[4]`)
abline(m1<-lm(Outputs_long$`beta.l[4]`~Outputs_long$pseudo_loading))

##----------------------------- MEANS ONLY
plot(protEfxAndTraits$pseudo_loading, protEfxAndTraits$effect)  
summary(lm(protEfxAndTraits$effect~protEfxAndTraits$pseudo_loading))

plot(pastEfxAndTraits$pseudo_loading, pastEfxAndTraits$effect)  
summary(lm(pastEfxAndTraits$effect~pastEfxAndTraits$pseudo_loading))

plot(homeEfxAndTraits$pseudo_loading, homeEfxAndTraits$effect)  
summary(lm(homeEfxAndTraits$effect~homeEfxAndTraits$pseudo_loading))

plot(sugarEfxAndTraits$pseudo_loading, sugarEfxAndTraits$effect)  
summary(lm(sugarEfxAndTraits$effect~sugarEfxAndTraits$pseudo_loading))

plot(shrubEfxAndTraits$pseudo_loading, shrubEfxAndTraits$effect)  
summary(lm(shrubEfxAndTraits$effect~shrubEfxAndTraits$pseudo_loading))
##-----------------------------------------------------------------------------------
## view mass plotted against betas
plot(log(shrubEfxAndTraits$Mass), shrubEfxAndTraits$effect)  
summary(lm(shrubEfxAndTraits$effect~log(shrubEfxAndTraits$Mass)))

plot(log(protEfxAndTraits$Mass), protEfxAndTraits$effect)        
summary(lm(protEfxAndTraits$effect~log(protEfxAndTraits$Mass)))

plot(log(pastEfxAndTraits$Mass), pastEfxAndTraits$effect)  
summary(lm(pastEfxAndTraits$effect~log(pastEfxAndTraits$Mass)))

plot(log(homeEfxAndTraits$Mass), homeEfxAndTraits$effect)  
summary(lm(homeEfxAndTraits$effect~log(homeEfxAndTraits$Mass)))

plot(log(sugarEfxAndTraits$Mass), sugarEfxAndTraits$effect)  
summary(lm(sugarEfxAndTraits$effect~log(sugarEfxAndTraits$Mass)))

