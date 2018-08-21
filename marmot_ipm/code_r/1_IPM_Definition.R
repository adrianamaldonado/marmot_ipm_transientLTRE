# !!! FIRST set your working directory to be wherever 'code_r' etc live !!!

# Load the imputed (for size) marmot dataset 
all.data <- read.csv(file.path(".", "data_raw", "all.data.csv"))
# read in the climate data -> the 'environment set'
load(file=file.path(".", "data_raw", "clim_data.rda"))

#============================================================================================
# (I) Define all global variables, parameters and demographic functions for computing the 
# kernel and read data here
#============================================================================================

# Define the range of years for which we have complete data 
transYrRange <- c(1976, 2012)

# Define body mass range for pups and everyone else in Jun and Aug
a0rangeJ <- c(07.5, 14.0)
a1rangeJ <- c(10.5, 16.5)
a0rangeA <- c(08.0, 14.0)
a1rangeA <- c(12.0, 18.0)

# Define the environment ste using range of years
eStateSet <- subset(rawClim, year >= transYrRange[1] & year <= transYrRange[2])
# clean up the names so they are a bit nicer
names(eStateSet)[names(eStateSet) == "mean_spT"]      <- "SprT"
names(eStateSet)[names(eStateSet) == "pwtemp"]        <- "WinT"  
names(eStateSet)[names(eStateSet) == "bareground"]    <- "BrGd"  
names(eStateSet)[names(rawClim) == "max_snowpack_sp"] <- "SwPk"
# use year for the row names and drop any variables that aren't part of the environment
rownames(eStateSet) <- eStateSet$year
eStateSet <- subset(eStateSet, select = -year)
eStateSet <- as.matrix(eStateSet)

# read in the fitted parameters and store them in a list object
read_derived <- function(f_name) {
  f_path <- file.path(".", "data_derived", f_name)
  read.csv(f_path, row.names=1, check.names=FALSE)
}
mParSet <- list()
mParSet[["GrowAJ"]] <- read_derived("1a.btCoefInt_ModGrowAJ.csv")
mParSet[["GrowJA"]] <- read_derived("1b.btCoefInt_ModGrowJA.csv")
mParSet[["RecrSz"]] <- read_derived( "2.btCoefInt_ModRecrSz.csv")
mParSet[["Surv"  ]] <- read_derived( "3.btCoefInt_ModSurv.csv"  )
mParSet[["Repr"  ]] <- read_derived( "4.btCoefInt_ModRepr.csv"  )
mParSet[["Recr"  ]] <- read_derived( "5.btCoefInt_ModRecr.csv"  )

# clean up the names so they are a bit nicer
lapply(mParSet, names)
mParSet <- lapply(mParSet, function(df) {
  names(df)[names(df) == "(Intercept)"]             <- "I"
  names(df)[names(df) == "mean_spT"]                <- "SprT"
  names(df)[names(df) == "pwtemp"]                  <- "WinT"  
  names(df)[names(df) == "bareground"]              <- "BrGd"  
  names(df)[names(df) == "max_snowpack_sp"]         <- "SwPk"
  names(df)[names(df) == "Af21"]                    <- "Idiff"
  names(df)[names(df) %in% c("Af21:zz","Af21:zz1")] <- "zdiff"
  names(df)[names(df) %in% c("zz","zz1","mom.zz")]  <- "z"
  return(df)
})
lapply(mParSet, names)

# clean up the dataframes
mParSet <- lapply(mParSet, function(df) {
  yr <- as.numeric(row.names(df))
  subset(df, yr >= transYrRange[1] & yr <= transYrRange[2])
})

# combine all the parameters
mParSet <- as.matrix(do.call("cbind", mParSet))
mParSetStore <- mParSet

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 1. Build the demographic functions that define the kernel components
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Functions to compute kernel components from the fitted models

## Ontogenic growth from Aug to Jun
# mean
mu_g_AJ_z1z <- function(z, A, eState, mPar) {
  eNames <- paste("GrowAJ", names(eState), sep=".")
  # linear predictor of z1
  z1mu <- mPar["GrowAJ.I"] + mPar["GrowAJ.z"] * z + mPar[eNames] %*% eState # multiply each environment state by the associated coef,then sum them up to get their contribution to the linear predictor
  if (A > 0) { # non pups individuals differ in their response (intercept)
    z1mu <- z1mu + mPar["GrowAJ.Idiff"] + mPar["GrowAJ.zdiff"] * z # adjust mean size of > 0
  }
  return(z1mu)
}
# mean and sd
g_AJ_z1z <- function(z1, z, A, eState, mPar)
{
  z1mu <- mu_g_AJ_z1z(z, A, eState, mPar)
  z1sd <- mPar["GrowAJ.sd"]
  return(dnorm(z1, z1mu, z1sd))
}

## Ontogenic growth from Jun to Aug
# mean
mu_g_JA_z1z <- function(z, A, eState, mPar) {
  eNames <- paste("GrowJA", names(eState), sep=".")
  # linear predictor of z1
  z1mu <- mPar["GrowJA.I"] + mPar["GrowJA.z"] * z + mPar[eNames] %*% eState
  if (A > 0) { # non pups individuals differ in their response (intercept)
    z1mu <- z1mu + mPar["GrowJA.Idiff"] + mPar["GrowJA.zdiff"] * z # adjust mean size of > 0
  }
  return(z1mu)
}
# mean and sd
g_JA_z1z <- function(z1, z, A, eState, mPar)
{
  z1mu <- mu_g_JA_z1z(z, A, eState, mPar)
  z1sd <- mPar["GrowJA.sd"]
  return(dnorm(z1, z1mu, z1sd))  
}

## Recruit size (mean and sd)
c_z1z <- function(z1, z, A, eState, mPar)
{
  eNames <- paste("RecrSz", names(eState), sep=".")
  z1mean <- mPar["RecrSz.I"] + mPar["RecrSz.z"] * z + mPar[eNames] %*% eState # mean size of recruits in Aug 
  z1sd <- mPar["RecrSz.sd"]        # sd about the mean 
  return(dnorm(z1, z1mean, z1sd))  # pdf for offspring size z1
}

## Survival (binary indicator)
s_z <- function(z, A, eState, mPar)
{
  eNames <- paste("Surv", names(eState), sep=".")
  nu <- mPar["Surv.I"] + mPar["Surv.z"] * z + mPar[eNames] %*% eState #  linear predictor
  return(1 / (1 + exp(-nu)))    # inv-logistic trans
}

## Reproduction (binary indicator)
pb_z <- function(z, A, eState, mPar)
{
  if (A > 0) { # The function is defined for individuals 1 yr or older
    eNames <- paste("Repr", names(eState), sep=".")
    nu <- mPar["Repr.I"] + mPar["Repr.z"] * z + mPar[eNames] %*% eState # linear predictor
    return(1 / (1 + exp(-nu)))  # inv-logistic trans
  } 
  return(rep(0, length(z)))     # no repro for newborn individuals (pups)
} 

## Recruitment (Number of female offspring)
nr_z <- function(z, A, eState, mPar)
{
  eNames <- paste("Recr", names(eState), sep=".")
  nu <- mPar["Recr.I"] + mPar["Recr.z"] * z + mPar[eNames] %*% eState # linear predictor
  return(exp(nu)) # inv-log trans
}


#=====================================================================================================
# (II) Build the kernel subcomponents and run dignostic plots
#=====================================================================================================

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 1. Define the reproduction kernel, F
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# conditional on survival (female has to survive first to be able to reproduce) and on the mom's mass in Aug previous year
F_z1z <- function (z1, z, A, eState, mPar) {
  return( s_z(z, A, eState, mPar) * pb_z(z, A, eState, mPar) * nr_z(z, A, eState, mPar) * c_z1z(z1, z, A, eState, mPar) )
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2. Define the survival-growth kernel, P
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## No survival-growth kernel function because of the seasonal growth model. It will be specified later
#  Step (1) winter survival-growth from Aug to Jun
#  Step (2) summer growth from Jun to Aug


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 3. Diagnostic plots
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## climatic means overall the years of study (mean state environment)
eStateMean <- apply(eStateSet, 2, mean)

## Growth plots
# 1. Aug -> Jun growth: age>0 individuals, mean environment
zvals <- seq(a1rangeA[1], a1rangeA[2], length=100) # vector that holds the potential size range for non-pup individuals
mu_z_years <- apply(mParSetStore, 1, # apply fun to the matrix that hold the parameters. Done by rows (1)
                    function(mPar) mu_g_AJ_z1z(zvals, A=1, eState=eStateMean, mPar=mPar)) # the linear predictor for non-pups using the range if zvals
# Plot the predicted sizes for each year (mu_z_years)
plot(0, 0, type="n", xlim=a1rangeA, ylim=c(11.5,18), xlab="z(Aug, t)", ylab="E[z(Jun, t+1)]")
apply(mu_z_years, 2, 
      function(mu_z) lines(zvals, mu_z, col=rgb(0,0,1,alpha=0.5)))
abline(0,1,lwd=2,lty=3)
lines(zvals, apply(mu_z_years, 1, mean), lwd=2) # mean size overall years
rug(subset(all.data, imputeNum==1 & AA>0)$zz)


# Jun -> Aug growth: age>0 individuals, mean environment
zvals <- seq(a1rangeJ[1], a1rangeJ[2], length=100) 
mu_z_years <- apply(mParSetStore, 1, 
                    function(mPar) mu_g_JA_z1z(zvals, A=1, eState=eStateMean, mPar=mPar))
# Plot the predicted sizes for each year (mu_z_years)
plot(0, 0, type="n", xlim=a1rangeJ, ylim=c(13, 17), xlab="z(t)", ylab="E[z(t+1)]")
apply(mu_z_years, 2, 
      function(mu_z) lines(zvals, mu_z, col=rgb(0,0,1,alpha=0.5)))
abline(0,1,lwd=2,lty=2)
lines(zvals, apply(mu_z_years, 1, mean), lwd=2)
rug(subset(all.data, imputeNum==1 & AA>0)$zz1)
rug(subset(all.data, imputeNum==1 & AA>0)$zz, side=2)

# survival: age>0 individuals, mean environment
zvals <- seq(a0rangeA[1], a1rangeA[2], length=100)
s_z_years <- apply(mParSetStore, 1, 
                   function(mPar) s_z(zvals, A=1, eState=eStateMean, mPar=mPar))
# Plot the predicted survival for each mass year
plot(0, 0, type="n", 
     xlim=c(a0rangeA[1], a1rangeA[2]), ylim=c(0, 1), xlab="z(t)", ylab="p(Survive)")
apply(s_z_years, 2, 
      function(s_z) lines(zvals, s_z, col=rgb(0, 0, 1, alpha=0.5)))
lines(zvals, apply(s_z_years, 1, mean), lwd=2)

# reproduction: age>0 individuals, mean environment
zvals <- seq(a1rangeA[1], a1rangeA[2], length=100)
pb_z_years <- apply(mParSetStore, 1, 
                    function(mPar) pb_z(zvals, A=1, eState=eStateMean, mPar=mPar))
# Plot the predicted reproduction for each mass year
plot(0, 0, type="n",
     xlim=c(a1rangeA[1], a1rangeA[2]), ylim=c(0, 1), xlab="z(t)", ylab="p(Reproduce)")
apply(pb_z_years, 2, 
      function(pb_z) lines(zvals, pb_z, col=rgb(0, 0, 1, alpha=0.5)))
lines(zvals, apply(pb_z_years, 1, mean), lwd=2, col="black")


#=====================================================================================================
# (III) Build the discretized year-specific Kernel, K, and the helper functions to iterate the model   
#=====================================================================================================

## Calculate the mesh points, mesh width and store with upper/lower bounds (integration parameters)
mk_intpar <- function(nm, l, u) {
  h <- (u - l) / nm                  # heigh of the mesh: max size - min size, divided by mesh size
  mesh  <-  l + ((1:nm) - 1/2) * h   # midpoint of the mesh
  return( list(mesh = mesh, h = h, nm = nm) )
}

## Create the year-especific parameter-environment sequence 
#  by random sampling with replacement the set of parameters (intercept and slopes) and climate variables
mk_env <- function(mParSet, eStateSet, nYrs) {
  im <- sample(nrow(mParSet),   nYrs, replace = TRUE)
  ie <- sample(nrow(eStateSet), nYrs, replace = TRUE)
  list(mParSet = mParSet[im,], eStateSet = eStateSet[ie,])
}


## Discretize kernel, contains F and P (year-specific kernel)
mk_IPM_iter_matrix <- function(eState, iPar, mPar) {
  F <- P <- list() # empty list to hold F and P
  with(iPar, {    
    # compute iteration matrix
    for (A in 0:1) {  # create P and F materices looping over age categories -- use parameters defined in iPar
      
      # F, reproduction (defined for non-pups)  -- outer computes the iteration matrix
      F[[A+1]] <<- outer(mesh, mesh, F_z1z, A = A, eState = eState, mPar = mPar) * h
      
      # P, survival-growth.
      # P1, for aug-jun growth
      g_AJ_iterM <- outer(mesh, mesh, 
                          function (z1, z, A, eState, mPar) {
                            s_z(z, A = A, eState = eState, mPar = mPar) * g_AJ_z1z(z1, z, A = A, eState = eState, mPar = mPar)
                          }, A = A, eState = eState, mPar = mPar) * h
      
      # P2, for jun-aug growth
      g_JA_iterM <- outer(mesh, mesh, g_JA_z1z, A = A, eState = eState, mPar = mPar) * h   
      # P is defined by growth from aug-jun and jun-aug
      P[[A+1]] <<- g_JA_iterM %*% g_AJ_iterM      
    }
  })
  return(list(F = F, P = P)) # Returns the kernel F and P as two different elements
}


#============================================================================================#
# (IV) Model iteration functions  
#============================================================================================#n 

# Iterate normal to get stable distribution size w(x) 
sim_IPM <- function(mParSet, eStateSet, iPar, nt0, normalize=TRUE) {
  # sanity check
  if (nrow(mParSet) != nrow(eStateSet)) # check matrix containing the coef for all demog fun and mtx containing clim 
    stop("Size of model parameter set and environment set must match")
  # 
  nSteps <- nrow(mParSet) # number of years in the model
  #
  if (missing(nt0)) { # if initial size distribution for the projection is not defined, use the following:
    nm <- iPar$nm # mesh size
    nt <- list(n0 = rep(1, nm)/nm, n1 = rep(1, nm) / (2*nm)) # initial size dist for pups and non-pups
  } else nt <- nt0 # use nt0 if it is defined. Should contain size dist for pups and non-pups
  #
  lamt <- numeric(nSteps)
  ntSim <- nt1 <- list()  # list to hold projected size distribution
  # loop to go through the # of sampled years to project the population
  ntSim[[1]] <- nt
  for (tt in seq_len(nSteps)) {
    # matrix iteration. M contains iterated values of F and P 
    M <- mk_IPM_iter_matrix(eStateSet[tt,], iPar, mParSet[tt,]) # use the parameters by row... specified by year
    # project time t + 1 and store projected F and P. Project for each age category pups and non-pups
    nt1[[1]] <- (M$F[[1]] %*% nt[[1]] + M$F[[2]] %*% nt[[2]])[,,drop=TRUE] # F[[1]] should be 0 b/c corresponds to pups
    nt1[[2]] <- (M$P[[1]] %*% nt[[1]] + M$P[[2]] %*% nt[[2]])[,,drop=TRUE]
    lamt[tt] <- sum(unlist(nt1))
    if (normalize) {
      nt1[[1]] <- nt1[[1]] / lamt[tt]      
      nt1[[2]] <- nt1[[2]] / lamt[tt]
    }
    ntSim[[tt+1]] <- nt <- nt1
  }
  if (normalize) return(list(ntSim=ntSim, lamt=lamt)) else return(ntSim)  
}

# ===============================================================================================
# (V) Plot the demographic functions
# ===============================================================================================

# Prepare the plotting window to include all the plots in a single frame
#png(file="FigGowFunc.png", width = 6, height = 6,units = "in", res = 300, pointsize=12) #save graph that goes in the paper
par(mfrow=c(2,2), mar=c(4,5,1,1)) #use the minimum margin area


## climatic means overall the years of study (mean state environment)
eStateMean <- apply(eStateSet, 2, mean)

## Growth plots
# 1. Aug -> Jun growth: age=0 individuals, mean environment
zvals <- seq(a0rangeA[1], a0rangeA[2], length=100) # vector that holds the potential size range for non-pup individuals
mu_z_years <- apply(mParSetStore, 1, # apply fun to the matrix that hold the parameters. Done by rows (1)
                    function(mPar) mu_g_AJ_z1z(zvals, A=0, eState=eStateMean, mPar=mPar)) # the linear predictor for non-pups using the range if zvals
# Plot the predicted sizes for each year (mu_z_years)
plot(0, 0, type="n", xlim=a0rangeA, ylim=c(7,14), axes=FALSE, tcl=0.5,xlab="", ylab="", cex.lab=1.2)#xlab="z(Aug, t)", ylab="E[z(Jun, t+1)]")
apply(mu_z_years, 2, 
      function(mu_z) lines(zvals, mu_z, col=rgb(128,128,128,maxColorValue=200))) #alpha=0.5,
abline(0,1,lwd=2,lty=3)
lines(zvals, apply(mu_z_years, 1, mean), lwd=2) # mean size overall years
rug(subset(all.data, imputeNum==1 & AA==0)$zz)
box(lwd="1.5", bty="o")
title(xlab="sqrt(mass) in August in year t", line= 2.2, cex.lab=1.2, mgp=c(2,1,0))
title(ylab="sqrt(nass) in June in year t + 1", line= 2.2, cex.lab=1.2)
# y axis
b <- seq(7,14, by= 1) # a0rangeJ[1]; a0rangeJ[2]
axis(side = 2, tck = -0.02,  at = b, labels= b, las =1, line = 0)
# x axes
b <- seq(8, 14, by= 1) # a0rangeA[1]; a0rangeA[2]
axis(side = 1, tck = -0.02, at = b,  labels= b) 
mtext("a.", side=3, line= -1.5, adj=0.1, font=2, cex=1.2)

# 2. Aug -> Jun growth: age>0 individuals, mean environment
zvals <- seq(a1rangeA[1], a1rangeA[2], length=100) # vector that holds the potential size range for non-pup individuals
mu_z_years <- apply(mParSetStore, 1, # apply fun to the matrix that hold the parameters. Done by rows (1)
                    function(mPar) mu_g_AJ_z1z(zvals, A=1, eState=eStateMean, mPar=mPar)) # the linear predictor for non-pups using the range if zvals
# Plot the predicted sizes for each year (mu_z_years)
plot(0, 0, type="n", xlim=a1rangeA, ylim=c(10,17), axes=FALSE, tcl=0.5,xlab="", ylab="", cex.lab=1.2)#xlab="z(Aug, t)", ylab="E[z(Jun, t+1)]")
apply(mu_z_years, 2, 
      function(mu_z) lines(zvals, mu_z, col=rgb(128,128,128,maxColorValue=200))) #alpha=0.5,
abline(0,1,lwd=2,lty=3)
lines(zvals, apply(mu_z_years, 1, mean), lwd=2) # mean size overall years
rug(subset(all.data, imputeNum==1 & AA>0)$zz)
box(lwd="1.5", bty="o")
title(xlab="sqrt(mass) in August in year t", line= 2.2, cex.lab=1.2, mgp=c(2,1,0))
title(ylab="sqrt(nass) in June in year t + 1", line= 2.2, cex.lab=1.3)
# y axis
b <- seq(10, 17, by= 1) # a1rangeJ[1]; a1rangeJ[2] 
axis(side = 2, tck = -0.02,  at = b, labels= b, las =1, line = 0)
# x axes
b <- seq(12, 18, by= 1) # a1rangeA[1]; a1rangeA[2]
axis(side = 1, tck = -0.02, at = b,  labels= b) 
mtext("b.", side=3, line= -1.5, adj=0.1, font=2, cex=1.2)

# 3.Jun -> Aug growth: age=0 individuals, mean environment
zvals <- seq(a0rangeJ[1], a0rangeJ[2], length=100) # vector that holds the potential size range for non-pup individuals
mu_z_years <- apply(mParSetStore, 1, # apply fun to the matrix that hold the parameters. Done by rows (1)
                    function(mPar) mu_g_JA_z1z(zvals, A=0, eState=eStateMean, mPar=mPar)) # the linear predictor for non-pups using the range if zvals
# Plot the predicted sizes for each year (mu_z_years)
plot(0, 0, type="n", xlim=a0rangeJ, ylim=c(13, 16), axes=FALSE, tcl=0.5,xlab="", ylab="", cex.lab=1.2)#xlab="z(Aug, t)", ylab="E[z(Jun, t+1)]")
apply(mu_z_years, 2, 
      function(mu_z) lines(zvals, mu_z, col=rgb(128,128,128,maxColorValue=200))) #alpha=0.5,
abline(0,1,lwd=2,lty=3)
lines(zvals, apply(mu_z_years, 1, mean), lwd=2) # mean size overall years
rug(subset(all.data, imputeNum==1 & AA==0)$zz1)
box(lwd="1.5", bty="o")
title(xlab="sqrt(mass) in June in year t + 1", line= 2.2, cex.lab=1.2, mgp=c(2,1,0))
title(ylab="sqrt(nass) in August in year t + 1", line= 2.2, cex.lab=1.2)
# y axis
b <- seq(13, 16, by= 1) # range(mu_z_years)
axis(side = 2, tck = -0.02,  at = b, labels= b, las =1, line = 0)
# x axes
b <- seq(7,14, by= 1) # a0rangeJ[1]; a0rangeJ[2]
axis(side = 1, tck = -0.02, at = b,  labels= b) 
mtext("c.", side=3, line= -1.5, adj=0.1, font=2, cex=1.2)

# 4.Jun -> Aug growth: age>0 individuals, mean environment
zvals <- seq(a1rangeJ[1], a1rangeJ[2], length=100) # vector that holds the potential size range for non-pup individuals
mu_z_years <- apply(mParSetStore, 1, 
                    function(mPar) mu_g_JA_z1z(zvals, A=1, eState=eStateMean, mPar=mPar))
# Plot the predicted sizes for each year (mu_z_years)
plot(0, 0, type="n", xlim=a1rangeJ, ylim=c(13,17), axes=FALSE, tcl=0.5,xlab="", ylab="", cex.lab=1.2)#xlab="z(Aug, t)", ylab="E[z(Jun, t+1)]")
apply(mu_z_years, 2, 
      function(mu_z) lines(zvals, mu_z, col=rgb(128,128,128,maxColorValue=200))) #alpha=0.5,
abline(0,1,lwd=2,lty=3)
lines(zvals, apply(mu_z_years, 1, mean), lwd=2) # mean size overall years
rug(subset(all.data, imputeNum==1 & AA>0)$zz1)
box(lwd="1.5", bty="o")
title(xlab="sqrt(mass) in June in year t + 1", line= 2.2, cex.lab=1.2, mgp=c(2,1,0))
title(ylab="sqrt(nass) in August in year t + 1", line= 2.2, cex.lab=1.2)
# y axis
b <- seq(13,17, by= 1) # range(mu_z_years)
axis(side = 2, tck = -0.02,  at = b, labels=b, las =1, line = 0)
# x axes
b <- seq(10, 17, by= 1) # a1rangeJ[1]; a1rangeJ[2]
axis(side = 1, tck = -0.02, at = b,  labels= b) 
mtext("d.", side=3, line= -1.5, adj=0.1, font=2, cex=1.2)

# close device
#dev.off()

## surv and repro plots
#Prepare the plotting window to include all the plots in a single frame
#png(file="FigSurvRepFunc.png", width = 6, height = 6,units = "in", res = 300, pointsize=12)#save graph that goes in the paper
par(mfrow=c(2,2), mar=c(4,5,1,1)) #use the minimum margin area

# 1. survival: age=0 individuals, mean environment
zvals <- seq(a0rangeA[1], a0rangeA[2], length=100)
s_z_years <- apply(mParSetStore, 1, 
                   function(mPar) s_z(zvals, A=0, eState=eStateMean, mPar=mPar))
plot(0, 0, type="n", 
     xlim=c(a0rangeA[1], a0rangeA[2]), ylim=c(0, 1), axes=FALSE, tcl=0.5,xlab="", ylab="", cex.lab=1.2)#xlab="z(Aug, t)", ylab="E[z(Jun, t+1)]")
apply(s_z_years, 2, 
      function(s_z) lines(zvals, s_z,  col=rgb(128,128,128,maxColorValue=200)))
lines(zvals, apply(s_z_years, 1, mean), lwd=2)
rug(subset(all.data, imputeNum==1 & AA==0)$zz)
box(lwd="1.5", bty="o")
title(xlab="sqrt(mass) in August in year t", line= 2.2, cex.lab=1.2, mgp=c(2,1,0))
title(ylab="Probability of survival", line= 2.3, cex.lab=1.2)
# y axis
b <- seq(0,1, by= 0.2) # range(s_z_years)
axis(side = 2, tck = -0.02,  at = b, labels=round(b, 1), las =1, line = 0)
# x axes
b <- seq(8, 14, by= 1) # a0rangeA[1]; a0rangeA[2]
axis(side = 1, tck = -0.02, at = b,  labels= b) 
mtext("a.", side=3, line= -1.5, adj=0.1, font=2, cex=1.2)


# 2. survival: age>0 individuals, mean environment
zvals <- seq(a1rangeA[1], a1rangeA[2], length=100)
s_z_years <- apply(mParSetStore, 1, 
                   function(mPar) s_z(zvals, A=1, eState=eStateMean, mPar=mPar))
plot(0, 0, type="n", 
     xlim=c(a1rangeA[1], a1rangeA[2]), ylim=c(0, 1), axes=FALSE, tcl=0.5,xlab="", ylab="", cex.lab=1.2)#xlab="z(Aug, t)", ylab="E[z(Jun, t+1)]")
apply(s_z_years, 2, 
      function(s_z) lines(zvals, s_z,  col=rgb(128,128,128,maxColorValue=200)))
lines(zvals, apply(s_z_years, 1, mean), lwd=2)
rug(subset(all.data, imputeNum==1 & AA == 1)$zz)
box(lwd="1.5", bty="o")
title(xlab="sqrt(mass) in August in year t", line= 2.2, cex.lab=1.2, mgp=c(2,1,0))
title(ylab="Probability of survival", line= 2.3, cex.lab=1.2)
# y axis
b <- seq(0,1, by= 0.2) # range(s_z_years)
axis(side = 2, tck = -0.02,  at = b, labels=round(b, 1), las =1, line = 0)
# x axes
b <- seq(12, 18, by= 1) # a1rangeA[1]; a1rangeA[2]
axis(side = 1, tck = -0.02, at = b,  labels= b) 
mtext("b.", side=3, line= -1.5, adj=0.1, font=2, cex=1.2)


# 3 reproduction: age>0 individuals, mean environment
zvals <- seq(a1rangeA[1], a1rangeA[2], length=100)
pb_z_years <- apply(mParSetStore, 1, 
                    function(mPar) pb_z(zvals, A=1, eState=eStateMean, mPar=mPar))
# Plot
plot(0, 0, type="n",
     xlim=c(a1rangeA[1], a1rangeA[2]), ylim=c(0, 1), axes=FALSE, tcl=0.5,xlab="", ylab="", cex.lab=1.2)
apply(pb_z_years, 2, 
      function(pb_z) lines(zvals, pb_z, col=rgb(128,128,128,maxColorValue=200)))
lines(zvals, apply(pb_z_years, 1, mean), lwd=2, col="black")
rug(subset(all.data, imputeNum==1 & AA == 1)$zz)
box(lwd="1.5", bty="o")
title(xlab="sqrt(mass) in August in year t", line= 2.2, cex.lab=1.2, mgp=c(2,1,0))
title(ylab="Probability of reproduction", line= 2.3, cex.lab=1.2)
# y axis
b <- seq(0,1, by= 0.2) # range(s_z_years)
axis(side = 2, tck = -0.02,  at = b, labels=round(b, 1), las =1, line = 0)
# x axes
b <- seq(12, 18, by= 1) # a1rangeA[1]; a1rangeA[2]
axis(side = 1, tck = -0.02, at = b,  labels= b) 
mtext("c.", side=3, line= -1.5, adj=0.1, font=2, cex=1.2)

# 3 recruitment: age>0 individuals, mean environment
zvals <- seq(a1rangeA[1], a1rangeA[2], length=100)
nr_z_years <- apply(mParSetStore, 1, 
                    function(mPar) nr_z(zvals, A=1, eState=eStateMean, mPar=mPar))
# Plot
plot(0, 0, type="n",
     xlim=c(a1rangeA[1], a1rangeA[2]), ylim=c(1, 3), axes=FALSE, tcl=0.5,xlab="", ylab="", cex.lab=1.2)#xlab="z(Aug, t)", ylab="E[z(Jun, t+1)]")
apply(nr_z_years, 2, 
      function(nr_z) lines(zvals, nr_z, col=rgb(128,128,128, maxColorValue=200)))
lines(zvals, apply(nr_z_years, 1, mean), lwd=2, col="black")
rug(subset(all.data, imputeNum==1 & AA == 1)$zz)
box(lwd="1.5", bty="o")
title(xlab="sqrt(mass) in August in year t", line= 2.2, cex.lab=1.2, mgp=c(2,1,0))
title(ylab="Recruitment", line= 2.3, cex.lab=1.2)
# y axis
b <- seq(1, 3, by= 0.5) # range(nr_z_years)
axis(side = 2, tck = -0.02,  at = b, labels=round(b, 1), las =1, line = 0)
# x axes
b <- seq(12, 18, by= 1) # a1rangeA[1]; a1rangeA[2]
axis(side = 1, tck = -0.02, at = b,  labels= b) 
mtext("d.", side=3, line= -1.5, adj=0.1, font=2, cex=1.2)

# close device
#dev.off()
par(mfrow=c(1,1))

# 4. recruitment size, for A >0
z1vals <- seq(a1rangeJ[1], a1rangeJ[2], length=100) 
zvals <- seq(a1rangeA[1], a1rangeA[2], length=100) 
mu_z_years <- apply(mParSetStore, 1, 
                    function(mPar) c_z1z(z1vals, zvals, A=1, eState=eStateMean, mPar=mPar))
# Plot the predicted sizes for each year (mu_z_years)
plot(0, 0, type="n", xlim=c(a1rangeJ[1], a1rangeA[2]), ylim=c(0, 1), axes=FALSE, tcl=0.5,xlab="", ylab="", cex.lab=1.2)#xlab="z(Aug, t)", ylab="E[z(Jun, t+1)]")
apply(mu_z_years, 2, 
      function(mu_z) lines(zvals, mu_z, col=rgb(128,128,128,maxColorValue=200))) 
abline(0,1,lwd=2,lty=3)
lines(zvals, apply(mu_z_years, 1, mean), lwd=2) # mean size overall years
rug(subset(all.data, imputeNum==1 & AA==1)$zz)
box(lwd="1.5", bty="o")
title(xlab="sqrt(mass) in August in year t", line= 2.2, cex.lab=1.2, mgp=c(2,1,0))
title(ylab="Recruit sqrt(mass) in year t + 1", line= 2.2, cex.lab=1.2)
# y axis
b <- seq(13, 16, by= 1) # range(mu_z_years)
axis(side = 2, tck = -0.02,  at = b, labels= b, las =1, line = 0)
# x axes
b <- seq(12, 18, by= 1) # a1rangeA[1]; a1rangeA[2]
axis(side = 1, tck = -0.02, at = b,  labels= b) 
mtext("e.", side=3, line= -1.5, adj=0.1, font=2, cex=1.2)
