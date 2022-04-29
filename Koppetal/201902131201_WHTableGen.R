#read in necessary packages
library(BDP2)
library(Hmisc)
library(clinfun)

# These are the fns that BDP2 uses for determining stopping rules.
# Note that they are slightly different with respect to v.boundary.

CritBoundaryE5<-
function (n, pE, critE, shape1 = 1, shape2 = 1) 
{
  v.boundary <- rep(n + 1, n)
  for (i in 1:n) {
    for (k in 0:i) {
      z <- pbeta(q = pE, shape1 = shape1 + k, shape2 = shape2 + 
                   i - k, lower.tail = FALSE)
      if (z >= critE) {
        v.boundary[i] <- k
        break
      }
      v.boundary[i] <- k
    }
    if (z < critE) 
      v.boundary[i] = i + 1
  }
  return(v.boundary)
}

CritBoundaryF5<-
function (n, pF, critF, shape1 = 1, shape2 = 1) 
{
  v.boundary <- rep(-1, n)
  for (i in 1:n) {
    for (k in 0:i) {
      z <- pbeta(q = pF, shape1 = shape1 + k, shape2 = shape2 + 
                   i - k, lower.tail = FALSE)
      if (z >= critF) 
        break
      v.boundary[i] <- k
    }
  }
  return(v.boundary)
}

# Parameterize the models from the analysis plan:

npts=100

# Beta priors:

shape1F=0.3
shape2F=0.7
shape1E=0.001 # 0.001 are pretty much uninformative priors
shape2E=0.001 # 0.001 are pretty much uninformative priors

results <- data.frame(shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
	              p0=0,p1=0,pF=0,pE=0,cF=0,cE=0,NEnrolledPts=1:1, EffStoppingBound = c(0,0), FutStoppingBound = c(0,0))

  for (p0 in seq(0.05, 0.30, by=0.05)) {
    for (p1 in seq(0.4, 0.80, by=0.2)) {
      for (cE in seq(0.95, 0.99, by=0.01)) {
      pF=p1
      pE=p0
      cF=0.01

      # These are 1 to n vectors with the stopping bounds using beta posteriors.
      # bF is the stopping boundary for futility, bE for efficacy
      bF <- CritBoundaryF5(n = npts, pF = pF, critF = cF, shape1 = shape1F, shape2 = shape2F)
      bE <- CritBoundaryE5(n = npts, pE = pE, critE = cE, shape1 = shape1E, shape2 = shape2E)

      newresults <- data.frame(shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
                               p0=p0,p1=p1,pF=pF,pE=pE,cF=cF,cE=cE,NEnrolledPts=1:npts, EffStoppingBound = bF, FutStoppingBound = bE)
      results <- rbind(results, newresults)
      }}}
write.csv(results, file="201902131201_WHTableGen_Out.csv", na = "")
	    
