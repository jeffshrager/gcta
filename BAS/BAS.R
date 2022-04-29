# Bayesian Adaptive Clinical Trial sim; adapted from the code provided
# around June 2013 by Chunyan Cai <Chunyan.Cai@uth.tmc.edu> from: Cai,
# et al. 2013 Adaptive Trials
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3867529/
# To run from R listener:
# source("BAS.R")
# run()

library(mvtnorm)
library(MASS)

set.seed(5)

# Basic params
step <- 10 # was 2000 # Number of simulated trials (Each round takes about 60 seconds on a 4G 2.66ghz i7 macbook pro.)
nmax <- 200 # sample size
nstart <- 32 # number of patients in the stage I burn-in period
Tf=10 # Follow up time
taccrual=30/12 # Accrual rate
  
# Priors params
g <- 10 # hyperparamter specified in the "Prior specification" section of CYJ2013
sigma <- sqrt(130)

# Early termination params
alpha <- 0.35
delta <- 10

# Sample size for numerical integration
n.sim <- 500 # was 10000 
# If this is too low you can end up with singular matrix (i.e., which is non-invertable), 
# This will report an errors as:
#    solve.default(t(xx) %*% xx + 1/g * diag(p)) : 
#    system is computationally singular: reciprocal condition number = 0
  
# True treatment effect for 16 combinations (scenario 1)
theta <- c(0, -25, -13, -12, -11, -20, -19, -18, -14, -15, -12, -10, -16, -17, -14, -21)

# Vector of treatments with highest efficay (i.e., lowest values in theta, above)
# Example: If there are 3 treatments (trt2, trt3, trt6) with highest efficacy probability, 
# then trt.eff <- c("2", "3", "6") (??? Why isn't this just computed from theta?)
trt.eff <- c("2")

# Main function, call by: 
run<-function(tracelevel=1, nsimdisplaylimit=100)
{
    setwd("/Users/jeffshrager/Desktop/cc/gcta/BAS")
    # These get matching globals
    NSimDisplayLimit <<- nsimdisplaylimit
    TraceLevel <<- tracelevel
	simu.res(step, g, sigma, alpha, delta, theta, nstart, trt.eff)
}

simu.res <- function(step, g, sigma, alpha, delta, theta, nstart, trt.eff)
{
  # ------ WARNING: VERY CONFUSING! ------ 
  # Read outcome file for initial patients/subjects, using +1 -1 to denote the design matrix.  
  # The file contains nstart rows, each being a combination of agents using +1/-1 in 5 columns (first is label), 
  # as: 1 1 1 1 1 // 2 -1 1 1 1 // 3 -1 -1 1 1 // and so on for each of the nstart patients, 
  # representing all combinations of agents, repeated as needed to make nstart subjects. 
  # This is used to build the model, X.
  # (??? Could this just be computed; is there 
  # any case where this model is otherwise than all combinations???)

  # This is at the same time clever and kludgy. It's using this external file that
  # is intented as the first nstart (x2 in this case, i.e., 32) to initialize the model
  # and also to provide the sourcing for each of the new patients added (up to nmax). 
  # So that when a new patient arrives, a row from these first n (16) gets pulled and 
  # replicated as a new row. So, where below see: "x <- rbind(x, A[m[nnn],])", the A came
  # (just below here) as the first 16 rows of X, so what the rbind is doing is adding a
  # new subject line to the model (X). 
  # ------ WARNING: VERY CONFUSING! ------ 

  data = matrix(scan("initial_diag.dat",skip=1),ncol=5,byrow=T)  
  X = data[,2:5] # Create the initial model (skipping label first col)
  
  # Looks like this overrides the call argument ???
  # At the moment these happen to be the same, fortunately.
  nstart <- dim(data)[1] 
  
  # Add intercept and interaction terms.
  # Intercept = 1, Curcumin=2, Minocycline=3, Modafinal=4, Buprion=5
  # set up interactions explicitly to avoid confusion with indices
  X = cbind(c(rep(1,nstart)),X)  # add intercept -- all 1
  # Cur*Min  in c6
  X = cbind(X,X[,2]*X[,3])
  # Cur*Mod in c7
  X = cbind(X,X[,2]*X[,4])
  # Cur*Bup in c8
  X = cbind(X,X[,2]*X[,5])
  # Min*Mod in c9
  X = cbind(X,X[,3]*X[,4])
  # Min*Bup in c10
  X = cbind(X,X[,3]*X[,5])
  # Mod*Bup in c11
  X = cbind(X,X[,4]*X[,5])
  # Cur*Min*Mod  in c12
  X = cbind(X,X[,2]*X[,3]*X[,4])
  # Cur*Min*Bup  in c13
  X = cbind(X,X[,2]*X[,3]*X[,5])
  # Cur*Mod*Bup  in c14
  X = cbind(X,X[,2]*X[,4]*X[,5])
  # Min*Mod*Bup  in c15
  X = cbind(X,X[,3]*X[,4]*X[,5])
  # Four way in c16
  X = cbind(X,X[,2]*X[,3]*X[,4]*X[,5])
  ..("Model (X) with intercept (1) and interaction terms",X,1)

  # Logging matrices. These are always STEP (i.e., number of runs) deep, by whatever wide.
  # Logs the posterior probability for each hypothesis,there are 16 models including null model
  mmpp <<- array(0, c(step, nmax, 16))
  # Logs the posterior probability pr(each treatment < placebo-delta)
  postsig <<- array(0, c(step, nmax, 16))
  # Logs the treatment assignment for each patient
  nnnn <<- matrix(0, step, nmax)   

  num.last <- vector("numeric", step)
  stopvalue <- vector("numeric", step)
  
  # Design matrix to indicate 16 combinations 
  A <- X[1:16, ]
 
  # True beta
  ..("True efficacies (theta)",theta,1)	
  beta <- solve(A)%*%theta # remember solve(m) is just inverse(m)
  ..("[partial result: solve.invert(A)]",solve(A),2)
  ..("True beta (beta <- solve.invert(A) %*% theta)",beta,1)

  # Display fixed parameters 
  ..("Prior hyperpameter (g)",g)
  ..("Total available patients (nmax)",nmax)
  ..("Followup time (months?) (Tf)",Tf)
  ..("Accural rate per month (I guess??? e.g, 30/12) (taccrual)",taccrual)
  ..("Samples for integration (n.sim)",n.sim)
  ..("Design matrix (A)",A)
  ..("True beta",beta)
  ..("sigma",sigma)
  ..("Model matrix (X in the main fn) (x)",X)
  ..("Number of burn-in subjects (nstart)",nstart)
  ..("delta",delta)
  ..("alpha",alpha)

  # Setup for timing
  print(proc.time())
  ptm <- proc.time()
  
  # Main adaptation loop
  for (ii in 1:step)
  {
    print(proc.time())
    print(proc.time()-ptm)
    ..("****************** Round",ii)
    res <- adaptdesign(nmax, Tf, taccrual, n.sim, A, beta, sigma, X, nstart, delta, alpha)

    # Store for logging and next time around
    mmpp[ii,,] <<- res$mpp
    ..("res$mpp",res$mpp,1)
    nnnn[ii,] <<- res$m
    ..("new nnn = res$m",res$m,1)
    postsig[ii,,] <<- res$effprob
    ..("res$effprob = new postsig",res$effprob,1)
    stopvalue[ii] <- res$stopflag
    ..("res$stopflag",res$stopflag,1)
    num.last[ii] <- res$nlast
    ..("res$nlast",res$nlast,1)
    
    # Backup occassionally
    if (ii %% 10==0) 
      {cat("(autosave!)")
       save(list=ls(),file = 'mvntest1.Rdata')
       }

  } # Close of main adaptaion loop
  
  save(list=ls(),file = 'mvntest1.Rdata')
  print(proc.time())
  print(proc.time()-ptm)

  # Summarize results
  
  # Select the treatment with the highest probability (that wasn't futile, i.e., stopped)
  trtdec <- table(apply(mmpp[stopvalue==0,200,], 1, which.max))
  ..("trtdec <- table(apply(mmpp[stopvalue==0,200,], 1, which.max))",trtdec,1)
    trtdec1 <- trtdec/step*100
  ..("trtdec1 <- trtdec/step*100",trtdec1,1)
  trtdec2 <- round(trtdec1,1)
  cat("selection % for each combination (trtdec2): ",trtdec2,"\n")
  cat("percentage of patients assigned to each combination: ",round((table(nnnn)/sum(num.last)*100),1),"\n")	
  cat("selection % for the efficacious combinations: ",trtdec2[trt.eff],"\n")
  cat("percentage of patients assigned to the efficacious combination: ",round((table(nnnn)/sum(num.last)*100)[trt.eff],1),"\n")
  cat("(See mmpp, nnnn, postsig for result dynamics.)")
}

adaptdesign <- function(nmax, Tf, taccrual, n.sim, A, beta, sigma, x, nstart, delta, alpha)
{
  ..(">>>>>>>> adaptdesign",,2)

  stopflag <- 0
  nnn <- nstart # This values is used (and destroyed) below, thus this odd replication
  
  # Outcome value for first nstart patients
  y = x%*%beta+rnorm(nnn, 0, sigma)
  ..("Outcome value for patients to date (y = x%*%beta+rnorm(nnn, 0, sigma))",y,2)

  # Calcuate the normalizing constant for each prior density function
  p <- 16
  betaprime <- mvrnorm(n.sim, mu=rep(0,p), Sigma=diag(p))
  ..("betaprime <- mvrnorm(n.sim, mu=rep(0,p), Sigma=diag( p ))",betaprime)
  c <- vector("numeric", 16)
  r <- vector("numeric", 16)
  sample <- betaprime

  # to calculate the theta for each model
  temp <- sample%*%t(A)
  ..("temp <- sample%*%t(A)",temp,2)
  # to store which theta is minimum
  rankc <- factor(apply(temp, 1, which.min), level=c(1:16))
  ..("rankc <- factor(apply(temp, 1, which.min), level=c(1:16))",rankc,2)
  c <- as.vector(table(rankc))
  ..("c <- as.vector(table(rankc))",c,2)
  c <- c/n.sim
  ..("normlized fraction each treatment combination won, of n.sim (c <- c/n.sim)",c)
  
  # to store the posterior probability after each patient for each model
  mpp <- matrix(0, nmax, 16) # posterior
  effprob <- matrix(0, nmax, 16)
  
  time.start=NULL
  
  # Skip past accrual of first nstart (32) subjects, incrementing the time by taccrual for each.
  for (i in 1:nstart)
  {
    time.start=c(time.start, 0)
    time.start=time.start+taccrual
  }
    
  # This pattern: time setup, then ncomp,xcomp,ycomp setup, then calling prob(..)
  # and then extracting mpp (new posterior) and effprob, is repeated three times. 
  # Here is the first call, then there is an if/then/else the then arm of which does 
  # this for all "inner" subjects, and the else arm does it for the last subject.
  # The ind.comp part followed by setting ncomp,xcomp,ycomp is choosing the subjects
  # who need to be observed. FFF Refactor!
  
  # The first time through this always ends up with 29 subjects, so don't be surprised 
  # when you end up with xcomp being only 29 rows long. It's because only 29 subjects are 
  # beyond the 10 month (default) observation time.

  # Pattern body (first subject):
  ind.comp <- time.start >= Tf	
  ncomp <- sum(ind.comp)
  xcomp <- x[1:ncomp,]
  ycomp <- y[1:ncomp]
  probres <- prob(n.sim, c, xcomp, ycomp, nnn, delta, A)
  mpp[nnn,] <- probres$mp
  effprob[nnn,] <- probres$effprob
  
  # to store the random uniform number from 0 and 1
  test <- vector('numeric', nmax)
  
  # m stores the treatment indicator assigned to each patient
  m <- vector('numeric', nmax)
  # initiate to store the treatment indicator for first 32 patients
  m[1:nstart] <- rep(1:16,nstart/16)

  if (nnn>=nstart)  # Check futility only after the burn in period.
                    # (Why would nnn ever be < nstart at this point?!)
    {stopflag <- futility.rule(effprob[nnn,], alpha)}

  # Central loop for nmax inner subjects (actually, nnn will start at nstart+1, usually 33).
  # Go subject by subject, incrementing the time by the taccrual for each subject accrued.
  while (stopflag==0 & nnn<nmax)
  {
    nnn <- nnn+1  

    ..("++++++++++++++++++++++ Determining treatment (by posterior mpp) for subject (nnn)",nnn,1)
    test[nnn] <- runif(1)
    ..("From (mpp[nnn-1,])",mpp[nnn-1,],1)
    # Select the treatment in accord with each one's prob. per the unirand just chosen (via the "sumup" trick).
    m[nnn] <- min(which(cumsum(mpp[nnn-1,]) >= test[nnn]))
    ..("selected treatment (m[nnn] <- min(which(cumsum(mpp[nnn-1,]) >= test[nnn])))",m[nnn],1)
    # Append the selected treatment row from the A design matrix to the X model matrix
    x <- rbind(x, A[m[nnn],])
    ..("x bound with new treatment rows from design (x <- rbind(x, A[m[nnn],]))",x,2)
    # Generate the new outcome by crossing the randomly chosen A row with true outcomes (beta) and adding noise
    ynew <- A[m[nnn],]%*%beta+rnorm(1, 0, sigma)
    ..("new outcome (ynew <- A[m[nnn],]%*%beta+rnorm(1, 0, sigma))",ynew,2)
    y <- c(y, ynew)
    ..("new y vector (y <- c(y, ynew))",y,2)
    ..(paste("patient: ",nnn,", treated with: ",m[nnn],", outcome: ",ynew,"\n"))
    # Setup for the next (or last) subject.
    if (nnn<nmax)
    {
      # Here's the pattern again: time setup, then ncomp,xcomp,ycomp setup, then calling prob(..)
      # this time for all "inner" subjects.
      time.start <- c(time.start, 0)
      time.start <- time.start+taccrual
      ind.comp <- time.start>=Tf
      ncomp <- sum(ind.comp)
      xcomp <- x[1:ncomp,]
      ycomp <- y[1:ncomp]
      probres <- prob(n.sim, c, xcomp, ycomp, nnn, delta, A)
    }
    else {
      # And here's the patterns for the last subject. (This only differs from the above in the time.start calc)
      time.start <- c(time.start, 0)
      time.start <- time.start+Tf
      ind.comp <- time.start>=Tf
      ncomp <- sum(ind.comp)
      xcomp <- x[1:ncomp,]
      ycomp <- y[1:ncomp]
      probres <- prob(n.sim, c, xcomp, ycomp, nnn, delta, A)
    }

    # Now rebind the results from one of the above prob() calls for the next time through.    
    mpp[nnn,] <- probres$mp
    effprob[nnn, ] <- probres$effprob
    
    # Check futility
    if (nnn>=nstart) # (only after the burn in period)
       {stopflag <- futility.rule(effprob[nnn,], alpha)}		
  }

  nlast <- nnn # Could be < nmax if we stopped early for futility
  m <- factor(m, level=c(1:16))
  ..("adaptdesign returns $mpp",mpp,1)
  ..("adaptdesign returns $m",m,1)
  ..("adaptdesign returns $effprob",effprob,1)
  ..("adaptdesign returns $nlast",nlast,1)
  ..("adaptdesign returns $stopflag",stopflag,1)
  ..("<<<<<<<<<< adaptdesign",,2)
  invisible(list(mpp=mpp, m=m, effprob=effprob, nlast=nlast,stopflag=stopflag))
}

# Called as: prob(n.sim, c, xcomp, ycomp, nnn, delta, A)
prob <- function(n.sim, c, x, y, nn, delta, A)
{
  ..(">>>>>>>>>> prob: Evaluating the posterior",,1)
  ..("normalizing constant for prior density, i.e., normlized fraction each treatment combination won, of n.sim ( c )",c,1)
  ..("design matrix (x from xcomp)",x,2)
  ..("outcomes (y from ycomp)",y,2)
  ..("Current sample size (nn from nnn)",nn,2)
  ..("delta",delta,2)
  ..("A",A,2)
   
  md <- vector('numeric', 16) # the marginal density proportion term
  mp <- vector("numeric",16) # the posterior probability for each model 
  xx <- x # Why is this necessary???
  p <- 16

  ..("Parameters for multivariate t distribution (p.357):",,2)
  V <- solve(t(xx)%*%xx+1/g*diag(p)) 
  ..("V <- solve.invert(t(xx)%*%xx+1/g*diag℗)",V,2)
  S <- as.numeric(t(y)%*%y-t(y)%*%xx%*%V%*%t(xx)%*%y)
  ..("S <- as.numeric(t(y)%*%y-t(y)%*%xx%*%V%*%t(xx)%*%y)",S,2)
  mu <- V%*%t(xx)%*%y
  ..("mu <- V%*%t(xx)%*%y",mu,2)
  
  ..("Getting n.sim samples from the multivariate t parameterized as above:",,2)
  sample <- rmvt(n.sim, V*S/nn, df=nn) 
  ..("sample <- rmvt(n.sim, V*S/nn, df=nn)",sample,2)
  # Add mu to each row (second arg = 1)
  sample <- t(apply(sample, 1, function(b) {b+mu}))
  ..("sample <- t(apply(sample, 1, function(b) {b+mu}))",sample,2)
  
  ..("Calculating the theta (treatment effect) for each model, matrix(n.sim*16)",,2)
  temp <- sample%*%t(A) 
  ..("temp <- sample%*%t(A)",temp,2)

  # Store which treatment is best for each simulation
  # See explanatory block comment below:
  rankr <- factor(apply(temp, 1, which.min), level=c(1:16))
  ..("rankr <- factor(apply(temp, 1, which.min), level=c(1:16))",rankr,1)
  r <- as.vector(table(rankr))
  ..("r <- as.vector(table(rankr))",r,1)
  r <- r/n.sim # first normalization by number of samples
  ..("r <- r/n.sim",r,1)
  md <- 1/c*r # invert and second nomralize (c) to end up with the marginal density
  # md is eq.2 (p. 355): "the marginal density of the observed data"
  # also described in in the eq. at the bottom of the righthand col. of pg. 356.
  ..("marginal density (md <- 1/c*r)",md,1)
  mp <- md/sum(md) # final normalization
  # ??? I think that this (mp) is the posterior by eq.5 (p. 355)
  ..("normalized md => marginal probability (mp <- md/sum(md))",mp,1)
  # Calculate the posterior probability pr(each treatment < placebo - delta)
  # Compare theta(all treatments)
  rankpost <- t(apply(temp, 1, function(b){return(as.numeric(b < (b[1] - delta)))}))
  # Store pr(each treatment < placebo-delta)
  rankprob <- colMeans(rankpost)
  ..("prob() returns posteriors: $mp",mp,1)
  ..("prob() returns posterior prob: (rankprob aka $effprob)",rankprob,1)
  ..("<<<<<<<<<< prob",,1)
  invisible(list(mp=mp, effprob=rankprob))
}

# =============================================================================
# Notes on the ranking algorithm, from above:
# This is based upon a tricky way to get the count of each min position.
# The which.min apply gives us a vector of the min position in each row (with second arg=1)
# Then the factor against all possible values (c:(1:16)) ends up with assigning each min position to a level
# Then the table of that gives you a count for each level, and finally the as.vector
# returns it to a c(..) vector, which we then normalize, marginalize, and finally normralize that:
# ==========  rankr <- factor(apply(temp, 1, which.min), level=c(1:16))  = "
# [1] 3  15 2  15 15 2  2  16 3  13 13 2  2  2  2  13 2  2  14 15 14 14 13 2  15 15 15 2  15 2  15 9  15 14 15 2  16 2
# ..
# [991] 2  10 13 2  2  14 15 15 2  15 (This will be as long as n.sim)
# Levels: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
# =========== The count for each level: r <- as.vector(table(rankr))
# [1]   1 370  29   4   0  10   2   8   2   7   0   0 153  51 250 113
# ========== fractionalized (normalized) by n.sim: r <- r/n.sim
# [1] 0.001 0.370 0.029 0.004 0.000 0.010 0.002 0.008 0.002 0.007 0.000 0.000 0.153 0.051 0.250 0.113
# ==========  marginalized (actually inverted) and adjusted by c normalizing constant density: md <- 1/c*r 
# [1] 0.01587302 7.55102041 0.34117647 0.05479452 0.00000000 0.16949153 0.03389831 0.14285714 0.02777778 ..
# ==========  normalized md => marginal probability: mp <- md/sum(md))
# [1] 0.0009137485 0.4346832005 0.0196402171 0.0031543098 0.0000000000 0.0097569752 0.0019513950 0.0082237362 0.0015990598 0.0068298826
# [11] 0.0000000000 0.0000000000 0.1314570373 0.0553938459 0.2321215870 0.0942750051
# So at the end this will be 16 long, and sum to 1.0...at least one hopes!
# =============================================================================

futility.rule <- function(effprob, alpha)
{
  ..(">>>>>>>> futility.rule",,2)
  ..("effprob",effprob,2)
  ..("alpha",alpha,2)
  bestind <- which.max(effprob[-1])+1
  ..("best treatment (bestind <- which.max(effprob[-1])+1)", bestind,2)
  ..("prob best treatment (effprob[bestind])", effprob[bestind],2)
  stopflag <- ifelse(effprob[bestind]<alpha, 1, 0) 
  ..("<<<<<<<< futility.rule returns (stopflag <- ifelse(effprob[bestind]<alpha), 1, 0)",stopflag,2)
  return(stopflag)
}

# To trace nothing, set tracelevel=-1

.. <- function(Text,Val="**NOVALUE**",tracelevellimit=0)
{
	if (TraceLevel >= tracelevellimit)
	{
  	    if(Val != "**NOVALUE**") {cat("========== ")}
		cat(Text)
    	if(Val != "**NOVALUE**")
    	{
    		if(((is.null(dim(Val)))||(max(dim(Val))<=NSimDisplayLimit)))
    			{
    			cat(":")
    			if ((!is.null(dim(Val)))||((is.null(dim(Val)))&&(length(Val)>1))){cat("\n")}
    		 	print(Val)
    		 	}
    		else {cat("  (dim: c(",dim(Val),") too big; limit=",NSimDisplayLimit,")\n")}
    	}
    	else {cat("\n")}
	}
	}
