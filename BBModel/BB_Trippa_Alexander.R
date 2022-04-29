# Below is a short note on several important variables/concepts. This info
# is repeated later in the relevant sections. See 257 and onward for output.
# SS = number of subjects
# W1 = vector to contain p values for each simulation 
# WW1 = vector to contain result of 1-compa(R[1,,1:2]) for each simulation
# RTTTT = 5 dimensional array of patient info for each simulation run
# MMMM = 3 dimensional array of biomarker info for each patient in each simulation
#   i.e. MMMM[2,1,] gives the patient's biomarker signature
# R = contains a single instance of the "outcome table" r
# RTT = "outcome table" for each patient

# the outcome table appears as such (from paper):
#                     +       _
# experimental arm  0   0   1   1
#      control arm  0   0   1   1
#
# columns 1 and 2 represent positive instance of this sub-biomarker
# columns 3 and 4 represent negative instance of this sub-biomarker
# columns 1 and 3 represent positive outcome
# columns 2 and 4 represent negative outcome
# there are 3 of these tables for each "biomarker signature", since it consists
# of three sub-biomarkers (e.g. a biomarker, or biomarker signature, 
# is shown as [+ - +] or [1 0 1])

# To trace nothing, set tracelevel=-1

TraceLevel = 0
NSimDisplayLimit = 10

.. <- function(Text,Val="**NOVALUE**",tracelevellimit=0)
{
	if (TraceLevel >= tracelevellimit)
	{
  	    if(Val != "**NOVALUE**") {cat(".. ")}
		cat(Text)
    	if(Val != "**NOVALUE**")
    	{
    		if(((is.null(dim(Val)))||(max(dim(Val))<=NSimDisplayLimit)))
    			{
    			cat("=")
    			if ((!is.null(dim(Val)))||((is.null(dim(Val)))&&(length(Val)>1))){cat("\n")}
    		 	print(Val)
    		 	}
    		else {cat("  (dim: c(",dim(Val),") too big; limit=",NSimDisplayLimit,")\n")}
    	}
    	else {cat("\n")}
	}
	}

#cluster.id = as.integer(commandArgs()[3])

cluster.id = 1 # just set it to 1 to get the script to run

set.seed(as.integer(cluster.id))

#++--
#1010
#T
#C

for(stst in 1:5){
  ..("*************** MASTER LOOP stst:",stst)  
  # Compa is called betacomparison in the paper: "Computation of the probability of a Beta(a,b) 
  # random variable to be larger than an independent Beta(a',b') random variable."
  # Compa is used in compa2, it takes in an array and calculates the logarithm of specific values
  # before taking the exponent of these calculations, and summing up a slice of them. 

  # Here's a whole compa i/o series:
  # ========== Compa >> entering with E:
  #      [,1] [,2]
  # [1,]    2    2
  # [2,]    3    3
  # ========== Compa xx at the top is:
  # [1] 0 0 0 0 0 0 0 0
  # ========== Compa: xx at the end is:
  # [1] -3.091042 -2.243745 -1.838279 -1.663926 -1.663926 -1.838279 -2.243745
  # [8] -3.091042
  # ========== Compa returns:[1] 0.5

  compa<-function(E){
    ..("Compa >> entering with E",E)
    # E is 2x2, starts out as 3,3,3,3 and then 7,3,3,3, and from there on is all 0s with one 1 sometimes.
    # These values keep increasing. 
    # ?? I guess these are the params of a beta function, but why 2x2 when Beta on needs 2 params ??
    # xx is  a vector of length 2 greater than the sum of the elements in the second row of E.
    # ?? What will xx end up representing ?? What's the second row ?? (What's the first row??)??
    xx=rep(0,sum(E[2,])+2) # create a vector of certain length filled with 0s
    ..("Compa xx at the top is",xx)
    for(i in 1:(sum(E[2,])+2)){ # iterate thru the vector
      x=log(1/(sum(E[2,])+2))
      if(i<=(i+E[1,1]-1)){ # always true unless E[1,1] <= 0
        x=x+sum(log(c(i:(i+E[1,1]-1))))
      } 
      if((sum(E[2,])+2-i)<=(sum(E[2,])+2-i+E[1,2]-1)){ # always true unless E[1,2] <= 0
        x=x+sum(log((c((sum(E[2,])+2-i+1):(sum(E[2,])+2-i+E[1,2]-1+1)))))
      } 
      if(sum(E[1,])>0){ # always true unless sum(E[1,]) <= 0
        x=x-sum(log((sum(E[2,])+2+1):(sum(E[2,])+2+sum(E[1,]-1+1))))
      } 
      xx[i]=x
    }
    if(E[1,1]>0&E[1,2]>0){ # always true?
      xx=xx+sum(log(1:sum(E[1,])))-sum(log(1:sum(E[1,1])))-sum(log(1:sum(E[1,2])))
    } 
    xx=xx+log(sum(E[1,])+1) 
    ..("Compa: xx at the end is",xx)
    result=sum(exp(xx)[1:(E[2,1]+1)])
    ..("Compa returns result",result)
    return(result) #return sum of exponents for a slice of the vector
  }
  
  # ?? Are the in line calls just for testing, or do they serve a function ??
  E=array(3,c(2,2)); 
  compa(E)
  E=array(3,c(2,2)); 
  E[1,1]=7; 
  compa(E) 

  # rbeta generates a random variate for the beta distribution of alpha and beta given
  # This code is standalone and I'm not sure where or if it's used elsewhere, the expression left of 
  # operator will be usually greater than that to the left
  # mean(rbeta(10000000,8,4)<rbeta(10000000,4,4))
  
  # called comparisonPositiveNegative in the paper
  # Description from paper: Posterior probabilities with the described model of treatment effects in the biomarker-positive and biomarker-negative
  # subpopulations, assuming a = b = 1
  # uses compa function (above)
  # Note that this function outputs a single value, while comparisonPositiveNegative (from the paper) returns
  # a vector of two values. 

for the control!

  compa2<-function(E){
    ..("Compa2 arg E",E)
    a1=1-compa(E[,1:2]) # probability of a Beta(a,b) random variable to be less than an independent Beta(a',b') random variable
    a2=1-compa(E[,3:4]) # probability of a Beta(a,b) random variable to be less than an independent Beta(a',b') random variable
    a=c(a1*a2,a1*(1-a2),a2*(1-a1),(1-a1)*(1-a2))*c(0.5*0.5,0.5*0.5,0.5*0.01,0.5*0.99) # comments below on what these calcs may mean
    ..("Compa2 a",a)
    # a1*a2 = prob. both Beta(a,b) random variables are less than independent Betas
    # a1*(1-a2) = prob. first Beta RV is less than ind. Beta RV and that the second Beta RV is greater than ind. Beta RV
    # a2*(1-a1) = opposite of the above
    # (1-a1)*(1-a2) = prob. both Beta RVs are greater than ind. Beta RVs
    # 0.333 = Pr(Ek+)
    # 0.667 = Pr(NEk+)
    # 0.5 = Pr(Ek-|Ek+) = Pr(NEk-|Ek+) 
    # 0.01 ~ 0 = Pr(Ek-|NEk+) 
    # 0.99 ~ 1 = Pr(NEk-|NEk+)
    # Ek+/- = effect in treatment arm k on biomarker k (biomarker k is pos or neg)
    # NEk+/- = no effect in treatment arm k on biomarker k
    # in the study text they use: a=c(a1*a2,a1*(1-a2),a2*(1-a1),(1-a1)*(1-a2))*c((1/3)*0.5,(1/3)*0.5,(2/3)*0.00,(2/3)*1)
    # note that different fractions are used, i.e. 0.5 instead of 1/3

    a=a/sum(a)
    a1=a[1]+a[2];
    a2=a[1]+a[3]; # essentially just a[1]
    a3=1-a[4]

    # in the paper a1, a2, and a3 are not calculated 
    # in the paper instead they simply return after the a=a/sum(a) calc: return(c(a[1]+a[2],a[1]))
    result = 5*(a3)^5.7+a1^5.7+a2^5.7
    ..("Compa2 returns resutl",result)
    return(result) # the paper returns a vector of two values, not a single value
  }
  
  E=cbind(E,E) # goes from 2x2 to 2x4 # ?? Are these things just for testing, or do they serve a function ??
  compa2(E)
  
  # Called prediction in the paper: "Predictive probabilities with the
  # described model of positive outcome for the next patient, assuming
  # the patient will be assigned to the experimental arm, in
  # biomarker-positive group and in the biomarker-negative group,
  # assuming a = b = 1."

  pred<-function(E){EE=E
    ..("Pred arg EE",EE)
    E=EE[,1:2]
    ..("Pred arg E after EE (minor) rearrangement",E)
    xx=rep(0,sum(E[2,])+2)
    for(i in 1:(sum(E[2,])+2)){
      x=log(1/(sum(E[2,])+2))
      if(i<=(i+E[1,1]-1)){
        x=x+sum(log(c(i:(i+E[1,1]-1))))
      }
      if((sum(E[2,])+2-i)<=(sum(E[2,])+2-i+E[1,2]-1)){
        x=x+sum(log((c((sum(E[2,])+2-i+1):(sum(E[2,])+2-i+E[1,2]-1+1)))))
      }
      if(sum(E[1,])>0){
        x=x-sum(log((sum(E[2,])+2+1):(sum(E[2,])+2+sum(E[1,]-1+1))))
      }
      xx[i]=x
    }
    if(E[1,1]>0&E[1,2]>0){
      xx=xx+sum(log(1:sum(E[1,])))-sum(log(1:sum(E[1,1])))-sum(log(1:sum(E[1,2])))
    }
    xx=xx+log(sum(E[1,])+1)
    # same as compa up to this point (compa instead returns a slice of the exponents)
    xx1=exp(xx)
  
    # nearly the same as compa but for different peice of E
    E=EE[,3:4]
    xx=rep(0,sum(E[2,])+2)
    for(i in 1:(sum(E[2,])+2)){
      x=log(1/(sum(E[2,])+2))
      if(i<=(i+E[1,1]-1)){
        x=x+sum(log(c(i:(i+E[1,1]-1))))
      }
      if((sum(E[2,])+2-i)<=(sum(E[2,])+2-i+E[1,2]-1)){
        x=x+sum(log((c((sum(E[2,])+2-i+1):(sum(E[2,])+2-i+E[1,2]-1+1)))))
      }
      if(sum(E[1,])>0){
        x=x-sum(log((sum(E[2,])+2+1):(sum(E[2,])+2+sum(E[1,]-1+1))))
      }
      xx[i]=x
    }
    if(E[1,1]>0&E[1,2]>0){
      xx=xx+sum(log(1:sum(E[1,])))-sum(log(1:sum(E[1,1])))-sum(log(1:sum(E[1,2])))
    }
    xx=xx+log(sum(E[1,])+1)
    xx2=exp(xx) # this line different from compa
    # block above same as compa but different slice of E and different exponent calculated
  
INSIGHT trial

    # code below is the first new code (not the same as compa)
    v1=c(1:(sum(EE[2,1:2])+2))%*%t(rep(1,(sum(EE[2,3:4])+2)))
    v2=t(c(1:(sum(EE[2,3:4])+2))%*%t(rep(1,(sum(EE[2,1:2])+2))))
    x=xx1%*%t(xx2)
  
    x=x*(v1>(EE[2,1]+1))*(v2>(EE[2,3]+1))*0.5*0.5+
      x*(v1>(EE[2,1]+1))*(v2<=(EE[2,3]+1))*0.5*0.5+
      x*(v1<=(EE[2,1]+1))*(v2>(EE[2,3]+1))*0.5*0.01+
      x*(v1<=(EE[2,1]+1))*(v2<=(EE[2,3]+1))*0.5*0.99

    x=x/sum(x)
  
    a1=sum(x*(c((1+EE[1,1]):(1+EE[1,1]+EE[2,1]+EE[2,2]+1))/((EE[1,1]+EE[1,2]+EE[2,1]+EE[2,2]+3))))
    a2=sum(t(x)*(c((1+EE[1,1+2]):(1+EE[1,1+2]+EE[2,1+2]+EE[2,2+2]+1))/((EE[1,1+2]+EE[1,2+2]+EE[2,1+2]+EE[2,2+2]+3))))
    ..("Pred returns c(a1,a2)",c(a1,a2))  
    return(c(a1,a2))
  } # Close Pred
  
  pred(E) # ?? I don't get why these are called along the way ??
  
  # Pred2 is similar to pred, the difference being that first row of
  # EE becomes the same as the second row. ????

It is the prediction for the control.

  pred2<-function(E){
    ..("Pred2 arg E",E)
    EE=E; 
    EE[1,]=E[2,];
    EE[2,]=E[1,];
    E=EE
    E=EE[,1:2] # ?? WTF! ??
    ..("Pred2 E after EE rearragement, E",E)  
    ..("Pred2 EE after EE rearragement, EE",EE)  
    xx=rep(0,sum(E[2,])+2)
    for(i in 1:(sum(E[2,])+2)){
      x=log(1/(sum(E[2,])+2))
      if(i<=(i+E[1,1]-1)){
        x=x+sum(log(c(i:(i+E[1,1]-1))))
      }
      if((sum(E[2,])+2-i)<=(sum(E[2,])+2-i+E[1,2]-1)){
        x=x+sum(log((c((sum(E[2,])+2-i+1):(sum(E[2,])+2-i+E[1,2]-1+1)))))
      }
      if(sum(E[1,])>0){
        x=x-sum(log((sum(E[2,])+2+1):(sum(E[2,])+2+sum(E[1,]-1+1))))
      }
      xx[i]=x
    }
    if(E[1,1]>0&E[1,2]>0){
      xx=xx+sum(log(1:sum(E[1,])))-sum(log(1:sum(E[1,1])))-sum(log(1:sum(E[1,2])))
    }
    xx=xx+log(sum(E[1,])+1)
    xx1=exp(xx)
  
    E=EE[,3:4]
    xx=rep(0,sum(E[2,])+2)
    for(i in 1:(sum(E[2,])+2)){
      x=log(1/(sum(E[2,])+2))
      if(i<=(i+E[1,1]-1)){
        x=x+sum(log(c(i:(i+E[1,1]-1))))
      }
      if((sum(E[2,])+2-i)<=(sum(E[2,])+2-i+E[1,2]-1)){
        x=x+sum(log((c((sum(E[2,])+2-i+1):(sum(E[2,])+2-i+E[1,2]-1+1)))))
      }
      if(sum(E[1,])>0){
        x=x-sum(log((sum(E[2,])+2+1):(sum(E[2,])+2+sum(E[1,]-1+1))))
      }
      xx[i]=x
    }
    if(E[1,1]>0&E[1,2]>0){
      xx=xx+sum(log(1:sum(E[1,])))-sum(log(1:sum(E[1,1])))-sum(log(1:sum(E[1,2])))
    }
    xx=xx+log(sum(E[1,])+1)
    xx2=exp(xx)
  
    v1=c(1:(sum(EE[2,1:2])+2))%*%t(rep(1,(sum(EE[2,3:4])+2)))
    v2=t(c(1:(sum(EE[2,3:4])+2))%*%t(rep(1,(sum(EE[2,1:2])+2))))
    x=xx1%*%t(xx2)
  
These are parmameters of the prior (prob of a tx in the biomarker + and cond on that there's a p of effecet in neg.

    x=x*(v1>(EE[2,1]+1))*(v2>(EE[2,3]+1))*0.5*0.99+         prob of suc in b.m. + = 0.5
      x*(v1>(EE[2,1]+1))*(v2<=(EE[2,3]+1))*0.5*0.01+        prob of - ...
      x*(v1<=(EE[2,1]+1))*(v2>(EE[2,3]+1))*0.5*0.5+
      x*(v1<=(EE[2,1]+1))*(v2<=(EE[2,3]+1))*0.5*0.5

    x=x/sum(x)
  
    a1=sum(x*(c((1+EE[1,1]):(1+EE[1,1]+EE[2,1]+EE[2,2]+1))/((EE[1,1]+EE[1,2]+EE[2,1]+EE[2,2]+3))))
    a2=sum(t(x)*(c((1+EE[1,1+2]):(1+EE[1,1+2]+EE[2,1+2]+EE[2,2+2]+1))/((EE[1,1+2]+EE[1,2+2]+EE[2,1+2]+EE[2,2+2]+3))))
    ..("Pred2 returns c(a1,a2)",c(a1,a2))
    return(c(a1,a2))
  } # Close pred2
  
  PAR1=0.4 # ??????????????????????? Used three times below 
  #runif(1,0,0.8)
  
  # These four variables contain the results of each run. 
  # ?? Where does the 2000 comes from ?? There appear to be 400 subject. 

  W1=rep(0,2000) # p values for each sim
  WW1=rep(0,2000) # result of 1-compa(R[1,,1:2]) per sim
  RTTTT<-array(0,c(2000,400,3,2,4)); # 5 dimensional array of patient info per sim
  MMMM<-array(0,c(2000,400,3)); # 3 dimensional array of biomarkers per patient per sim
        # i.e. MMMM[2,1,] gives the patient's biomarker signature
  
  # the outcome table appears as such (from paper):
  #                     +       _
  # experimental arm  0   0   1   1
  #      control arm  0   0   1   1
  #
  # columns 1 and 2 represent positive instance of this sub-biomarker
  # columns 3 and 4 represent negative instance of this sub-biomarker
  # columns 1 and 3 represent positive outcome
  # columns 2 and 4 represent negative outcome
  # there are 3 of these tables for each "biomarker signature", since it consists
  # of three sub-biomarkers (e.g. a biomarker, or biomarker signature, 
  # is shown as [+ - +] or [1 0 1])
  
  #SS=400 # number of subjects
  SS=4 # number of subjects # <===================================================================================== TESTING ONLY !!!!!!!!!!!!
  ncycles = 3 # Jeff might have added this, extracting it from the for loop (for no obvious reason -- there's LOTS of weirder constants around!)
  ..("************************************ WARNING SS=4 for testing ************************************ ",SS)
  for(kk in 1:ncycles){ # kk is the cycle number
    ..("********** MAIN INNER CYCLE kk",kk)
    R<-array(0,c(3,2,4)) # on the first pass R consists of zeros, contains a single instance of the "outcome table" 
    RTT<-array(0,c(SS,3,2,4)) # the "outcome table" for each patient
    v2=rep(0,3)
    for(ii in 1:SS){
      ..("***** SUBJECT CYCLE ii",ii)
      mk=c(rbinom(3,1,0.5)); # generating the biomarkers
      mk[3]=rbinom(3,1,0.5-0.05*(as.integer(cluster.id)-1)); # not sure what this line does ??
      ..("  biomarkers for this subject (mk)",mk)
      MMMM[kk,ii,]=mk # storing the biomarkers in MMMM
      u=rep(0,4) # uncertainty reduction???
      for(i in 1:3){
        ..("inner loop i",i)
        E=R[i,,]
        u[i]=compa2(E)
        v=pred(E);
        v=v[2-mk[i]];
        E1=E;
        E2=E
        E1[1,1+2*(1-mk[i])]=E1[1,1+2*(1-mk[i])]+1
        E2[1,2+2*(1-mk[i])]=E2[1,2+2*(1-mk[i])]+1
        u[i]=-(u[i]-v*compa2(E1)-(1-v)*compa2(E2))
        v2[i]=v
	}
	 ..("At end of first of paired loops, E",E)
	 ..("  E1",E1)
	 ..("  E1",E2)
	 ..("  u",u)
	 ..("  v",v)
	 ..("  v2",v2)

       u[4]=0
       for(i in 1:3){
         ..("Inner 1:3 loop with i",i)
         E=R[i,,]
         ..("  and E",E)

         uu=compa2(E)
         v=pred2(E);v=v[2-mk[i]];E1=E;E2=E
         E1[2,1+2*(1-mk[i])]=E1[2,1+2*(1-mk[i])]+1
         E2[2,2+2*(1-mk[i])]=E2[2,2+2*(1-mk[i])]+1
         u[4]=u[4]-(uu-v*compa2(E1)-(1-v)*compa2(E2))}

	 ..("At end of paired 3 loops, E",E)
	 ..("  uu",uu)
	 ..("  v",v)
	 ..("  E1",E1)
	 ..("  E1",E2)
                                                       
       Rnd=rep(0,4)
       Und=0

       # This loop appears to just randomize rnd and und, but the inner loop seems useless ??

       for(i in 1:100){
         # ?? kkkk is counting something, but is never used -- what's this for ??
         kkkk=0
         while(kkkk<3){rand=runif(4,0,1);
          kkkk=0
          for(i1 in 1:3){
            for(i2 in i1:3){
              if(i1<i2){
                if(((rand[i1]-rand[i2])*(v2[i1]-v2[i2]))>=0){ 
                  kkkk=kkkk+1 
                }
              }
            }
          }
         }
         rand= rand/sum(rand)
         if(sum(rand*u)>Und){
           Rnd=rand;Und=sum(rand*u)
         }
       }

       Rnd=Rnd/sum(Rnd)
       if(!(stst==1 | stst==3 | stst==5)){
         Rnd=Rnd^0
	 Rnd=Rnd/sum(Rnd)
       }
                                                       
       # ?? What does this series of TR-conditioned actions do ???

This decides which arm the patient goes into:

       TR=sum(rmultinom(1,1,Rnd)*c(1:4))
       ou=rbinom(1,1,PAR1)

       ..("TR entering TR condition cascade",TR)

       if(TR<4){
        if(stst==1){
          if(TR==1 & mk[TR]==1){
            ou=rbinom(1,1,(PAR1+0.2))
          }
        }
        if(stst>1){
          if(TR==1 ){
            ou=rbinom(1,1,(PAR1+0.2))
          }
        }
        R[TR,1,2*(1-mk[TR])+(1-ou)+1]=R[TR,1,2*(1-mk[TR])+(1-ou)+1]+1
       }
                                                       
       if(TR==4){
         for(i in 1:3){
           R[i,2,2*(1-mk[i])+(1-ou)+1]=R[i,2,2*(1-mk[i])+(1-ou)+1]+1
         }
       }
                                                       
       RTT[ii,,,]=R # taking R and storing it in RTT for the iith patient
       ..("Final R for this subject", R)
    } # End of per-subject loop
    
    RTTTT[kk,,,,]=RTT # taking RTT and storing it in RTTTT for the kkth simulation
    E=R[1,,1:2]
    mm=fisher.test(E,alternative="greater");
    pv=as.vector((mm[1])$p.value)
    W1[kk]=pv
    WW1[kk]=1-compa(R[1,,1:2])
    print(kk)
    # the file path was changed to this new format from original format on 3/27/17
    if(stst==1){
      save(list=ls(), file = paste0(getwd(),"/eigthA.N.",as.integer(cluster.id),".RData"))
    }
    if(stst==2){
      save(list=ls(), file = paste0(getwd(),"/eigthB.N.",as.integer(cluster.id),".RData"))
    }
    if(stst==3){
      save(list=ls(), file = paste0(getwd(),"/eigthA1.N.",as.integer(cluster.id),".RData"))
    }
    if(stst==4){
      save(list=ls(), file = paste0(getwd(),"/eigthB1.N.",as.integer(cluster.id),".RData"))
    }
    if(stst==5){
      save(list=ls(), file = paste0(getwd(),"/eigthC.N.",as.integer(cluster.id),".RData"))
    }
    
    
  }
}

