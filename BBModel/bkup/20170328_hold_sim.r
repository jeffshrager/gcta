#cluster.id = as.integer(commandArgs()[3])
cluster.id = 2 # Lorenzo says that this can be any integer -- I'm assuming that it's a run-id, and serves as a random seed for the run.
nmajorloops = 3 # 300 for a full run (but that's extremely slow!)

set.seed(as.integer(cluster.id))

#++--
#1010
#T
#C

for(stst in 1:5){

  compa<-function(E){
    print("----->>>>> Entering Compa(E)")
    print(E)
    xx=rep(0,sum(E[2,])+2)
    for(i in 1:(sum(E[2,])+2)){
      x=log(1/(sum(E[2,])+2))
      if(i<=(i+E[1,1]-1)){x=x+sum(log(c(i:(i+E[1,1]-1))))}
      if((sum(E[2,])+2-i)<=(sum(E[2,])+2-i+E[1,2]-1)){x=x+sum(log((c((sum(E[2,])+2-i+1):(sum(E[2,])+2-i+E[1,2]-1+1)))))}
      if(sum(E[1,])>0){x=x-sum(log((sum(E[2,])+2+1):(sum(E[2,])+2+sum(E[1,]-1+1))))}
      xx[i]=x
      }
    if(E[1,1]>0&E[1,2]>0){xx=xx+sum(log(1:sum(E[1,])))-sum(log(1:sum(E[1,1])))-sum(log(1:sum(E[1,2])))}
    xx=xx+log(sum(E[1,])+1)
    print("< Compa returns:")
    result = sum(exp(xx)[1:(E[2,1]+1)])
    print(result)
    print("          <<<<<---------")
    return(result)
    }

  E=array(3,c(2,2))
  compa(E)

  E=array(3,c(2,2))
  E[1,1]=7
  compa(E)

  mean(rbeta(10000000,8,4)<rbeta(10000000,4,4)) # ??? This doesn't seem to do anything, but take a lot of time! ???

  compa2<-function(E){
    print("----->>>>> Entering Compa2(E)")
    print(E)
    a1=1-compa(E[,1:2])
    a2=1-compa(E[,3:4])
    a=c(a1*a2,a1*(1-a2),a2*(1-a1),(1-a1)*(1-a2))*c(0.5*0.5,0.5*0.5,0.5*0.01,0.5*0.99)
    a=a/sum(a)
    a1=a[1]+a[2]
    a2=a[1]+a[3]
    a3=1-a[4]
    result = 5*(a3)^5.7+a1^5.7+a2^5.7
    print("< Compa2 returns:")
    print(result)
    print("          <<<<<---------")
    return(result)
    }

  E=cbind(E,E)
  compa2(E)

  pred<-function(E){
    print("----->>>>> Entering pred(E)")
    print(E)
    EE=E
    E=EE[,1:2]
    xx=rep(0,sum(E[2,])+2)

    for(i in 1:(sum(E[2,])+2)){
      x=log(1/(sum(E[2,])+2))
      if(i<=(i+E[1,1]-1)){x=x+sum(log(c(i:(i+E[1,1]-1))))}
      if((sum(E[2,])+2-i)<=(sum(E[2,])+2-i+E[1,2]-1)){x=x+sum(log((c((sum(E[2,])+2-i+1):(sum(E[2,])+2-i+E[1,2]-1+1)))))}
      if(sum(E[1,])>0){x=x-sum(log((sum(E[2,])+2+1):(sum(E[2,])+2+sum(E[1,]-1+1))))}
      xx[i]=x
      }

    if(E[1,1]>0&E[1,2]>0){xx=xx+sum(log(1:sum(E[1,])))-sum(log(1:sum(E[1,1])))-sum(log(1:sum(E[1,2])))}
    xx=xx+log(sum(E[1,])+1)
    xx1=exp(xx)

    E=EE[,3:4]
    xx=rep(0,sum(E[2,])+2)

    for(i in 1:(sum(E[2,])+2)){
      x=log(1/(sum(E[2,])+2))
      if(i<=(i+E[1,1]-1)){x=x+sum(log(c(i:(i+E[1,1]-1))))}
      if((sum(E[2,])+2-i)<=(sum(E[2,])+2-i+E[1,2]-1)){x=x+sum(log((c((sum(E[2,])+2-i+1):(sum(E[2,])+2-i+E[1,2]-1+1)))))}
      if(sum(E[1,])>0){x=x-sum(log((sum(E[2,])+2+1):(sum(E[2,])+2+sum(E[1,]-1+1))))}
      xx[i]=x
      }

    if(E[1,1]>0&E[1,2]>0){xx=xx+sum(log(1:sum(E[1,])))-sum(log(1:sum(E[1,1])))-sum(log(1:sum(E[1,2])))}
    xx=xx+log(sum(E[1,])+1)
    xx2=exp(xx)

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

    result=c(a1,a2)
    print("< Pred returns:")
    print(result)
    print("          <<<<<---------")
    return(result)
    } # Close pred()

  pred(E) # ???? Again, the result of this doesn't go anywhere ????

  pred2<-function(E){

    EE=E
    EE[1,]=E[2,]
    EE[2,]=E[1,]
    E=EE
    E=EE[,1:2]
    xx=rep(0,sum(E[2,])+2)

    for(i in 1:(sum(E[2,])+2)){
      x=log(1/(sum(E[2,])+2))
      if(i<=(i+E[1,1]-1)){x=x+sum(log(c(i:(i+E[1,1]-1))))}
      if((sum(E[2,])+2-i)<=(sum(E[2,])+2-i+E[1,2]-1)){x=x+sum(log((c((sum(E[2,])+2-i+1):(sum(E[2,])+2-i+E[1,2]-1+1)))))}
      if(sum(E[1,])>0){x=x-sum(log((sum(E[2,])+2+1):(sum(E[2,])+2+sum(E[1,]-1+1))))}
      xx[i]=x
      }

    if(E[1,1]>0&E[1,2]>0){xx=xx+sum(log(1:sum(E[1,])))-sum(log(1:sum(E[1,1])))-sum(log(1:sum(E[1,2])))}
    xx=xx+log(sum(E[1,])+1)
    xx1=exp(xx)

    E=EE[,3:4]
    xx=rep(0,sum(E[2,])+2)

    for(i in 1:(sum(E[2,])+2)){
      x=log(1/(sum(E[2,])+2))
      if(i<=(i+E[1,1]-1)){x=x+sum(log(c(i:(i+E[1,1]-1))))}
      if((sum(E[2,])+2-i)<=(sum(E[2,])+2-i+E[1,2]-1)){x=x+sum(log((c((sum(E[2,])+2-i+1):(sum(E[2,])+2-i+E[1,2]-1+1)))))}
      if(sum(E[1,])>0){x=x-sum(log((sum(E[2,])+2+1):(sum(E[2,])+2+sum(E[1,]-1+1))))}
      xx[i]=x
      }

    if(E[1,1]>0&E[1,2]>0){xx=xx+sum(log(1:sum(E[1,])))-sum(log(1:sum(E[1,1])))-sum(log(1:sum(E[1,2])))}
    xx=xx+log(sum(E[1,])+1)
    xx2=exp(xx)

    v1=c(1:(sum(EE[2,1:2])+2))%*%t(rep(1,(sum(EE[2,3:4])+2)))
    v2=t(c(1:(sum(EE[2,3:4])+2))%*%t(rep(1,(sum(EE[2,1:2])+2))))
    x=xx1%*%t(xx2)

    x=x*(v1>(EE[2,1]+1))*(v2>(EE[2,3]+1))*0.5*0.99+
      x*(v1>(EE[2,1]+1))*(v2<=(EE[2,3]+1))*0.5*0.01+
      x*(v1<=(EE[2,1]+1))*(v2>(EE[2,3]+1))*0.5*0.5+
      x*(v1<=(EE[2,1]+1))*(v2<=(EE[2,3]+1))*0.5*0.5
    x=x/sum(x)

    a1=sum(x*(c((1+EE[1,1]):(1+EE[1,1]+EE[2,1]+EE[2,2]+1))/((EE[1,1]+EE[1,2]+EE[2,1]+EE[2,2]+3))))
    a2=sum(t(x)*(c((1+EE[1,1+2]):(1+EE[1,1+2]+EE[2,1+2]+EE[2,2+2]+1))/((EE[1,1+2]+EE[1,2+2]+EE[2,1+2]+EE[2,2+2]+3))))

    return(c(a1,a2))
    }

  ####################
  ####################

  PAR1=0.4
  #runif(1,0,0.8)

  ####################
  ####################

    W1=rep(0,2000)
    WW1=rep(0,2000)
    RTTTT<-array(0,c(2000,400,3,2,4))
    MMMM<-array(0,c(2000,400,3))

    for(kk in 1:nmajorloops){
      R<-array(0,c(3,2,4))
      SS=400
      RTT<-array(0,c(SS,3,2,4))
      v2=rep(0,3)

      for(ii in 1:SS){

        # Create the mutation array, as 1/0 1/0 1/0
        mk=c(rbinom(3,1,0.5)) 
        mk[3]=rbinom(3,1,0.5-0.05*(as.integer(cluster.id)-1)) # !!! This seems to be wrong (or the above line is) ???

        MMMM[kk,ii,]=mk
        u=rep(0,4)

        for(i in 1:3){
          E=R[i,,]
          u[i]=compa2(E)
          v=pred(E)
          v=v[2-mk[i]]
          E1=E
          E2=E
          E1[1,1+2*(1-mk[i])]=E1[1,1+2*(1-mk[i])]+1
          E2[1,2+2*(1-mk[i])]=E2[1,2+2*(1-mk[i])]+1
          u[i]=-(u[i]-v*compa2(E1)-(1-v)*compa2(E2))
          v2[i]=v
          }

        u[4]=0

        for(i in 1:3){
          E=R[i,,]
          uu=compa2(E)
          v=pred2(E)
          v=v[2-mk[i]]
          E1=E
          E2=E
          E1[2,1+2*(1-mk[i])]=E1[2,1+2*(1-mk[i])]+1
          E2[2,2+2*(1-mk[i])]=E2[2,2+2*(1-mk[i])]+1
          u[4]=u[4]-(uu-v*compa2(E1)-(1-v)*compa2(E2))
          }

        Rnd=rep(0,4)
        Und=0

        for(i in 1:100){
          kkkk=0
          while(kkkk<3){
            rand=runif(4,0,1)
            kkkk=0
            for(i1 in 1:3){for(i2 in i1:3){if(i1<i2){if( ((rand[i1]-rand[i2])*(v2[i1]-v2[i2]))>=0){ kkkk=kkkk+1 }}}}
            }
          rand= rand/sum(rand)
          if(sum(rand*u)>Und){Rnd=rand;Und=sum(rand*u)}
	  }

        Rnd=Rnd/sum(Rnd)
        if(!(stst==1 | stst==3 | stst==5)){Rnd= Rnd^0;Rnd= Rnd/sum(Rnd)}

        TR=sum(rmultinom(1,1,Rnd)*c(1:4))
        ou=rbinom(1,1,PAR1)

        if(TR<4){
          if(stst==1){if(TR==1 & mk[TR]==1){ou=rbinom(1,1,(PAR1+0.2))}}
          if(stst>1){if(TR==1 ){ou=rbinom(1,1,(PAR1+0.2))}}
          R[TR,1,2*(1-mk[TR])+(1-ou)+1]=R[TR,1,2*(1-mk[TR])+(1-ou)+1]+1
          }

        if(TR==4){for(i in 1:3){R[i,2,2*(1-mk[i])+(1-ou)+1]=R[i,2,2*(1-mk[i])+(1-ou)+1]+1}}

        RTT[ii,,,]=R
        }

      RTTTT[kk,,,,]=RTT
      E=R[1,,1:2]
      mm=fisher.test(E,alternative="greater")
      pv=as.vector((mm[1])$p.value)
      W1[kk]=pv
      WW1[kk]=1-compa(R[1,,1:2])

      print(kk) # This is just a tracer, I think.

      resultsdir = "/Users/jeffshrager/Desktop/cc/gcta/BrianAlexaderSim/results/"
      if(stst==1){save(list=ls(), file = paste0(resultsdir,"eigthA.N.",as.integer(cluster.id),".RData"))}
      if(stst==2){save(list=ls(), file = paste0(resultsdir,"eigthB.N.",as.integer(cluster.id),".RData"))}
      if(stst==3){save(list=ls(), file = paste0(resultsdir,"eigthA1.N.",as.integer(cluster.id),".RData"))}
      if(stst==4){save(list=ls(), file = paste0(resultsdir,"eigthB1.N.",as.integer(cluster.id),".RData"))}
      if(stst==5){save(list=ls(), file = paste0(resultsdir,"eigthC.N.",as.integer(cluster.id),".RData"))}

      } # Close: for(kk in 1:nmajorloops)...
    } # Close: for(stst in 1:5)...
