library(rjags)
library(R2jags)
library(pse)

## Simulates the competition dynamics
compet1 <- function(rA, kA, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea=0.25*90, A0=50, P0=50, obs=TRUE){
    A <- A0
    P <- P0
    lambA <- lambP <- yA <- yP <- c()
    for(t in 1:(length(dt))){
        muA <- A + (rA * A * (1 - ((A+aPA*P)/kA)))*dt[t]
        muP <- P + (rP * P * (1 - ((P+aAP*A)/kP)))*dt[t]
        A <- ifelse(muA>0, rpois(1, muA),0)
        P <- ifelse(muP>0, rpois(1, muP), 0)
        lambA[t] <- A*sampArea/totArea
        lambP[t] <- P*sampArea/totArea
        nA <- rpois(1,lambA[t])
        nP <- rpois(1,lambP[t])
        yA[t] <- rbinom(1, nA, pA)
        yP[t] <- rbinom(1, nP, pP)
    }
    if(obs)
        cbind(yA,yP)
    else
        cbind(lambA, lambP)
}

## Checking if Nash equilibirum occurs from different combinations of parameters
## For each combination nrep repetitions of the simulated dynamics is ran and then a payoff matrix is built
## from wich Nash equilibrium is checked

sim.Nash <- function(obj1, obj2, obj3, data, dt, sampArea=10*pi*0.25^2,
                     totArea=0.25*90, A0=50, P0=50, ncomb=10, nrep=100, obs=FALSE){
    if(missing(dt))
        dt <- diff(c(0,as.numeric(colnames(data[[1]]))))
    Nash <- c()
    for(k in 1:ncomb){
        LL <- matrix(NA, nrep, 2)
        LH <- matrix(NA, nrep, 2)
        HL <- matrix(NA, nrep, 2)
        HH <- matrix(NA, nrep, 2)
        i <- sample(1:length(obj1$BUGSoutput$sims.list$rA), 1)
        j <- sample(1:length(obj2$BUGSoutput$sims.list$r), 1)
        z <- sample(1:length(obj3$BUGSoutput$sims.list$r), 1)
        rAi <- obj2$BUGSoutput$sims.list$r[j]
        kAi <- obj2$BUGSoutput$sims.list$k[j]
        rPi <- obj3$BUGSoutput$sims.list$r[z]
        kPi <- obj3$BUGSoutput$sims.list$k[z]
        rA <- obj1$BUGSoutput$sims.list$rA[i]
        kA <- obj1$BUGSoutput$sims.list$kA[i]
        rP <- obj1$BUGSoutput$sims.list$rP[i]
        kP <- obj1$BUGSoutput$sims.list$kP[i]
        pA <- obj1$BUGSoutput$sims.list$pA[i]
        pP <- obj1$BUGSoutput$sims.list$pP[i]
        aPA <- obj1$BUGSoutput$sims.list$aPA[i]
        aAP <- obj1$BUGSoutput$sims.list$aAP[i]
        ## Stochastic
        for(j in 1:nrep){
            LL[j,] <- compet1(rA, kA, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea, A0, P0, obs=obs)[length(dt),]>0
            LH[j,] <- compet1(rA, kA, pA, rPi, kPi, pP, aPA, aAP, dt, sampArea, totArea, A0, P0, obs=obs)[length(dt),]>0
            HL[j,] <- compet1(rAi, kAi, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea, A0, P0, obs=obs)[length(dt),]>0
            HH[j,] <- compet1(rAi, kAi, pA, rPi,kPi, pP, aPA, aAP, dt, sampArea, totArea, A0, P0, obs=obs)[length(dt),]>0
        }
        #browser()
        LLp <- apply(LL,2,sum)/nrep
        LHp <- apply(LH,2,sum)/nrep
        HLp <- apply(HL,2,sum)/nrep
        HHp <- apply(HH,2,sum)/nrep
        Nash[k] <-
            LLp[1]>=HLp[1]&LHp[1]>=HHp[1] &
            #LLp[1]>LHp[1]&HLp[1]>HHp[1] &
            #LLp[2]>HLp[2]&LHp[2]>HHp[2] &
            LLp[2]>=LHp[2]&HLp[2]>=HHp[2]
        ## Lotka-Volterra, no stochasticity
        ## LLe <- c( (kA - aAP*kP)/(1-aAP*aPA), (kP - aPA*kA)/(1-aAP*aPA) )
        ## LHe <- c( (kA - aAP*kPi)/(1-aAP*aPA), (kPi - aPA*kA)/(1-aAP*aPA) )
        ## HLe <- c( (kAi - aAP*kP)/(1-aAP*aPA), (kP - aPA*kAi)/(1-aAP*aPA) )
        ## HHe <- c( (kAi - aAP*kPi)/(1-aAP*aPA), (kPi - aPA*kAi)/(1-aAP*aPA) )
        ## LL.sc <- 1/aPA < kA/kP & kA/kP < aAP
        ## LH.sc <- 1/aPA < kA/kPi & kA/kPi < aAP
        ## HL.sc <- 1/aPA < kAi/kP & kAi/kP < aAP
        ## HH.sc <- 1/aPA < kAi/kPi & kAi/kPi < aAP
        ## A.payoff.e <- A.payoff.e +
        ##     matrix(c(
        ##     (LLe[1]>0&LLe[2]>0&LL.sc)|(LLe[1]>0&LLe[2]<=0),
        ##     (LHe[1]>0&LHe[2]>0&LH.sc)|(LHe[1]>0&LHe[2]<=0),
        ##     (HLe[1]>0&HLe[2]>0&HL.sc)|(HLe[1]>0&HLe[2]<=0),
        ##     (HHe[1]>0&HHe[2]>0&HH.sc)|(HHe[1]>0&HHe[2]<=0)),
        ##     2, 2)
        ## P.payoff.e <- P.payoff.e +
        ##     matrix(c(
        ##     (LLe[2]>0&LLe[1]>0&LL.sc)|(LLe[2]>0&LLe[1]<=0),
        ##     (LHe[2]>0&LHe[1]>0&LH.sc)|(LHe[2]>0&LHe[1]<=0),
        ##     (HLe[2]>0&HLe[1]>0&HL.sc)|(HLe[2]>0&HLe[1]<=0),
        ##     (HHe[2]>0&HHe[1]>0&HH.sc)|(HHe[2]>0&HHe[1]<=0)),
        ##     2, 2)
        
    }
    return(Nash)
}
