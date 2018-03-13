library(zoo)
library(xts)
library(caTools)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rjags)
library(R2jags)

## Cria matrizes de experimentos x tempo
sel.ind <- function(data, batch.number = 1, time.limit = 20){
    filter(data, Time < time.limit, batch==batch.number) %>%
        mutate(replicate=rep(letters[1:3],length(unique(Time))))%>%
        select(Time, replicate, counts) %>%
        spread(Time, value=counts) %>%
        select(-replicate) %>%
        as.matrix()
}

## Cria uma lista com as matrizes de experimentos x tempo com as contagens
## para os experimentos com duas especies
sel.exp <- function(data = comp, batch.number = 1, time.limit = 20){
    Arc <- filter(data, Time < time.limit, batch==batch.number) %>%
        mutate(replicate=rep(letters[1:3],length(unique(Time))))%>%
        select(Time, replicate, arc) %>%
        spread(Time, value=arc) %>%
        select(-replicate) %>%
        as.matrix()
    Pyx <- filter(data, Time < time.limit&batch==batch.number) %>%
        mutate(replicate=rep(letters[1:3],length(unique(Time))))%>%
        select(Time, replicate, pyx) %>%
        spread(Time, value=pyx) %>%
        select(-replicate) %>%
        as.matrix()
    list(Arc, Pyx)
}

## Calculates the overlap of two posterior distributions that are in a object of class "rjags"
post.overl <- function(obj1, obj2, par.name){
    d1 <- density(obj1$BUGSoutput$sims.list[[par.name]])
    d2 <- density(obj2$BUGSoutput$sims.list[[par.name]])
    f1 <- approxfun(d1, yleft=0, yright=0)
    f2 <- approxfun(d2, yleft=0, yright=0)
    f3 <- Vectorize(function(x) min(f1(x),f2(x)))
    if(min(d1$x)>max(d2$x)|min(d2$x)>max(d1$x))
        x <- 0
    else
        x <- integrate(f3, min(c(d1$x, d2$x)), max(d1$x,d2$x))$value
    return(x)
    }

## Plots posteriors of a given parameter for each experiment
post.plot <- function(lista, par.name, transf=FALSE, legend=TRUE, ...){
    if(transf)
        lista <- lapply(lista, function(x) lapply(x,log))
    lista.dens <- lapply(lista, function(x)density(x[[par.name]]))
    x.minimos <- sapply(lista.dens, function(x)min(x$x))
    x.maximos <- sapply(lista.dens, function(x)max(x$x))
    x.limites <- c(min(x.minimos), max(x.maximos))
    y.minimos <- sapply(lista.dens, function(x)min(x$y))
    y.maximos <- sapply(lista.dens, function(x)max(x$y))
    y.limites <- c(min(y.minimos), max(y.maximos))
    plot(lista.dens[[1]], xlim = x.limites, ylim = y.limites, main=par.name, ...)
    if(length(lista.dens)>1){
        for(i in 2:length(lista.dens))
            lines(lista.dens[[i]], col=i)
        if(legend)
            legend("topright", names(lista.dens), lty=1,
                   col=1:length(lista), bty="n")
    }
}

## Gets the mean and 95% CI of the posterior distributions of the parameters that do not vary in Time
dem.pars <- function(object){
    tmp <- object$BUGSoutput$summary
    df <- data.frame(rA=tmp[grep("rA",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")],
                     rP=tmp[grep("rP",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")],
                     kA=tmp[grep("kA",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")],
                     kP=tmp[grep("kP",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")],
                     aPA=tmp[grep("aPA",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")],
                     aAP=tmp[grep("aAP",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")],
                     pA=tmp[grep("^pA",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")],
                     pP=tmp[grep("^pP",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")])
    t(df)
    }
    
## simulates a logistic dynamic as modelled 
logist1 <- function(r, k, p, dt, sampArea, totArea=0.25*90, N0=50){
    tot <- N0
    n <- c()
    for(t in 1:(length(dt))){
        mu <- tot + (r * tot * (1 - tot/k))*dt[t]
        tot <- rpois(1, mu)
        lamb <- tot*sampArea/totArea
        N <- rpois(1,lamb)
        n[t] <- rbinom(1, N, p)
    }
n    
}

## sample from posterior of simulated logistic
## from estimated lambda
post.logist <- function(obj, nrep=1000){
    results <- matrix(ncol=ncol(obj$BUGSoutput$sims.list$lamb), nrow=nrep)
    for(j in 1:nrep){
        i <- sample(1:length(obj$BUGSoutput$sims.list$r), 1)
        lam <- obj$BUGSoutput$sims.list$lamb[i,]
        p <- obj$BUGSoutput$sims.list$p[i]
        N <- rpois(length(lam), lam)
        results[j,] <- rbinom(length(lam), N, p)
    }
    results
}

## Simulates a logistic from their parameters
sim.logist <- function(obj, data, sampArea, totArea=0.25*90, N0=50, nrep=1000){
    dt <- diff(c(0,as.numeric(colnames(data))))
    results <- matrix(ncol=length(dt), nrow=nrep)
    for(j in 1:nrep){
        i <- sample(1:length(obj$BUGSoutput$sims.list$r), 1)
        r <- obj$BUGSoutput$sims.list$r[i]
        k <- obj$BUGSoutput$sims.list$k[i]
        p <- obj$BUGSoutput$sims.list$p[i]
        results[j,] <- logist1(r, k, p, dt, sampArea, totArea, N0)
    }
    results
}

## Simulates the competition dynamics
compet1 <- function(rA, kA, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea=0.25*90, A0=50, P0=50){
    A <- A0
    P <- P0
    yA <- yP <- c()
    for(t in 1:(length(dt))){
        muA <- A + (rA * A * (1 - ((A+aPA*P)/kA)))*dt[t]
        muP <- P + (rP * P * (1 - ((P+aAP*A)/kP)))*dt[t]
        A <- rpois(1, muA)
        P <- rpois(1, muP)
        lambA <- A*sampArea/totArea
        lambP <- P*sampArea/totArea
        nA <- rpois(1,lambA)
        nP <- rpois(1,lambP)
        yA[t] <- rbinom(1, nA, pA)
        yP[t] <- rbinom(1, nP, pP)
    }
    cbind(yA,yP)    
}

## Posterior of competition, from expected posterior of expected population sizes
## Run simulations of competition dynamics, with parameters taken form the posterior distribution
post.compet <- function(obj, nrep=1000){
    nobs <- ncol(obj$BUGSoutput$sims.list$lambA)
    results <- array( dim=c(nobs, 2, nrep))
    for(j in 1:nrep){
        i <- sample(1:nrow(obj$BUGSoutput$sims.list$lambA), 1)
        lambA <- obj$BUGSoutput$sims.list$lambA[i,]
        lambP <- obj$BUGSoutput$sims.list$lambP[i,]
        pA <- obj$BUGSoutput$sims.list$pA[i]
        pP <- obj$BUGSoutput$sims.list$pP[i]
        nA <- rpois(nobs,lambA)
        nP <- rpois(nobs,lambP)
        results[,1,j] <- rbinom(nobs, nA, pA)
        results[,2,j] <- rbinom(nobs, nP, pP)
    }
    results
}

## Run simulations of competition dynamics, with parameters taken from the posterior distribution
sim.compet <- function(obj, data, dt, sampArea=10*pi*0.25^2, totArea=0.25*90, A0=50, P0=50, nrep=1000){
    if(missing(dt))
        dt <- diff(c(0,as.numeric(colnames(data[[1]]))))
    results <- array( dim=c(length(dt), 2, nrep))
    for(j in 1:nrep){
        i <- sample(1:length(obj$BUGSoutput$sims.list$rA), 1)
        rA <- obj$BUGSoutput$sims.list$rA[i]
        kA <- obj$BUGSoutput$sims.list$kA[i]
        pA <- obj$BUGSoutput$sims.list$pA[i]
        rP <- obj$BUGSoutput$sims.list$rP[i]
        kP <- obj$BUGSoutput$sims.list$kP[i]
        pP <- obj$BUGSoutput$sims.list$pP[i]
        aPA <- obj$BUGSoutput$sims.list$aPA[i]
        aAP <- obj$BUGSoutput$sims.list$aAP[i]    
        results[,,j] <- compet1(rA, kA, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea, A0, P0)
    }
    results
}

## Simulates the competition dynamics modelled
## compet1 <- function(A0, P0, rA, rP, kA, kP, aPA, aAP, dt, nIntervals, nSamples, sampArea, totArea){
##     Arc <- Pyx <- lambA <- lambP <- muA <- muP <- c()
##     Arc[1] <- A0
##     Pyx[1] <- P0
##     lambA[1] <- sampArea*Arc[1]/totArea
##     lambP[1] <- sampArea*Pyx[1]/totArea
## 	for(t in 2:nIntervals) {
## 	      ## muA[t-1] <- exp(rA*dt[t-1]) * Arc[t-1] * ( 1 - ((Arc[t-1] + aPA*Pyx[t-1])/kA) )
## 	      ## muP[t-1] <- exp(rP*dt[t-1]) * Pyx[t-1] * ( 1 - ((Pyx[t-1] + aAP*Arc[t-1])/kP) )
##               muA[t-1] <- Arc[t-1] + (rA* Arc[t-1] * ( 1 - ((Arc[t-1] + aPA*Pyx[t-1])/kA)))*dt[t-1]
## 	      muP[t-1] <- Pyx[t-1] + (rP * Pyx[t-1] * ( 1 - ((Pyx[t-1] + aAP*Arc[t-1])/kP) ))*dt[t-1]
## 	      Arc[t] <- rpois(1, muA[t-1])
## 	      Pyx[t] <- rpois(1, muP[t-1] )
## 	      lambA[t] <- Arc[t] * sampArea/totArea
## 	      lambP[t] <- Pyx[t] * sampArea/totArea	
## 	}
##     data.frame(Arc, Pyx, lambA, lambP)
##     }

## zoo object with predicted and observed values
## For logistic
summary.df <- function(fit.obj, data){
    tmp <- fit.obj$BUGSoutput$summary
    tmp2 <- post.logist(fit.obj)
    zoo(
        data.frame(
            obs1=data[1,],
            obs2=data[2,],
            obs3=data[3,],
            obsm = apply(data,2,mean),
        esp=tmp[grep("lamb",rownames(tmp)),"mean"][-1],
        esp.low=tmp[grep("lamb",rownames(tmp)),"2.5%"][-1],
        esp.up=tmp[grep("lamb",rownames(tmp)),"97.5%"][-1],
        post.low=apply(tmp2,2,quantile,0.025)[-1],
        post.up=apply(tmp2,2,quantile,0.975)[-1]),
        order.by=as.numeric(colnames(data)))
}

## for competition
summary.df2 <- function(fit.obj, experim, k=1){
    tmp <- fit.obj$BUGSoutput$summary
    tmp2 <- post.compet(fit.obj)
    zoo(
        data.frame(obsA=runmean(apply(experim[[1]],2,mean),k),
                   obsP=runmean(apply(experim[[2]],2,mean),k),
                   espA=tmp[grep("lambA",rownames(tmp)),"mean"][-1],
                   espA.low=tmp[grep("lambA",rownames(tmp)),"2.5%"][-1],
                   espA.up=tmp[grep("lambA",rownames(tmp)),"97.5%"][-1],
                   espP=tmp[grep("lambP",rownames(tmp)),"mean"][-1],
                   espP.low=tmp[grep("lambP",rownames(tmp)),"2.5%"][-1],
                   espP.up=tmp[grep("lambP",rownames(tmp)),"97.5%"][-1],
                   projA.low=apply(tmp2,c(1,2),quantile,0.025, type=9)[,1][-1],
                   projA.up=apply(tmp2,c(1,2),quantile,0.975, type=9)[,1][-1],
                   projP.low=apply(tmp2,c(1,2),quantile,0.025, type=9)[,2][-1],
                   projP.up=apply(tmp2,c(1,2),quantile,0.975, type=9)[,2][-1]),
        order.by=as.numeric(colnames(experim[[1]])))
}

## Only projected trajectories for competition
## for competition
comp.proj <- function(fit.obj, experim, dt, ...){
    tmp <- sim.compet(fit.obj, dt = dt, data = experim, ...)
    zoo(
        data.frame(
            espA=apply(tmp,c(1,2),mean)[,1][-1],
            espP=apply(tmp,c(1,2),mean)[,2][-1],
            espA.low=apply(tmp,c(1,2),quantile,0.025, type=9)[,1][-1],
            espA.up=apply(tmp,c(1,2),quantile,0.975, type=9)[,1][-1],
            espP.low=apply(tmp,c(1,2),quantile,0.025, type=9)[,2][-1],
            espP.up=apply(tmp,c(1,2),quantile,0.975, type=9)[,2][-1]),
        order.by=cumsum(c(0,dt)))
}

## Dada uma media e um parametro sigma da lognormal retorna o parametro mu
mu.ln <- function(x,lsd) log(x)-(lsd^2)/2

## Plots
## For logistics
p1 <- function(zoo.obj){
    zoo.obj %>%
        fortify() %>%
        ggplot(aes(Index, obsm)) +
        geom_point(size=2) +
        geom_ribbon(aes(ymin=post.low, ymax=post.up), alpha=0.25, fill="blue") +
        scale_y_continuous(name="Number of cells")+
        scale_x_continuous(name="Time (days)")+
        theme_bw()
    }

## for competition
## Predicted total abundances
p2 <- function(zoo.obj){
    fortify(zoo.obj) %>%
        ggplot(aes(Index, espA)) +
        geom_line(colour="red", lwd=1) +
        geom_ribbon(aes(ymin=espA.low, ymax=espA.up), alpha=0.25, fill="red") +
        geom_line(aes(Index,espP),colour="blue", lwd=1) +
        geom_ribbon(aes(ymin=espP.low, ymax=espP.up), alpha=0.25, fill="blue") +
        scale_y_continuous(name="Expected density of cells (1/cm2)")+
        scale_x_continuous(name="Time (days)")+
        theme_bw()
}

## Predicted counts and credible intervals (including detection)
p3 <- function(zoo.obj){
    fortify(zoo.obj) %>%
        ggplot(aes(Index, obsA)) +
        geom_point(colour="red", size=2) +
        geom_ribbon(aes(ymin=projA.low, ymax=projA.up), alpha=0.25, fill="red") +
        geom_point(aes(Index,obsP),colour="blue", size=2) +
        geom_ribbon(aes(ymin=projP.low, ymax=projP.up), alpha=0.25, fill="blue") +
        scale_y_continuous(name="Number of cells")+
        scale_x_continuous(name="Time (days)")+
        theme_bw()
    }