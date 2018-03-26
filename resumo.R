## ----setup, echo=FALSE, warning=FALSE, message=F-------------------------
library(knitr); library(zoo); library(xts);
library(caTools);
library(ggplot2)#; library(cowplot)
library(gridExtra)
library(dplyr); library(tidyr);
library(rjags); library(R2jags);
library(mcmcplots); library(superdiag)
opts_chunk$set(fig.align = 'center', fig.show = 'hold',
               fig.height = 4, warning = FALSE, message = FALSE,
               error = FALSE, echo=FALSE)
options(formatR.arrow = TRUE,width = 90, cache=TRUE, scipen=999)
source("functions.R") ## auxiliary functions 

## ----load model fits, eval=FALSE, echo=TRUE------------------------------
## load("cluster/arc1.RData")
## load("cluster/arc2.RData")
## load("cluster/arc3.RData")
## load("cluster/pyx1.RData")
## load("cluster/pyx2.RData")
## load("cluster/pyx3.RData")
## load("cluster/compfit1.RData")
## load("cluster/compfit2.RData")
## load("cluster/compfit3.RData")
## load("cluster/compfit4.RData")
## load("cluster/compfit5.RData")
## load("cluster/compfit6.RData")
## load("cluster/compfit7.RData")

## ----load data single-species experiments--------------------------------
## Tidy data table with all experiments
pyx.ind <- read.csv2("../dados_giulia/pyx_ind.csv") %>%
    gather(Pyx1A:Pyx3C, key=id, value=counts) %>%
    mutate(batch=as.factor(substr(id,4,4)))
arc.ind <- read.csv2("../dados_giulia/arc_ind.csv") %>%
    gather(Arcella1A:Arcella3C, key=id, value=counts) %>%
    mutate(batch=as.factor(substr(id,8,8)))
## Separate matrices for each experiment
pyx1 <- sel.ind(pyx.ind, 1)
pyx2 <- sel.ind(pyx.ind, 2)
pyx3 <- sel.ind(pyx.ind, 3)
arc1 <- sel.ind(arc.ind, 1)
arc2 <- sel.ind(arc.ind, 2)
arc3 <- sel.ind(arc.ind, 3)

## ----single species obs x pred, fig.height=12----------------------------
pP1 <- p1(summary.df(fit.P1, data=pyx1)) + ggtitle("Pyx #1")
pP2 <- p1(summary.df(fit.P2, data=pyx2))+ ggtitle("Pyx #2")
pP3 <- p1(summary.df(fit.P3, data=pyx3))+ ggtitle("Pyx #3")
## Arcella
pA1 <- p1(summary.df(fit.A1, data=arc1)) + ggtitle("Arc #1")
pA2 <- p1(summary.df(fit.A2, data=arc2))+ ggtitle("Arc #2")
pA3 <- p1(summary.df(fit.A3, data=arc3))+ ggtitle("Arc #3")
grid.arrange(pP1,  pA1, pP2, pA2, pP3, pA3)

## ----logistic posteriors, fig.height=13.5--------------------------------
lista.single <- list(
    Arc1=fit.A1$BUGSoutput$sims.list,
    Arc2=fit.A2$BUGSoutput$sims.list,
    Arc3=fit.A3$BUGSoutput$sims.list,
    Pyx1=fit.P1$BUGSoutput$sims.list,
    Pyx2=fit.P2$BUGSoutput$sims.list,
    Pyx3=fit.P3$BUGSoutput$sims.list)
par(mfrow=c(3,1), lwd=1.5)
post.plot(lista.single,par.name="p",  legend=FALSE, xlab="Parameter value")
post.plot(lista.single,par.name="k", legend=FALSE, xlab="Parameter value")
post.plot(lista.single,par.name="r", xlab="Parameter value")
par(mfrow=c(1,1))

## ----competicao leitura dados--------------------------------------------
comp <- read.csv2("../dados_giulia/comp_counts_sintese_final.csv")
experim1 <- sel.exp(comp, batch = 1)
experim2 <- sel.exp(comp, batch = 2)
experim3 <- sel.exp(comp, batch = 3)
experim4 <- sel.exp(comp, batch = 4)
experim5 <- sel.exp(comp, batch = 5)
experim6 <- sel.exp(comp, batch = 6)
experim7 <- sel.exp(comp, batch = 7)

## ----competition obs x predicted counts, fig.height=13.5-----------------
## Graficos
pe1 <- p3(summary.df2(fit1, experim1)) + ggtitle("Experiment #1") +
    annotate("text", x = c(8,9.5), y = c(70,7.5),
             label = c("Pyxidicula","Arcella"),
             color=c("blue","red"), size = 4)
pe2 <- p3(summary.df2(fit2, experim2)) + ggtitle("Experiment #2")
pe3 <- p3(summary.df2(fit3, experim3)) + ggtitle("Experiment #3")
pe4 <- p3(summary.df2(fit4, experim4)) + ggtitle("Experiment #4")
pe5 <- p3(summary.df2(fit5, experim5)) + ggtitle("Experiment #5")
pe6 <- p3(summary.df2(fit6, experim6)) + ggtitle("Experiment #6")
pe7 <- p3(summary.df2(fit7, experim7)) + ggtitle("Experiment #7")
grid.arrange(pe1,pe2,pe3,pe4,pe5,pe6,pe7, ncol=2)

## ----competition posteriors, fig.height=12-------------------------------
lista2 <- list(
    Exp1=fit1$BUGSoutput$sims.list,
    Exp2=fit2$BUGSoutput$sims.list,
    Exp3=fit3$BUGSoutput$sims.list,
    Exp4=fit4$BUGSoutput$sims.list,
    Exp5=fit5$BUGSoutput$sims.list,
    Exp6=fit6$BUGSoutput$sims.list,
    Exp7=fit7$BUGSoutput$sims.list)
par(mfrow=c(3,2), lwd=1.5)
post.plot(lista2,par.name="kA", transf=TRUE, legend=FALSE, xlab="Log(Parameter value)")
post.plot(lista2,par.name="kP", transf=TRUE, legend=FALSE, xlab="Log(Parameter value)")
post.plot(lista2,par.name="rA", legend=FALSE, xlab="Parameter value")
post.plot(lista2,par.name="rP", legend=FALSE, xlab="Parameter value")
post.plot(lista2,par.name="aAP", legend=FALSE, xlab="Parameter value")
post.plot(lista2,par.name="aPA")
par(mfrow=c(1,1))

## ----compet simul--------------------------------------------------------
sim.comb <- sapply(expand.grid(paste("fit",1:7, sep=""),
                        paste("fit.A",1:3, sep=""),
                        paste("fit.P",1:3, sep="")),
                     as.character)

nrep <- round(5000/nrow(sim.comb))
sim.res.e <- matrix(NA, 2, nrep*nrow(sim.comb))
for(i in 1:nrow(sim.comb)){
    intervalo <- (i*nrep-(nrep-1)):(i*nrep)
    sim.res.e[,intervalo] <- do.call("sim.compet2",
                                    args=list(obj1=eval(parse(text=sim.comb[i,1])), obj2=eval(parse(text=sim.comb[i,2])),
                                              obj3=eval(parse(text=sim.comb[i,3])),
                                              data=experim1, nrep=nrep))$estocast[ncol(experim1[[1]]),,]
    }

## Repetindo para 60 dias
nrep <- round(5000/nrow(sim.comb))
sim.res <- matrix(NA, 2, nrep*nrow(sim.comb))
sim.res.a <- matrix(NA, nrep*nrow(sim.comb), 3)
dt <- rep(1/5, 300)
for(i in 1:nrow(sim.comb)){
    intervalo <- (i*nrep-(nrep-1)):(i*nrep)
    tmp <- do.call("sim.compet2",
                   args=list(obj1=eval(parse(text=sim.comb[i,1])), obj2=eval(parse(text=sim.comb[i,2])),
                             obj3=eval(parse(text=sim.comb[i,3])),
                             data=experim1, nrep=nrep, dt=dt))
    sim.res[,intervalo] <- tmp$estocast[length(dt),,]
    sim.res.a[intervalo, ] <- tmp$analitico
    }

## ----tabela extincoes----------------------------------------------------
t1 <- table(sim.res.e[1,]>0, sim.res[2,]>0)/(nrep*nrow(sim.comb))
colnames(t1) <- c("N(Pyx) = 0", "N(Pyx) > 0")
rownames(t1) <- c("N(Arc) = 0", "N(Arc) > 0")
kable(t1, digits=3,
      caption="Proporção de simulações do modelo de competição em a população de cada espécie chegou ao fim de 11.3 dias com tamanhos maiores que zero.")

## ----tabela extincoes 60 dias--------------------------------------------
t1m <- table(sim.res[1,]>0, sim.res[2,]>0)/(nrep*nrow(sim.comb))
colnames(t1m) <- c("N(Pyx) = 0", "N(Pyx) > 0")
rownames(t1m) <- c("N(Arc) = 0", "N(Arc) > 0")
kable(t1m, digits=3,
      caption="Proporção de simulações do modelo de competição em a população de cada espécie chegou ao fim de 60 dias com tamanhos maiores que zero.")

## ----solucao analitica---------------------------------------------------
t1a <- table(sim.res.a[,1]>0, sim.res.a[,2]>0)/(nrep*nrow(sim.comb))
nce <- sum(sim.res.a[,1]>0 & sim.res.a[,2]>0 & !sim.res.a[,3])
t1a[1,1] <- t1a[1,1]-nce/(nrep*nrow(sim.comb))
t1a[1,2] <- t1a[1,2]+nce/(nrep*nrow(sim.comb)*2)
t1a[2,1] <- t1a[2,1]+nce/(nrep*nrow(sim.comb)*2)
colnames(t1a) <- c("N(Pyx) = 0", "N(Pyx) > 0")
rownames(t1a) <- c("N(Arc) = 0", "N(Arc) > 0")
kable(t1a, digits=3,
      caption="Proporção das combinações de parâmetros estimados em que  a população de cada espécie persiste, sob o modelo determinístico de competição.")

## ----compet simul payoffs, eval=FALSE------------------------------------
## ## As duas especies baixam o r
## sim.comb2 <- unique(sim.comb[,1])
## nrep <- round(5000/length(sim.comb2))
## sim.res2 <- matrix(NA, 2, nrep*length(sim.comb2))
## sim.res2.a <- matrix(NA, nrep*length(sim.comb2),3)
## for(i in 1:length(sim.comb2)){
##     intervalo <- (i*nrep-(nrep-1)):(i*nrep)
##     tmp <- do.call("sim.compet", args=list(obj=eval(parse(text=sim.comb2[i])), data=experim1, nrep=nrep, dt=dt))
##     sim.res2[,intervalo] <- tmp$estocast[length(dt),,]
##     sim.res2.a[intervalo,] <- tmp$analitico
##     }
## ## Apenas Arcella baixa o r
## sim.comb3 <- sapply(expand.grid(paste("fit",1:7, sep=""),
##                         paste("fit.A",1:3, sep="")),
##                      as.character)
## nrep <- round(5000/nrow(sim.comb3))
## sim.res3 <- matrix(NA,2,nrep*nrow(sim.comb3))
## sim.res3.a <- matrix(NA,nrep*nrow(sim.comb3),3)
## for(i in 1:nrow(sim.comb3)){
##     intervalo <- (i*nrep-(nrep-1)):(i*nrep)
##     tmp <- do.call("sim.compet3",
##                                     args=list(obj1=eval(parse(text=sim.comb3[i,1])),
##                                               obj2=eval(parse(text=sim.comb3[i,2])),
##                                               data=experim1, nrep=nrep, dt=dt))
##     sim.res3[,intervalo] <- tmp$estocast[length(dt),,]
##     sim.res3.a[intervalo,] <- tmp$analitico
##     }
## 
## ## Apenas Pixydiculla baixa o r
## sim.comb4 <- sapply(expand.grid(paste("fit",1:7, sep=""),
##                         paste("fit.P",1:3, sep="")),
##                      as.character)
## 
## nrep <- round(5000/nrow(sim.comb4))
## sim.res4 <- matrix(NA,2,nrep*nrow(sim.comb4))
## sim.res4.a <- matrix(NA,nrep*nrow(sim.comb4),3)
## for(i in 1:nrow(sim.comb4)){
##     intervalo <- (i*nrep-(nrep-1)):(i*nrep)
##     tmp <- do.call("sim.compet4",
##                                     args=list(obj1=eval(parse(text=sim.comb4[i,1])),
##                                               obj2=eval(parse(text=sim.comb4[i,2])),
##                                               data=experim1, nrep=nrep, dt=dt))
##     sim.res4[,intervalo] <- tmp$estocast[length(dt),,]
##     sim.res4.a[intervalo,] <- tmp$analitico
##     }

## ----tabela payoffs estocast---------------------------------------------
## Arcella
A.payoff <- matrix(
    c(
        sum(sim.res2[1,]>0)/dim(sim.res2)[2], # A baixo, P baixo
        sum(sim.res4[1,]>0)/dim(sim.res4)[2], # A baixo, P alto
        sum(sim.res3[1,]>0)/dim(sim.res3)[2], # A alto, P baixo
        sum(sim.res[1,]>0)/dim(sim.res)[2] # A alto, P alto
        ),
    ncol=2,
    dimnames=list(c("P r low", "P r high"), c("A r low", "A r high"))
)
## Pyxidicolla
P.payoff <- matrix(
    c(
        sum(sim.res2[2,]>0)/dim(sim.res2)[2], # A baixo, P baixo
        sum(sim.res4[2,]>0)/dim(sim.res4)[2], # A baixo, P alto
        sum(sim.res3[2,]>0)/dim(sim.res3)[2], # A alto, P baixo
        sum(sim.res[2,]>0)/dim(sim.res)[2] # A alto, P alto
        ),
    ncol=2,
    dimnames=list(c("P r low", "P r high"), c("A r low", "A r high"))
)

## ----arcella payoff estocast---------------------------------------------
kable(A.payoff, digits=4, caption="Arcella - probabilidades estimadas de persistência pelo modelo com estocasticidade")

## ----Pixydiculla payoff estocast-----------------------------------------
kable(P.payoff, digits=4, caption="Pixydiculla - probabilidades estimadas de persistência pelo modelo com estocasticidade")

## ----tabela payoffs no estocast------------------------------------------
## Arcella
A.payoff2 <- matrix(
    c(
        (sum(sim.res2.a[,1]>0) - sum(sim.res2.a[,1]>0 & sim.res2.a[,2]>0 & !sim.res2.a[,3])/2)/dim(sim.res2.a)[1], # A baixo, P baixo
        (sum(sim.res4.a[,1]>0) - sum(sim.res4.a[,1]>0 & sim.res4.a[,2]>0 & !sim.res4.a[,3])/2)/dim(sim.res4.a)[1], # A baixo, P alto
        (sum(sim.res3.a[,1]>0) - sum(sim.res3.a[,1]>0 & sim.res3.a[,2]>0 & !sim.res3.a[,3])/2)/dim(sim.res3.a)[1], # A alto, P baixo
        (sum(sim.res.a[,1]>0)- sum(sim.res.a[,1]>0 & sim.res.a[,2]>0 & !sim.res.a[,3])/2)/dim(sim.res.a)[1] # A alto, P alto
        ),
    ncol=2,
    dimnames=list(c("P r low", "P r high"), c("A r low", "A r high"))
)
## Pyxidicolla
P.payoff2 <- matrix(
    c(
        (sum(sim.res2.a[,2]>0) - sum(sim.res2.a[,1]>0 & sim.res2.a[,2]>0 & !sim.res2.a[,3])/2)/dim(sim.res2.a)[1], # A baixo, P baixo
        (sum(sim.res4.a[,2]>0) - sum(sim.res4.a[,1]>0 & sim.res4.a[,2]>0 & !sim.res4.a[,3])/2)/dim(sim.res4.a)[1], # A baixo, P alto
        (sum(sim.res3.a[,2]>0) - sum(sim.res3.a[,1]>0 & sim.res3.a[,2]>0 & !sim.res3.a[,3])/2)/dim(sim.res3.a)[1], # A alto, P baixo
        (sum(sim.res.a[,2]>0)- sum(sim.res.a[,1]>0 & sim.res.a[,2]>0 & !sim.res.a[,3])/2)/dim(sim.res.a)[1] # A alto, P alto
        ),
    ncol=2,
    dimnames=list(c("P r low", "P r high"), c("A r low", "A r high"))
)

## ----arcella payoff no estocast------------------------------------------
kable(A.payoff2, digits=4, caption="Arcella - probabilidades estimadas de persistência pelo modelo sem estocasticidade")

## ----Pixydiculla payoff no estocast--------------------------------------
kable(P.payoff2, digits=4, caption="Pixydiculla - probabilidades estimadas de persistência pelo modelo sem estocasticidade")

## ----posterior payoffs, eval=FALSE---------------------------------------
## sim.comb <- sapply(expand.grid(paste("fit",1:7, sep=""),
##                         paste("fit.A",1:3, sep=""),
##                         paste("fit.P",1:3, sep="")),
##                    as.character)
## nrep <- round(5000/nrow(sim.comb))
## A.sim.pay <- P.sim.pay <- A.sim.pay.e <- P.sim.pay.e <- array(dim=c(nrow(sim.comb),2,2))
## dt <- rep(1/5, 300)
## for(i in 1:nrow(sim.comb)){
##     tmp <- do.call("sim.Nash",
##                    args=list(obj1=eval(parse(text=sim.comb[i,1])), obj2=eval(parse(text=sim.comb[i,2])),
##                              obj3=eval(parse(text=sim.comb[i,3])),
##                              data=experim1, nrep=nrep, dt=dt))
##     A.sim.pay[i,,] <- tmp$A.payoff
##     P.sim.pay[i,,] <- tmp$P.payoff
##     A.sim.pay.e[i,,] <- tmp$A.payoff.e
##     P.sim.pay.e[i,,] <- tmp$P.payoff.e
##     }
## apply(A.sim.pay, c(2,3), mean)
## apply(P.sim.pay, c(2,3), mean)
## apply(A.sim.pay.e, c(2,3), mean)
## apply(P.sim.pay.e, c(2,3), mean)

