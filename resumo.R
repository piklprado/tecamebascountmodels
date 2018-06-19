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
pA1 <- p1(summary.df(fit.A1, data=arc1), col="red") + ggtitle("Arc #1")
pA2 <- p1(summary.df(fit.A2, data=arc2), col="red")+ ggtitle("Arc #2")
pA3 <- p1(summary.df(fit.A3, data=arc3), col="red")+ ggtitle("Arc #3")
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

## ----compet simul HH-----------------------------------------------------
##  Simulacao da competicao com ambas as especies com valores de r altos
## (tomados das posteriores das culturas monoespecíficas)
sim.comb <- sapply(expand.grid(paste("fit",c(1,2,4,6,7), sep=""),
                        paste("fit.A",1:3, sep=""),
                        paste("fit.P",1:3, sep="")),
                     as.character)

nrep <- round(5000/nrow(sim.comb))
sim.res.e <- matrix(NA, 2, nrep*nrow(sim.comb)) ## para guardar resultados da simulacao estocastica
sim.res.e.s <- matrix(NA, 2, nrep*nrow(sim.comb)) ## para guardar resultados da simulacao sem estocasticidade
sim.res.e.a <- matrix(NA, nrep*nrow(sim.comb), 3) ## para guardar resultados dos valores no equilibrio analitico
for(i in 1:nrow(sim.comb)){
    intervalo <- (i*nrep-(nrep-1)):(i*nrep)
    tmp <- do.call("sim.compet2",
                   args=list(obj1=eval(parse(text=sim.comb[i,1])), obj2=eval(parse(text=sim.comb[i,2])),
                             obj3=eval(parse(text=sim.comb[i,3])),
                             data=experim1, nrep=nrep, dt=dt))
    sim.res.e[,intervalo] <- tmp$estocast[length(dt),,]
    sim.res.e.s[,intervalo] <- tmp$sem.estocast[length(dt),,]
    sim.res.e.a[intervalo, ] <- tmp$analitico
    }

## Repetindo para 30 dias
sim.res <- matrix(NA, 2, nrep*nrow(sim.comb)) ## para guardar resultados da simulacao estocastica
sim.res.s <- matrix(NA, 2, nrep*nrow(sim.comb)) ## para guardar resultados da simulacao sem estocasticidade
sim.res.a <- matrix(NA, nrep*nrow(sim.comb), 3) ## para guardar resultados dos valores no equilibrio analitico
dt <- rep(1/3, 90)
for(i in 1:nrow(sim.comb)){
    intervalo <- (i*nrep-(nrep-1)):(i*nrep)
    tmp <- do.call("sim.compet2",
                   args=list(obj1=eval(parse(text=sim.comb[i,1])), obj2=eval(parse(text=sim.comb[i,2])),
                             obj3=eval(parse(text=sim.comb[i,3])),
                             data=experim1, nrep=nrep, dt=dt))
    sim.res[,intervalo] <- tmp$estocast[length(dt),,]
    sim.res.s[,intervalo] <- tmp$sem.estocast[length(dt),,]
    sim.res.a[intervalo, ] <- tmp$analitico
    }

## ----compet simul LL-----------------------------------------------------
##  Simulacao da competicao com ambas as especies com valores de r baixos
## (tomados das posteriores das culturas monoespecíficas)
sim.comb.ll <- paste("fit",c(1,2,4,6,7), sep="")
nrep.ll <- round(5000/length(sim.comb.ll))
sim.res.e.ll <- matrix(NA, 2, nrep.ll*length(sim.comb.ll))
sim.res.a.ll <- matrix(NA, nrep.ll*length(sim.comb.ll), 3) ## para guardar resultados da solucao analitica no equilibrio
sim.res.e.s.ll <- matrix(NA, 2, nrep.ll*length(sim.comb.ll)) ## simulacao analitica no intervalo
for(i in 1:length(sim.comb.ll)){
    intervalo <- (i*nrep.ll-(nrep.ll-1)):(i*nrep.ll)
    tmp <- do.call("sim.compet",
                   args=list(obj=eval(parse(text=sim.comb.ll[i])), 
                             data=experim1, nrep=nrep.ll))
    sim.res.e.ll[,intervalo] <- tmp$estocast[ncol(experim1[[1]]),,]
    sim.res.e.s.ll[,intervalo] <- tmp$sem.estocast[ncol(experim1[[1]]),,]
    sim.res.a.ll[intervalo, ] <- tmp$analitico
    }
## Para 30 dias
sim.res.30.e.ll <- matrix(NA, 2, nrep.ll*length(sim.comb.ll))
sim.res.30.a.ll <- matrix(NA, nrep.ll*length(sim.comb.ll), 3) ## para guardar resultados da solucao analitica no equilibrio
sim.res.30.e.s.ll <- matrix(NA, 2, nrep.ll*length(sim.comb.ll)) ## simulacao analitica no intervalo
dt <- rep(1/3, 90)
for(i in 1:length(sim.comb.ll)){
    intervalo <- (i*nrep.ll-(nrep.ll-1)):(i*nrep.ll)
    tmp <- do.call("sim.compet",
                   args=list(obj=eval(parse(text=sim.comb.ll[i])), 
                             nrep=nrep.ll, dt = dt))
    sim.res.30.e.ll[,intervalo] <- tmp$estocast[ncol(experim1[[1]]),,]
    sim.res.30.e.s.ll[,intervalo] <- tmp$sem.estocast[ncol(experim1[[1]]),,]
    sim.res.30.a.ll[intervalo, ] <- tmp$analitico
    }

## ----tabela extincoes----------------------------------------------------
t1 <- table(sim.res.e[1,]>0, sim.res.e[2,]>0)/(nrep*nrow(sim.comb))
colnames(t1) <- c("N(Pyx) = 0", "N(Pyx) > 0")
rownames(t1) <- c("N(Arc) = 0", "N(Arc) > 0")
kable(t1, digits=3,
      caption="Proporção de simulações do modelo de competição estocástico com valores de r altos em que a população de cada espécie chegou ao fim de 11.3 dias com tamanhos maiores que zero.")

## ----tabela extincoes 30 dias--------------------------------------------

t1m <- table(sim.res[1,]>0, sim.res[2,]>0)/(nrep*nrow(sim.comb))
colnames(t1m) <- c("N(Pyx) = 0", "N(Pyx) > 0")
rownames(t1m) <- c("N(Arc) = 0", "N(Arc) > 0")
kable(t1m, digits=3,
      caption="Proporção de simulações do modelo de competição em a população de cada espécie chegou ao fim de 30 dias com tamanhos maiores que zero.")

## ----solucao analitica---------------------------------------------------
t1a <- table(sim.res.a[,1]>0, sim.res.a[,2]>0)/(nrep*nrow(sim.comb))
nce <- sum(sim.res.a[,1]>0 & sim.res.a[,2]>0 & !sim.res.a[,3])
t1a[2,2] <- t1a[1,1]-nce/(nrep*nrow(sim.comb))
t1a[1,2] <- t1a[1,2]+nce/(nrep*nrow(sim.comb)*2)
t1a[2,1] <- t1a[2,1]+nce/(nrep*nrow(sim.comb)*2)
colnames(t1a) <- c("N(Pyx) = 0", "N(Pyx) > 0")
rownames(t1a) <- c("N(Arc) = 0", "N(Arc) > 0")
kable(t1a, digits=3,
      caption="Proporção de cálculos do equilíbrio determinístico com $r$ altos em que a população de cada espécie tem tamanhos maiores que zero.")

## ----tabela extincoes sem estocast 11,3 dias-----------------------------
t1m <- table(factor(sim.res.e.s[1,]>0, levels=c("FALSE","TRUE")),
             factor(sim.res.e.s[2,]>0 , levels=c("FALSE", "TRUE")))/(nrep*nrow(sim.comb))
colnames(t1m) <- c("N(Pyx) = 0", "N(Pyx) > 0")
rownames(t1m) <- c("N(Arc) = 0", "N(Arc) > 0")
kable(t1m, digits=3,
      caption="Proporção de simulações do modelo de competição com valores de $r$ alto e sem estocasticidade em que a população de cada espécie chegou ao fim do experimento (11,3 dias) com tamanhos maiores que zero." )

## ----tabela extincoes LL-------------------------------------------------
t1 <- table(sim.res.e.ll[1,]>0, sim.res.e.ll[2,]>0)/(nrep.ll*length(sim.comb.ll))
colnames(t1) <- c("N(Pyx) = 0", "N(Pyx) > 0")
rownames(t1) <- c("N(Arc) = 0", "N(Arc) > 0")
kable(t1, digits=3,
      caption="Proporção de simulações do modelo de competição estocástico com valores de r baixos em que a população de cada espécie chegou ao fim de 11.3 dias com tamanhos maiores que zero.")

## ----tabela extincoes LL 30 dias-----------------------------------------
t1 <- table(sim.res.30.e.ll[1,]>0, sim.res.30.e.ll[2,]>0)/(nrep.ll*length(sim.comb.ll))
colnames(t1) <- c("N(Pyx) = 0", "N(Pyx) > 0")
rownames(t1) <- c("N(Arc) = 0", "N(Arc) > 0")
kable(t1, digits=3,
      caption="Proporção de simulações do modelo de competição estocástico com valores de r baixos em que a população de cada espécie chegou ao fim de 30 dias com tamanhos maiores que zero.")

## ----solucao analitica LL------------------------------------------------
t1b <- table(sim.res.a.ll[,1]>0, sim.res.a.ll[,2]>0)/(nrep.ll*length(sim.comb.ll))
nce <- sum(sim.res.a.ll[,1]>0 & sim.res.a.ll[,2]>0 & !sim.res.a.ll[,3])
t1b[2,2] <- t1b[2,2]-nce/(nrep.ll*length(sim.comb.ll))
t1b[1,2] <- t1b[1,2]+nce/(nrep.ll*length(sim.comb.ll)*2)
t1b[2,1] <- t1b[2,1]+nce/(nrep.ll*length(sim.comb.ll)*2)
colnames(t1b) <- c("N(Pyx) = 0", "N(Pyx) > 0")
rownames(t1b) <- c("N(Arc) = 0", "N(Arc) > 0")
kable(t1b, digits=3,
      caption="Proporção de cálculos do equilíbrio determinístico com $r$ baixos em que a população de cada espécie tem tamanhos maiores que zero.")

## ----tabela extincoes sem estocast no intervalo do experimento-----------
t1m <- table(factor(sim.res.e.s.ll[1,]>0, levels=c("FALSE","TRUE")),
             factor(sim.res.e.s.ll[2,]>0 , levels=c("FALSE", "TRUE")))/(nrep.ll*length(sim.comb.ll))
colnames(t1m) <- c("N(Pyx) = 0", "N(Pyx) > 0")
rownames(t1m) <- c("N(Arc) = 0", "N(Arc) > 0")
kable(t1m, digits=3,
      caption="Proporção de simulações do modelo de competição sem estocasticidade com valores de r baixos em que a população de cada espécie chegou ao fim de 11.3 dias com tamanhos maiores que zero.")

## ----posterior payoffs, eval=FALSE---------------------------------------
## ## Codigo para rodar as matrizes de payoff a posteriori
## ## Um pouco lento, está paralelizado
## source("postNash.R")

## ----prob strategies pop persist-----------------------------------------
f1 <- function(x, i, j, k){
    apply(x[[i]][,j:k], 2, mean)
}
pers.m1 <- matrix(
    c(apply(sapply(post.Nash.short, f1, i=2, j=1, k=12), 1, mean)[1:4],
      apply(results3, 2, mean)),
    byrow=TRUE, ncol=4,
    dimnames=list(c("Equilíbrio de Nash", "Equilíbrio geral"), c("LL", "LH", "HL", "HH"))
)
kable(pers.m1, digits=3, caption="Posterior probabilities of Nash equilibrium for each combination of strategies. Probabilities calculated from simulations of the competition dynamics ran till 11.3 days. Payoff: probability of population persistence at the end of simulations.")

## ----prob strategies pop persist 30 days---------------------------------
f1 <- function(x, i, j, k){
    apply(x[[i]][,j:k], 2, mean)
}
pers.m2 <- matrix(
    c(apply(sapply(post.Nash, f1, i=2, j=1, k=12), 1, mean)[1:4],
      apply(results, 2, mean)),
    byrow=TRUE, ncol=4,
    dimnames=list(c("Equilíbrio de Nash", "Equilíbrio geral"), c("LL", "LH", "HL", "HH"))
)
kable(pers.m2, digits=3, caption="Posterior probabilities of Nash equilibrium for each combination of strategies. Probabilities calculated from simulations of the competition dynamics ran till 30 days. Payoff: probability of population persistence at the end of simulations.")

## ----probs persistencias em cada comb de estrateg 11.3 dias--------------
apply(sapply(post.Nash.short, f1, i=2, j=1, k=12), 1, mean)[5:12] %>%
    matrix(, ncol=2, byrow=TRUE,
           dimnames=list(c("LL","LH", "HL", "HH"), c("Arcella", "Pyxidiculla"))) %>%
    kable( digits=3, caption="Posterior probabilities of population persistence for each combination of strategies. Probabilities calculated from simulations of the competition dynamics ran till 11.3 days.")

## ----probs persistencias em cada comb de estrateg 30 dias----------------
apply(sapply(post.Nash, f1, i=2, j=1, k=12), 1, mean)[5:12] %>%
    matrix(, ncol=2, byrow=TRUE,
           dimnames=list(c("LL","LH", "HL", "HH"), c("Arcella", "Pyxidiculla"))) %>%
    kable( digits=3, caption="Posterior probabilities of population persistence for each combination of strategies. Probabilities calculated from simulations of the competition dynamics ran till 30 days.")

