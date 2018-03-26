library(parallel)
load("arc1.RData")
load("arc2.RData")
load("arc3.RData")
load("pyx1.RData")
load("pyx2.RData")
load("pyx3.RData")
load("compfit1.RData")
load("compfit2.RData")
load("compfit3.RData")
load("compfit4.RData")
load("compfit5.RData")
load("compfit6.RData")
load("compfit7.RData")
source("postfunctions.R")
sim.comb <- sapply(expand.grid(paste("fit",1:7, sep=""),
                        paste("fit.A",1:3, sep=""),
                        paste("fit.P",1:3, sep="")),
                   as.character)
dt <- rep(1/3, 90)
cl1 <- makePSOCKcluster(4)
clusterExport(cl1, c("dt", "sim.comb", "sim.Nash", "compet1",
                     paste("fit",1:7, sep=""),
                     paste("fit.A",1:3, sep=""),
                     paste("fit.P",1:3, sep="")))
sim.Nash.par <- function(i, nomes.objs,...){
    dots <- list(...)
    do.call("sim.Nash",
            args=c(
                list(obj1=eval(parse(text=nomes.objs[i,1])),
                     obj2=eval(parse(text=nomes.objs[i,2])),
                     obj3=eval(parse(text=nomes.objs[i,3]))),
                dots)
            )
}

## 60 days
post.Nash <- parLapply(cl1, X=1:nrow(sim.comb), fun=sim.Nash.par,
                       nomes.obj = sim.comb, dt=dt, data=experim1, nsamp=20, nrep=100)
## Only the time of the experiment
post.Nash.short <- parLapply(cl1, X=1:nrow(sim.comb), fun=sim.Nash.par,
                             nomes.obj = sim.comb, data=experim1, nsamp=20, nrep=100)
## 120 days
post.Nash.120 <- parLapply(cl1, X=1:nrow(sim.comb), fun=sim.Nash.par,
                       nomes.obj = sim.comb, dt=rep(1/3, 360), data=experim1, nsamp=20, nrep=100)
stopCluster(cl1)
save.image()
save(post.Nash.120, file="postNash120.RData")
save(post.Nash, file="postNash.RData")
save(post.Nash.short, file="postNashShort.RData")
## Proportion of Nash equilibrium for each possibility
f1 <- function(x, i, j, k){
    apply(x[[i]][,j:k], 2, mean)
}
##  60days
apply(sapply(post.Nash, f1, i=2, j=1, k=4), 1, mean)
apply(sapply(post.Nash, f1, i=1, j=1, k=4), 1, mean)
## Mean values, short term
apply(sapply(post.Nash.short, f1, i=2, j=5, k=12), 1, mean)

## Tentando este pacote, mas não entendi
## library(GPGame)
## eq <- "NE"
## #eq <- "NKSE"
## (integcontrol <- generate_integ_pts(n.s=c(2,2), d=2, nobj=2, equilibrium=eq, gridtype="cartesian"))
## integ.pts <- integcontrol$integ.pts
## (Z <- matrix(post.Nash[[56]][[2]][4,c(5,9,7,11,6,10,8,12)], ncol=2))
## (Z <- matrix(c(1, 2, -1, 0, 1, -1, 2, 0), ncol=2)) # Dilema do prisioneiro
## trueEq <- getEquilibrium(Z = Z, equilibrium = eq, nobj = 2, n.s = c(2,2),
##                                return.design = TRUE, expanded.indices = integcontrol$expanded.indices)
## ## Equilibrium
## trueEq$NE
## print(integ.pts[trueEq$NE,])
## ## Payoffs
## trueEq$NEPoff

## Com função que criei
## Verificando o dilema do prisioneiro
check.eq(1, -1, 2, 0, 1, 2, -1, 0)
## Testando com dados da simulacao
index <- c(5, 7, 9, 11, 6, 8, 10, 12)
j <- 1
i <- 8
post.Nash[[i]][[2]][j,c(1:4, index)]
check.eq(lista=as.list(post.Nash[[i]][[2]][j,5:12]))

## Varrendo as listas
## 1. Simulacoes por 60 dias
## por prob persistencia
results <- matrix(NA, nrow=length(post.Nash)*nrow(post.Nash[[1]][[1]]), ncol=4,
                  dimnames=list(NULL, c("LL","LH","HL","HH")))
k <- 1
for(i in 1:length(post.Nash)){
    for(j in 1:nrow(post.Nash[[1]][[1]])){
        results[k,] <-check.eq(lista = as.list(post.Nash[[i]][[2]][j,5:12]))
        k <- k+1
        }
}
## Por tamanho pop
results2 <- matrix(NA, nrow=length(post.Nash)*nrow(post.Nash[[1]][[1]]), ncol=4,
                   dimnames=list(NULL, c("LL","LH","HL","HH")))
k <- 1
for(i in 1:length(post.Nash)){
    for(j in 1:nrow(post.Nash[[1]][[1]])){
        results2[k,] <-check.eq(lista = as.list(post.Nash[[i]][[1]][j,5:12]))
        k <- k+1
        }
}

## 2. Simulacoes pela duracao do experimento
## por prob persistencia
results3 <- matrix(NA, nrow=length(post.Nash.short)*nrow(post.Nash.short[[1]][[1]]), ncol=4,
                  dimnames=list(NULL, c("LL","LH","HL","HH")))
k <- 1
for(i in 1:length(post.Nash.short)){
    for(j in 1:nrow(post.Nash.short[[1]][[1]])){
        results3[k,] <-check.eq(lista = as.list(post.Nash.short[[i]][[2]][j,5:12]))
        k <- k+1
        }
}
## Por tamanho pop
results4 <- matrix(NA, nrow=length(post.Nash.short)*nrow(post.Nash.short[[1]][[1]]), ncol=4,
                   dimnames=list(NULL, c("LL","LH","HL","HH")))
k <- 1
for(i in 1:length(post.Nash.short)){
    for(j in 1:nrow(post.Nash.short[[1]][[1]])){
        results4[k,] <-check.eq(lista = as.list(post.Nash.short[[i]][[1]][j,5:12]))
        k <- k+1
        }
}

## 3. Simulacoes por 120 dias
## por prob persistencia
results5 <- matrix(NA, nrow=length(post.Nash.120)*nrow(post.Nash.120[[1]][[1]]), ncol=4,
                  dimnames=list(NULL, c("LL","LH","HL","HH")))
k <- 1
for(i in 1:length(post.Nash.120)){
    for(j in 1:nrow(post.Nash.120[[1]][[1]])){
        results5[k,] <-check.eq(lista = as.list(post.Nash.120[[i]][[2]][j,5:12]))
        k <- k+1
        }
}
## Por tamanho pop
results6 <- matrix(NA, nrow=length(post.Nash.120)*nrow(post.Nash.120[[1]][[1]]), ncol=4,
                   dimnames=list(NULL, c("LL","LH","HL","HH")))
k <- 1
for(i in 1:length(post.Nash.120)){
    for(j in 1:nrow(post.Nash.120[[1]][[1]])){
        results6[k,] <-check.eq(lista = as.list(post.Nash.120[[i]][[1]][j,5:12]))
        k <- k+1
        }
}


## Strict Nash equilibrium and multiple equilibria
f1 <- function(x, i, j, k){
    apply(x[[i]][,j:k], 2, mean)
}

## Simulacoes pela duracao do experimento
## Persistence probability
apply(sapply(post.Nash.short, f1, i=2, j=1, k=12), 1, mean)
apply(results3, 2, mean)
## Relative population size
apply(sapply(post.Nash.short, f1, i=1, j=1, k=12), 1, mean)
apply(results4, 2, mean)

## Simulacoes por 60 dias
## Presistence prob
apply(sapply(post.Nash, f1, i=2, j=1, k=12), 1, mean)
apply(results, 2, mean)
## Pop size
apply(sapply(post.Nash, f1, i=1, j=1, k=12), 1, mean)
apply(results2, 2, mean)

## Simulacoes por 120 dias
## Persit prob
apply(sapply(post.Nash.120, f1, i=2, j=1, k=12), 1, mean)
apply(results5, 2, mean, na.rm=TRUE)
## Pop size
apply(sapply(post.Nash.120, f1, i=1, j=1, k=12), 1, mean)
apply(results6, 2, mean)
