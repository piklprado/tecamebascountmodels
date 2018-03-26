
## Simulacoes antigas para calcular a probabilidade de cada outcome
## Substituido pelos codigos que calulam tabelas de payoff a posteriori
## ver postNash.R

## As duas especies baixam o r
sim.comb2 <- unique(sim.comb[,1])
nrep <- round(5000/length(sim.comb2))
sim.res2 <- matrix(NA, 2, nrep*length(sim.comb2))
sim.res2.a <- matrix(NA, nrep*length(sim.comb2),3)
for(i in 1:length(sim.comb2)){
    intervalo <- (i*nrep-(nrep-1)):(i*nrep)
    tmp <- do.call("sim.compet", args=list(obj=eval(parse(text=sim.comb2[i])), data=experim1, nrep=nrep, dt=dt))
    sim.res2[,intervalo] <- tmp$estocast[length(dt),,]
    sim.res2.a[intervalo,] <- tmp$analitico
    }
## Apenas Arcella baixa o r
sim.comb3 <- sapply(expand.grid(paste("fit",1:7, sep=""),
                        paste("fit.A",1:3, sep="")),
                     as.character)
nrep <- round(5000/nrow(sim.comb3))
sim.res3 <- matrix(NA,2,nrep*nrow(sim.comb3))
sim.res3.a <- matrix(NA,nrep*nrow(sim.comb3),3)
for(i in 1:nrow(sim.comb3)){
    intervalo <- (i*nrep-(nrep-1)):(i*nrep)
    tmp <- do.call("sim.compet3",
                                    args=list(obj1=eval(parse(text=sim.comb3[i,1])),
                                              obj2=eval(parse(text=sim.comb3[i,2])),
                                              data=experim1, nrep=nrep, dt=dt))
    sim.res3[,intervalo] <- tmp$estocast[length(dt),,]
    sim.res3.a[intervalo,] <- tmp$analitico
    }

## Apenas Pixydiculla baixa o r
sim.comb4 <- sapply(expand.grid(paste("fit",1:7, sep=""),
                        paste("fit.P",1:3, sep="")),
                     as.character)

nrep <- round(5000/nrow(sim.comb4))
sim.res4 <- matrix(NA,2,nrep*nrow(sim.comb4))
sim.res4.a <- matrix(NA,nrep*nrow(sim.comb4),3)
for(i in 1:nrow(sim.comb4)){
    intervalo <- (i*nrep-(nrep-1)):(i*nrep)
    tmp <- do.call("sim.compet4",
                                    args=list(obj1=eval(parse(text=sim.comb4[i,1])),
                                              obj2=eval(parse(text=sim.comb4[i,2])),
                                              data=experim1, nrep=nrep, dt=dt))
    sim.res4[,intervalo] <- tmp$estocast[length(dt),,]
    sim.res4.a[intervalo,] <- tmp$analitico
    }
