library(plyr)
library(dplyr)
library(R2jags)
load("cluster/arc1.RData")
load("cluster/arc2.RData")
load("cluster/arc3.RData")
load("cluster/pyx1.RData")
load("cluster/pyx2.RData")
load("cluster/pyx3.RData")
load("cluster/compfit1.RData")
load("cluster/compfit2.RData")
load("cluster/compfit3.RData")
load("cluster/compfit4.RData")
load("cluster/compfit5.RData")
load("cluster/compfit6.RData")
load("cluster/compfit7.RData")

f1 <- function(obj, index){
    x <- as.data.frame(summary(as.mcmc(obj)[,c("rA","rP","kA","kP","aAP", "aPA")], quantile=c(0.025,0.5,0.975))[[2]])
    names(x) <- c("lwr", "median", "upr")
    x$type <- "both"
    x$experiment <- index
    x$variable <- substr(rownames(x),1,1)
    x$sp <- substr(rownames(x),2,2)
    x
}

f2 <- function(obj, index, sp){
    x <- as.data.frame(summary(as.mcmc(obj)[,c("r","k")], quantile=c(0.025,0.5,0.975))[[2]])
    names(x) <- c("lwr", "median", "upr")
    x$type <- "single"
    x$experiment <- index
    x$variable <- rownames(x)
    x$sp <- sp
    x
    }

tmp1 <- f1(fit1, 1)
lista <- list(fit1,fit2, fit3, fit4, fit5, fit6, fit7)
lista2 <- list(fit.A1,fit.A2, fit.A3, fit.P1,fit.P2, fit.P3)
sps <- rep(c("A","P"), each=3)
exps <- rep(1:3, each=2)
for(i in 2:7)
    tmp1 <- rbind(tmp1, f1(lista[[i]], i))
for(i in 1:6)
    tmp1 <- rbind(tmp1, f2(lista2[[i]], exps[i], sps[i]))

## The figure with parameters estimates and credibility intervals for all experiments
##pdf("fig5_alt.pdf", width = 10.5)
pdf("fig5_alt.pdf")

##par(mar=c(5,6,4,0), mfrow=c(1,3), las=1)
par(fig=c(0.1,0.4,0,1), mar=c(5,1,4,1), new=TRUE, las=1)
## r
espaco <- 0.1
espaco2 <- espaco*3
################################################################################
## r
## Open frame
plot(1:10, 1:10,
     xlim = range(c(tmp1$upr[tmp1$variable=="r"], tmp1$lwr[tmp1$variable=="r"])),
     type="n", axes=FALSE, xlab="Intrinsic growth ratio (r)", ylab="")
## Competition experiments, Arcella
with(subset(tmp1, sp=="A"&variable=="r"&type=="both"),{
    points(median, (1:7)-espaco, col="red", pch=19)
    segments(lwr, (1:7)- espaco, upr, (1:7)-espaco, col="red")
}
)
## Competition experiments, Pyx
with(subset(tmp1, sp=="P"&variable=="r"&type=="both"),{
    points(median, (1:7)+espaco, col="blue", pch=19)
    segments(lwr, (1:7) + espaco, upr, (1:7)+espaco, col="blue")
}
)
## Single experiments, Arc
with(subset(tmp1, sp=="A"&variable=="r"&type=="single"),{
    points(median, 7.6 +(1:3)*espaco2, col="red", pch=19)
    segments(lwr, 7.6 +(1:3)*espaco2, upr, 7.6 +(1:3)*espaco2, col="red")
}
)
## Single experim, Pyx
with(subset(tmp1, sp=="P"&variable=="r"&type=="single"),{
    points(median, 9.1 +(1:3)*espaco2, col="blue", pch=19.1)
    segments(lwr, 9.1 +(1:3)*espaco2, upr, 9.1+(1:3)*espaco2, col="blue")
}
)
## Axes, scales ...
axis(2, at=c(1:7), labels=paste("#",1:7), tick=FALSE)
axis(1)
box()
abline(h=(1:7)+0.5,
       lty=rep(c(2,1),c(6,1)),
       col=rep(c("gray","black"), c(6,1)))
par(las=0)
mtext("Competition", at=4, side=2, line=2.75, cex=1.5)
mtext("Mono-specific", at=9, side=2, line=2.75, cex=1.5)
################################################################################
## K
## Open frame
#par(mar=c(5,2,4,2), las=1)
par(fig=c(0.4,0.7,0,1), new=TRUE, las=1)
plot(1:10, 1:10,
     xlim = range(c(tmp1$upr[tmp1$variable=="k"], tmp1$lwr[tmp1$variable=="k"])),
     type="n", axes=FALSE, xlab="Carrying capacity (K)", ylab="", log="x")
## Competition experiments, Arcella
with(subset(tmp1, sp=="A"&variable=="k"&type=="both"),{
    points(median, (1:7)-espaco, col="red", pch=19)
    segments(lwr, (1:7)- espaco, upr, (1:7)-espaco, col="red")
}
)
## Competition experiments, Pyx
with(subset(tmp1, sp=="P"&variable=="k"&type=="both"),{
    points(median, (1:7)+espaco, col="blue", pch=19)
    segments(lwr, (1:7) + espaco, upr, (1:7)+espaco, col="blue")
}
)
## Single experiments, Arc
with(subset(tmp1, sp=="A"&variable=="k"&type=="single"),{
    points(median, 7.6 +(1:3)*espaco2, col="red", pch=19)
    segments(lwr, 7.6 +(1:3)*espaco2, upr, 7.6 +(1:3)*espaco2, col="red")
}
)
## Single experim, Pyx
with(subset(tmp1, sp=="P"&variable=="k"&type=="single"),{
    points(median, 9.1 +(1:3)*espaco2, col="blue", pch=19.1)
    segments(lwr, 9.1 +(1:3)*espaco2, upr, 9.1+(1:3)*espaco2, col="blue")
}
)
## Axes, scales ...
axis(1)
box()
abline(h=(1:7)+0.5,
       lty=rep(c(2,1),c(6,1)),
       col=rep(c("gray","black"), c(6,1)))
################################################################################
## coefficients
## Open frame
par(fig=c(0.7,1,0,0.78), las=1, new=TRUE)
box()
par(fig=c(0.7,1,0,1), las=1, new=TRUE)
plot(1:10, 1:10,
     xlim = range(c(tmp1$upr[tmp1$variable=="a"], tmp1$lwr[tmp1$variable=="a"])),
     type="n", axes=FALSE, xlab="Competition coefficients", ylab="", log="x")
## Competition experiments, Arcella
with(subset(tmp1, sp=="A"&variable=="a"&type=="both"),{
    points(median, (1:7)-espaco, col="red", pch=19)
    segments(lwr, (1:7)- espaco, upr, (1:7)-espaco, col="red")
}
)
## Competition experiments, Pyx
with(subset(tmp1, sp=="P"&variable=="a"&type=="both"),{
    points(median, (1:7)+espaco, col="blue", pch=19)
    segments(lwr, (1:7) + espaco, upr, (1:7)+espaco, col="blue")
}
)
## Legend
legend(x = 0.05, y = 10, c("A. intermedia", "P. operculata"),
       pch=19, lty=1, col=c("red", "blue"), bty="n")
## Axes, scales ...
axis(1)
abline(h=(1:6)+0.5,
       lty=2,
       col="gray")

dev.off()
################################################################################
