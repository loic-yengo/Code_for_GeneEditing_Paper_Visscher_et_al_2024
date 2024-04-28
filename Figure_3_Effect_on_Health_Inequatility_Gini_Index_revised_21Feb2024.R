setwd("~/Desktop/Papers/Gene-editing/Resubmission - Code - Figure/")
library(DescTools)

Data <- rbind(
  AD=c(K=0.05,Fold=7.966382),
  MDD=c(K=0.15,Fold=1.743556),
  SCZ=c(K=0.01,Fold=10.703546),
  T2D=c(K=0.10,Fold=60.884997),
  CAD=c(K=0.06,Fold=32.193018)
)
Data <- as.data.frame(Data)

FoldChangeInGini <- function(m,p,K,N){
  # m is the shift in disease liability for edited genomes
  # p is the proportion of edited genomes
  # K is the disease prevalence
  # N is the population size set to 1M (does not influence much the results)

  if(p==0){
    return(1)
  }else{
    U  <- seq(1/N,1-1/N,length=N)
    x  <- qnorm(U)
    t  <- qnorm(1-K)
    R0 <- pnorm(x,mean=-t,sd=1,lower.tail = F) # Risk distribution in the population
    G0 <- Gini(R0, unbiased=FALSE) # Gini index of risk in unedited population
    n  <- round(p * N) # size of edited sub-population
    xE <- rnorm(n=n,mean=m,sd=1) # liability of edited genomes (should it be the same variance?? (maybe 1-q^2))
    RE <- c(R0[-(1:n)],pnorm(xE,mean=-t,sd=1,lower.tail = F)) # Risk distribution in mixture population
    GE <- Gini(RE, unbiased=FALSE)
    return(GE/G0)
  }
}

## Scenarios
Cols <- c("dodgerblue","coral1","goldenrod")
ps   <- seq(0.0,1,by=0.05)   # values for p
Ks   <- c(0.01,0.05,0.1,0.2) # values for K
N    <- 1e6


Data$t  <- qnorm(1-Data$K)
Data$m  <- sapply(1:nrow(Data),function(k) uniroot(f=function(x) Data[k,"K"] / pnorm(Data[k,"t"],mean=x,lower.tail = F)-Data[k,"Fold"],c(-5,0))$root )
Data$Ke <- Data$K / Data$Fold #pnorm(Data$t,mean=Data$m,lower.tail = F)
K_e     <- K / c(Ke_1,Ke_2,Ke_3)

Rmat  <- do.call("cbind",lapply(1:nrow(Data),function(k){
  sapply(ps, function(p) FoldChangeInGini(m=Data[k,"m"],p=p,K=Data[k,"K"],N=N) )
}))

# plot
ColsD <- c(AD="dodgerblue",MDD="coral1",SCZ="goldenrod",T2D="darkred",CAD="lightgreen")

png("Figure_3_Effect_on_Health_Inequatility_Gini_Index_revised_21Feb2024.png",width=2500,height=2200,res=300)
par(mar=c(5,5,3,2))
matplot(ps,Rmat,type="l",lwd=3,col=ColsD,axes=FALSE,
        ylim=c(0,1.2),lty=1,
        main="",cex.lab=1.2,
        xlab="Proportion of edited genomes in the population",
        ylab="Fold change in Gini Index of disease risk")
axis(1)
axis(2)
abline(h=1,col="grey",lty=2)
abline(v=ps[apply(Rmat,2,which.max)],col=ColsD,lty=3)
legend(0,0.35,legend=
         paste0(c("Alzheimer's","MDD","Schizophrenia",
                  "Type 2 Diabetes","Coronary Artery Disease"),
                paste0(": K' = ",round(Data$K/Data$Fold*100,1),"%")),
       title="Disease Prevalence Among\nEdited Genomes (at 10 loci)",lwd=3,
       box.lty=0,col=ColsD,lty=1,title.font=2)

# legend(0.6,0.4,legend=c("Alzheimer's (5%)","MDD (15%)","Schizophrenia (1%)",
#                        "Type 2 Diabetes (10%)","Coronary Artery Disease (6%)"),col=ColsD,pch=15:19,
#        box.lty=0,cex=0.9,title="Disease (Prevalence)",title.adj = .1,title.col = 2)



dev.off()

## Why sample size?
## Directions: more distributions? Binary? WGS data? More ancestry?
## Anything specific about these false positives?
## False positives: height / BMI from larger GWASs?


