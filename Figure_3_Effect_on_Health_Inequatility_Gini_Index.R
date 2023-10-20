setwd("~/Desktop/Papers/Gene-editing/Resubmission - Code - Figure/")
library(DescTools)
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

## Relationship between K and Gini
# G  <- function(K,N){
#   U  <- seq(1/N,1-1/N,length=N)
#   x  <- qnorm(U)
#   t  <- qnorm(1-K)
#   R0 <- pnorm(x,mean=-t,sd=1,lower.tail = F) # Risk distribution in the population
#   G0 <- Gini(R0, unbiased=FALSE)
#   return(G0)
# }
# Ks <- c(0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5)
# Gs <- sapply(Ks,function(K) G(K,N=1e5))
# 
# par(mar=c(5,5,3,2))
# plot(Ks,Gs,pch=19,log="x",ylim=c(0.2,1),cex.lab=1.2,
#      xlab="Disease Prevalence",axes=FALSE,type='l',lwd=3,
#      ylab="Gini coefficient of risk distribution")
# axis(1,at=Ks)
# axis(2)
# abline(h=1/3,lwd=2,lty=2,col=2)

png("Figure_3_Effect_on_Health_Inequatility_Gini_Index.png",width=3000,height=3000,res=300)
op <- par(mfrow=c(2,2))
for(K in Ks){
  t    <- qnorm(1-K)
  m1   <- uniroot(f=function(x) K / pnorm(t,mean=x,lower.tail = F)-10,c(-5,0))$root
  m2   <- uniroot(f=function(x) K / pnorm(t,mean=x,lower.tail = F)-100,c(-5,0))$root
  m3   <- uniroot(f=function(x) K / pnorm(t,mean=x,lower.tail = F)-1000,c(-5,0))$root
  
  Ke_1 <- pnorm(t,mean=m1,lower.tail = F) # 10x lower
  Ke_2 <- pnorm(t,mean=m2,lower.tail = F) # 100x lower
  Ke_3 <- pnorm(t,mean=m3,lower.tail = F) # 1000x lower
  K_e <- K / c(Ke_1,Ke_2,Ke_3)
    
  R_1   <- sapply(ps, function(p) FoldChangeInGini(m=m1,p=p,K=K,N=N) )
  R_2   <- sapply(ps, function(p) FoldChangeInGini(m=m2,p=p,K=K,N=N) )
  R_3   <- sapply(ps, function(p) FoldChangeInGini(m=m3,p=p,K=K,N=N) )
  Rmat  <- cbind(R_1,R_2,R_3)
  ## plot
  par(mar=c(5,5,3,2))
  matplot(ps,Rmat,type="l",lwd=2,col=Cols,axes=FALSE,
          ylim=c(0,1.2),lty=1:3,
          main=paste0("Disease prevalence: K = ",K," in the population"),
          xlab="Proportion of edited genomes in the population",
          ylab="Fold change in Gini Index of disease risk")
  axis(1)
  axis(2)
  abline(h=1,col="grey",lty=2)
  abline(v=ps[c(which.max(R_1),which.max(R_2),which.max(R_3))],col=Cols,lty=3)
  legend(0,0.5,legend=paste0("K' = K / ",round(K_e)),
         title="Prevalence Among\nEdited Genomes",lwd=2,
         box.lty=0,col=Cols,lty=1:3)
}
par(op)
dev.off()
