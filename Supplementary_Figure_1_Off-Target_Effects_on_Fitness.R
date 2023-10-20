setwd("~/Desktop/Papers/Gene-editing/Resubmission - Code - Figure/")
Cols <- c("dodgerblue","coral1","goldenrod","lightgreen","maroon")
## Reduction in fitness
## NC paper: https://www.nature.com/articles/s41467-022-28244-5
## y = 1 - m * pOff * s

fitness <- function(m,pOff,s) 1 - m * pOff * s
ss      <- c(1e-4,2e-4,5e-4,1e-3,2e-3,5e-3,1e-2)
ms      <- c(10,20,50,100)
pOff    <- c(0.05,0.1,0.2,0.5)

scenarios <- expand.grid(s=ss,m=ms,pOff=pOff)
scenarios$fitness <- fitness(m=scenarios$m,pOff=scenarios$pOff,s=scenarios$s)

png("Supplementary_Figure_1_Off-Target_Effects_on_Fitness.png",width=3000,height=2500,res=300)
op <- par(mfrow=c(2,2))
for(i in 1:length(ms)){
  d <- scenarios[which(scenarios$m==ms[i]),]
  z <- NULL
  for(j in 1:length(pOff)){
    z <- cbind(z,d[which(d$pOff==pOff[j]),"fitness"])    
  }
  par(mar=c(5,5,3,2))
  matplot(ss,z,pch=10,type="l",log="x",col=Cols,lwd=2,lty=1,
          main=paste0("Number of loci edited: m = ",ms[i]),
          xlab="Mean selection coefficient at off-target site",
          ylab="Relative fitness among edited genomes",axes=FALSE,ylim=c(0.5,1))
  axis(1,at=ss)
  axis(2)
  legend(ss[1],0.9,legend=paste0(100*pOff,"%"),
         title="Prob. Off-Target Effect",col=Cols,lwd=2,lty=1,box.lty=0)
}
par(op)
dev.off()

