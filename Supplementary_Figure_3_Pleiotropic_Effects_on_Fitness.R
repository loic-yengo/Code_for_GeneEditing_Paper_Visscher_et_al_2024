setwd("~/Desktop/Papers/Gene-editing/Resubmission - Code - Figure/")
Cols <- c("dodgerblue","coral1","goldenrod","lightgreen","maroon")
gs <- c(-0.01,-0.02,-0.05,-0.1)
xs <- seq(-4,4,len=100)
fs <- do.call("cbind",lapply(gs,function(g) exp(g*xs^2)))

png("Supplementary_Figure_3_Pleiotropic_Effects_on_Fitness.png",width=2000,height=2000,res=300)
par(mar=c(5,5,3,2))
matplot(xs,fs,col=Cols,pch=19,type="l",lwd=4,lty=1,cex.lab=1.1,
        xlab="Mean phenotypic changes in edited genomes",
        ylab="Fitness of edited genomes relative to unedited genomes",axes=FALSE,ylim=c(0.2,1))
axis(1)
axis(2)
legend(-2,0.5,legend=-gs,
            title=expression(paste("Strength of stabilizing selection (",gamma,")")),
       col=Cols,lwd=2,lty=1,box.lty=0)
dev.off()
