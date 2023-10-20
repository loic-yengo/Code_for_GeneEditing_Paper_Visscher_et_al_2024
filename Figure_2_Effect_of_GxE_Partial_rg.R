setwd("~/Desktop/Papers/Gene-editing/Resubmission - Code - Figure/")
MaxAlleles <- 10
ColsD <- c(ad="dodgerblue",mdd="coral1",scz="goldenrod",t2d="darkred",cad="lightgreen")

## Diseases
ad  <- read.table(file= "sumstats/ad_marioni2018_ukbf.txt", header=TRUE)
mdd <- read.table(file= "sumstats/mdd_howard2019_ukbf.txt", header=TRUE)
scz <- read.table(file= "sumstats/scz_ng2014_ukbf.txt", header=TRUE) # need to add T2D and CAD later
t2d <- read.table(file="sumstats/t2d_Vujkovic2020.txt",h=T,stringsAsFactors = F) # https://www.nature.com/articles/s41588-020-0637-y#Sec45
t2d <- t2d[which(t2d$P<5e-16),]
t2d <- t2d[,c("rsid","EAF","Beta")]
colnames(t2d) <- colnames(ad)
## CAD: https://www.medrxiv.org/content/10.1101/2021.05.24.21257377v1.full.pdf
cad <- read.table("sumstats/cad_aragam2021.txt",h=T,stringsAsFactors = F)
cad <- cad[which(cad[,"P"]<5e-14),]

plotDisease <- function(dt,K,rg){
  t    <- qnorm(1-K)
  z    <- dnorm(t)
  ## Risk Lowering Allelic Frequency
  RLAF <- ifelse(dt[,"Beta"]>0,dt[,"freq_A1"],1-dt[,"freq_A1"])
  BETA <- abs(dt[,"Beta"]) * K * (1 - K) / z
  contribution_to_mean <- 2 * RLAF * BETA * rg
  Y    <- cumsum(sort(contribution_to_mean,decreasing = TRUE))
  K_edited <- pnorm(t,mean=-Y,sd=1,lower.tail = FALSE) # same variance?
  return(K_edited)
}

K_ad     <- 0.05
K_mdd    <- 0.15
K_scz    <- 0.01
K_t2d    <- 0.10
K_cad    <- 0.06

prevTable <- c(T2D=K_t2d,CAD=K_cad,SCZ=K_scz,MDD=K_mdd,AD=K_ad)

bigPlot <- function(rg){
  p_ad     <- plotDisease(dt=ad,K=K_ad,rg=rg)
  p_mdd    <- plotDisease(dt=mdd,K=K_mdd,rg=rg)
  p_scz    <- plotDisease(dt=scz,K=K_scz,rg=rg)
  p_t2d    <- plotDisease(dt=t2d,K=K_t2d,rg=rg)
  p_cad    <- plotDisease(dt=t2d,K=K_cad,rg=rg)
  
  nsnpmax  <- max(c(length(p_ad),length(p_mdd),length(p_scz),length(p_t2d),length(p_cad)))
  diseaseTab <- cbind.data.frame(T2D=c(p_t2d,rep(NA,nsnpmax-length(p_t2d))),
                                 CAD=c(p_cad,rep(NA,nsnpmax-length(p_cad))),
                                 SCZ=c(p_scz,rep(NA,nsnpmax-length(p_scz))),
                                 MDD=c(p_mdd,rep(NA,nsnpmax-length(p_mdd))),
                                 AD=c(p_ad,rep(NA,nsnpmax-length(p_ad))))
  ymin  <- 1
  ymax  <- 50
  par(mar=c(5,5,3,1))
  plot(c(0,10),c(1,ymax),type="n",cex.lab=1.2,log="y",
       main=paste0("Genetic correlation with future environment: r(g) = ",rg),
       xlab="Number of Edited Allele(s)\n(Ranked from largest to smallest expected effect size)",
       ylab="Expected Fold Reduction in Prevalence among Edited Genomes",
       axes=FALSE)
  axis(1,line=-0.5)
  axis(2,at=c(1,2,10,20,30,40,50),labels=c("1x","2x","10x","20x","30x","40x","50x"),las=2)
  
  
  for(i in 1:ncol(diseaseTab)){
    x <- c(1,prevTable[i] / diseaseTab[1:MaxAlleles,i])
    points(0:MaxAlleles,x,pch=14+i,col=ColsD[casefold(colnames(diseaseTab)[i],upper = FALSE)],type="b")
    lines(0:MaxAlleles,x,col=ColsD[casefold(colnames(diseaseTab)[i],upper = FALSE)],lwd=2)
    Ymax <- max(x,na.rm = T)
    Xmax <- which.max(x)
    segments(x0=0,y0=Ymax,x1=Xmax,y1=Ymax,col=ColsD[casefold(colnames(diseaseTab)[i],upper = FALSE)],lty=2)
  }
  if(rg==1){
    legend(0,40,legend=c("Alzheimer's (5%)","MDD (15%)","Schizophrenia (1%)",
                         "Type 2 Diabetes (10%)","Coronary Artery Disease (6%)"),col=ColsD,pch=15:19,
           box.lty=0,cex=1.1,title="Disease (Prevalence)",title.adj = .1,title.col = 2)
  }
  legend(0,50,legend="Higest Reduction",lty=2,box.lty=0)
  
}

png("Figure_2_Effect_of_GxE_Partial_rg.png",width=4000,height=4000,res=300)
op <- par(mfrow=c(2,2))
for(rg in c(1,0.9,0.75,0.5)){
  bigPlot(rg)
}
par(op)
dev.off()
