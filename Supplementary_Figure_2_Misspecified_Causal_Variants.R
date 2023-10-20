setwd("~/Desktop/Papers/Gene-editing/Resubmission - Code - Figure/")
## Editing here aims at reducing risk
## QT are FG, LDL, TG, SBP and DBP (x5)
## Diseases are AD, MDD, SCZ, T2D, CAD (x5)

inputFiles <- c(
  FG="sumstats/FG_36snps.txt",
  GLGC="sumstats/GLGC_Graham2021.txt",
  ICBP="sumstats/ICBP_Warren2022.txt",
  AD="sumstats/ad_marioni2018_ukbf.txt",
  MDD="sumstats/mdd_howard2019_ukbf.txt",
  T2D="sumstats/t2d_Vujkovic2020.txt",# https://www.nature.com/articles/s41588-020-0637-y#Sec45
  CAD="sumstats/cad_aragam2021.txt",# https://www.medrxiv.org/content/10.1101/2021.05.24.21257377v1.full.pdf   
  SCZ="sumstats/scz_ng2014_ukbf.txt"
)

qtTraits <- sort(c("FG","LDL","TG","SBP","DBP"))
diseases <- sort(c("AD","MDD","SCZ","T2D","CAD"))
K        <- c(AD=0.05,MDD=0.15,SCZ=0.01,T2D=0.1,CAD=0.06)
nLociMax <- 10

## QUANTITATIVE TRAITS
fg    <- read.table(inputFiles["FG"],h=T,stringsAsFactors = F,sep="\t")
Vp_fg <- mean( fg$N * 2 * fg$EAF * (1-fg$EAF) * (fg$SE^2) )
sd_fg <- sqrt(Vp_fg)
FG    <- fg[,c("SNP","EAF","BETA")]
FG[,"BETA"] <-  FG[,"BETA"] / sd_fg

glgc <- read.table(inputFiles["GLGC"],h=T,stringsAsFactors = F,sep="\t")
glgc <- glgc[which(glgc$Ancestry=="EUR"),]
ldl  <- glgc[which(glgc$Lipid=="LDL"),]
tg   <- glgc[which(glgc$Lipid=="logTG"),]

sd_ldl <- sqrt( mean( ldl$N * 2 * ldl$AF * (1-ldl$AF) * (ldl$Effectsize_SD^2) ) )
sd_tg  <- sqrt( mean( tg$N * 2 * tg$AF * (1-tg$AF) * (tg$Effectsize_SD^2) ) )

LDL <- ldl[,c("rs_dbSNP150","AF","Effectsize")]; colnames(LDL) <- c("SNP","EAF","BETA")
TG  <- tg[,c("rs_dbSNP150","AF","Effectsize")]; colnames(TG) <- c("SNP","EAF","BETA")

#LDL <- LDL[which(LDL$EAF>0.001 & LDL$EAF<0.999),]
#TG  <- TG[which(TG$EAF>0.001 & TG$EAF<0.999),]

icbp <- read.table(inputFiles["ICBP"],h=T,stringsAsFactors = F,sep="\t",comment.char = "!")
icbp <- icbp[which(icbp$CHR%in%(1:22)),]

sbp <- na.omit(icbp[which(icbp$Trait=="SBP"),])
dbp <- na.omit(icbp[which(icbp$Trait=="DBP"),])

sd_sbp <- sqrt( mean( sbp$Neff * 2 * sbp$EAF * (1-sbp$EAF) * (sbp$SE^2) ) )
sd_dbp <- sqrt( mean( dbp$Neff * 2 * dbp$EAF * (1-dbp$EAF) * (dbp$SE^2) ) )

SBP <- sbp[,c("rsID","EAF","Effect")]; colnames(SBP) <- c("SNP","EAF","BETA"); SBP[,"BETA"] <- SBP[,"BETA"] / sd_sbp
DBP <- dbp[,c("rsID","EAF","Effect")]; colnames(DBP) <- c("SNP","EAF","BETA"); DBP[,"BETA"] <- DBP[,"BETA"] / sd_dbp

# 1, 2 or 5
cmb1 <- combn(nLociMax,1)
cmb2 <- combn(nLociMax,2)
cmb5 <- combn(nLociMax,5)

gmb1 <- apply(cmb1,2,toString)
gmb2 <- apply(cmb2,2,toString)
gmb5 <- apply(cmb5,2,toString)

ReductionQt <- NULL
for(trait in sort(qtTraits)){
  dt <- get(trait)
  dt$TIAF <- ifelse(dt[,"BETA"]>0,dt[,"EAF"],1-dt[,"EAF"])
  dt$Contribution <- 2 * abs(dt$BETA) * dt$TIAF
  dt <- dt[order(dt$Contribution,decreasing = TRUE),]
  dt$OrderEditing <- 1:nrow(dt)
  dt   <- dt[which(dt$OrderEditing<=nLociMax),]
  St   <- sum(dt[,"Contribution"])
  dt_1 <- sapply(1:ncol(cmb1),function(j) sum(dt[-cmb1[,j],"Contribution"]) / St)
  dt_2 <- sapply(1:ncol(cmb2),function(j) sum(dt[-cmb2[,j],"Contribution"]) / St)
  dt_5 <- sapply(1:ncol(cmb5),function(j) sum(dt[-cmb5[,j],"Contribution"]) / St)
  ReductionQt <- rbind(ReductionQt,
                       rbind(
                         cbind.data.frame(Trait=trait,Nmissed=1,FracOfMaxEff=dt_1),
                         cbind.data.frame(Trait=trait,Nmissed=2,FracOfMaxEff=dt_2),
                         cbind.data.frame(Trait=trait,Nmissed=5,FracOfMaxEff=dt_5)
                         )
  )
}

## COMMON DISEASES
AD  <- read.table(inputFiles["AD"], header=TRUE)
MDD <- read.table(inputFiles["MDD"], header=TRUE)
SCZ <- read.table(inputFiles["SCZ"], header=TRUE) 
T2D <- read.table(inputFiles["T2D"],h=T,stringsAsFactors = F);T2D <- T2D[,c("rsid","EAF","Beta")];colnames(T2D) <- colnames(AD)
CAD <- read.table(inputFiles["CAD"],h=T,stringsAsFactors = F)

ReductionDs <- NULL
for(disease in sort(diseases)){
  dt      <- get(disease)
  t       <- qnorm(1-K[disease])
  z       <- dnorm(t)
  dt$TIAF <- ifelse(dt[,"Beta"]>0,dt[,"freq_A1"],1-dt[,"freq_A1"])
  dt$BETA <- abs(dt[,"Beta"]) * K[disease] * (1 - K[disease]) / z
  dt$Contribution <- 2 * dt$TIAF * dt$BETA
  dt <- dt[order(dt$Contribution,decreasing = TRUE),]
  dt$OrderEditing <- 1:nrow(dt)
  dt <- dt[which(dt$OrderEditing<=nLociMax),]
  t  <- qnorm(1-K[disease])
  K_edited <- pnorm(t,mean=-sum(dt[,"Contribution"]),sd=1,lower.tail = FALSE) # same variance?

  dt_1 <- sapply(1:ncol(cmb1),function(j) K_edited / pnorm(t,mean=-sum(dt[-cmb1[,j],"Contribution"]),sd=1,lower.tail = FALSE) )
  dt_2 <- sapply(1:ncol(cmb2),function(j) K_edited / pnorm(t,mean=-sum(dt[-cmb2[,j],"Contribution"]),sd=1,lower.tail = FALSE) )
  dt_5 <- sapply(1:ncol(cmb5),function(j) K_edited / pnorm(t,mean=-sum(dt[-cmb5[,j],"Contribution"]),sd=1,lower.tail = FALSE) )
  ReductionDs <- rbind(ReductionDs,
                       rbind(
                         cbind.data.frame(Disease=disease,Nmissed=1,FracOfMaxEff=dt_1),
                         cbind.data.frame(Disease=disease,Nmissed=2,FracOfMaxEff=dt_2),
                         cbind.data.frame(Disease=disease,Nmissed=5,FracOfMaxEff=dt_5)
                       )
  )
}


png("Supplementary_Figure_2_Misspecified_Causal_Variants.png",width=5000,height=2000,res=300)
op <- par(mfrow=c(1,2))
## PLOT DISEASE
ColsDs <- c(AD="dodgerblue",CAD="lightgreen",MDD="coral1",SCZ="goldenrod",T2D="darkred")

par(mar=c(5,7,3,1))
boxplot(FracOfMaxEff~Disease+Nmissed,data=ReductionDs,cex.lab=1.2,col=ColsDs,
        main="Common Disease",
        axes=FALSE,xlab="Number of missed causal variants out of 10 loci targeted",
        ylab="Fraction of maximum efficiency\nrelative to if all 10 causal variants were edited",
        ylim=c(0,1))
axis(1,at=c(3,8,13),labels=c("n = 1\n(10 combinations)","n = 2\n(45 combinations)","n = 5\n(252 combinations)"),tick=FALSE)
axis(2)
abline(v=c(5.5,10.5),col="grey",lty=2)
abline(h=1-c(1,2,5)/nLociMax,col="lightblue",lty=3)
k <- 0
for(n in c(1,2,5)){
  for(disease in sort(diseases)){
    k <- k + 1
    y <- ReductionDs[which(ReductionDs$Disease==disease & ReductionDs$Nmissed==n),"FracOfMaxEff"]
    x <- rep(k,length(y))
    points(x,y,pch=19,col=ColsDs[disease],cex=0.8)
  }
}
legend(0,0.3,legend=diseases,col=ColsDs,pch=19,box.lty=0,cex=1.1,lwd=2,lty=0)


## PLOT QUANTITATIVE TRAITS
ColsQt <- c(DBP="lightgreen",FG="dodgerblue",LDL="coral1",SBP="darkred",TG="goldenrod")

par(mar=c(5,7,3,1))
boxplot(FracOfMaxEff~Trait+Nmissed,data=ReductionQt,cex.lab=1.2,col=ColsQt,
        main="Quantitative Disease Risk Factors",
        axes=FALSE,xlab="Number of missed causal variants out of 10 loci targeted",
        ylab="Fraction of maximum efficiency\nrelative to if all 10 causal variants were edited",
        ylim=c(0,1))
axis(1,at=c(3,8,13),labels=c("n = 1\n(10 combinations)","n = 2\n(45 combinations)","n = 5\n(252 combinations)"),tick=FALSE)
axis(2)
abline(v=c(5.5,10.5),col="grey",lty=2)
abline(h=1-c(1,2,5)/nLociMax,col="lightblue",lty=3)
k <- 0
for(n in c(1,2,5)){
  for(trait in qtTraits){
    k <- k + 1
    y <- ReductionQt[which(ReductionQt$Trait==trait & ReductionQt$Nmissed==n),"FracOfMaxEff"]
    x <- rep(k,length(y))
    points(x,y,pch=19,col=ColsQt[trait],cex=0.8)
  }
}
legend(11,1,legend=qtTraits,col=ColsQt,pch=19,box.lty=0,cex=1.1,lwd=2,lty=0)
par(op)
dev.off()

## Numerical values

## Disease
meanReductDisease <- aggregate(FracOfMaxEff~Disease+Nmissed,data=ReductionDs,FUN=mean) # mean reduction
sdReductDisease   <- aggregate(FracOfMaxEff~Disease+Nmissed,data=ReductionDs,FUN=sd) # SD of reduction

## Quantitative traits
meanReductQt <- aggregate(FracOfMaxEff~Trait+Nmissed,data=ReductionQt,FUN=mean) # mean reduction
sdReductQt   <- aggregate(FracOfMaxEff~Trait+Nmissed,data=ReductionQt,FUN=sd) # SD of reduction


