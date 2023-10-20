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

qtTraits <- c("FG","LDL","TG","SBP","DBP")
diseases <- c("AD","MDD","SCZ","T2D","CAD")
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

for(trait in qtTraits){
  dt <- get(trait)
  dt$TIAF <- ifelse(dt[,"BETA"]>0,dt[,"EAF"],1-dt[,"EAF"])
  dt$Contribution <- 2 * abs(dt$BETA) * dt$TIAF
  dt$q2   <- 2 * dt$TIAF * (1 - dt$TIAF) * (dt$BETA^2)
  dt <- dt[order(dt$Contribution,decreasing = TRUE),]
  dt$OrderEditing <- 1:nrow(dt)
  dt <- dt[which(dt$OrderEditing<=nLociMax),]
  assign(trait,dt)
}

## COMMON DISEASES
AD  <- read.table(inputFiles["AD"], header=TRUE)
MDD <- read.table(inputFiles["MDD"], header=TRUE)
SCZ <- read.table(inputFiles["SCZ"], header=TRUE) 
T2D <- read.table(inputFiles["T2D"],h=T,stringsAsFactors = F);T2D <- T2D[,c("rsid","EAF","Beta")];colnames(T2D) <- colnames(AD)
CAD <- read.table(inputFiles["CAD"],h=T,stringsAsFactors = F)

for(disease in diseases){
  dt <- get(disease)
  t    <- qnorm(1-K[disease])
  z    <- dnorm(t)
  dt$TIAF <- ifelse(dt[,"Beta"]>0,dt[,"freq_A1"],1-dt[,"freq_A1"])
  dt$BETA <- abs(dt[,"Beta"]) * K[disease] * (1 - K[disease]) / z
  dt$Contribution <- 2 * dt$TIAF * dt$BETA
  dt$q2 <- 2 * dt$TIAF * (1 - dt$TIAF) * (dt$BETA^2)
  dt <- dt[order(dt$Contribution,decreasing = TRUE),]
  dt$OrderEditing <- 1:nrow(dt)
  dt <- dt[which(dt$OrderEditing<=nLociMax),]
  assign(disease,dt)
}

png("Figure_1_Top_10_loci_edited.png",width=5000,height=2800,res=300)
op <- par(mfrow=c(1,2))
## PLOT DISEASE
diseaseTab <- matrix(NA,nrow=nLociMax,ncol=length(diseases))
colnames(diseaseTab) <- diseases

diseaseTab_min <- matrix(NA,nrow=nLociMax,ncol=length(diseases))
colnames(diseaseTab_min) <- diseases

diseaseTab_max <- matrix(NA,nrow=nLociMax,ncol=length(diseases))
colnames(diseaseTab_max) <- diseases

for(disease in diseases){
  dt <- get(disease)
  Y  <- cumsum(dt[,"Contribution"])
  t  <- qnorm(1-K[disease])
  K_edited <- pnorm(t,mean=-Y,sd=1,lower.tail = FALSE) # same variance?
  se    <- sqrt( cumsum(dt[,"q2"]) )
  Y_min <- Y-se
  Y_max <- Y+se
  K_min <- pnorm(t,mean=-Y_min,sd=1,lower.tail = FALSE)
  K_max <- pnorm(t,mean=-Y_max,sd=1,lower.tail = FALSE)
  diseaseTab[,disease] <- K_edited
  diseaseTab_min[,disease] <- K_min
  diseaseTab_max[,disease] <- K_max
}

ymin  <- 1
ymax  <- 100
ColsD <- c(AD="dodgerblue",MDD="coral1",SCZ="goldenrod",T2D="darkred",CAD="lightgreen")
par(mar=c(5,5,3,1))
plot(c(0,nLociMax),c(1,ymax),type="n",cex.lab=1.2,log="y",
     xlab="Number of Edited Loci\n(Ranked from largest to smallest expected effect size)",
     ylab="Expected Fold Reduction in Prevalence among Edited Genomes",
     axes=FALSE)
axis(1,line=-0.5)
axis(2,at=c(1,10,20,30,40,50,60),labels=c("1x","10x","20x","30x","40x","50x","60x"),las=2)

for(i in 1:ncol(diseaseTab)){
  disease <- diseases[i]
  x <- c(1,K[disease] / diseaseTab[,disease])
  points(0:nLociMax,x,pch=14+i,col=ColsD[disease],type="b")
  lines(0:nLociMax,x,col=ColsD[disease],lwd=2)
  Ymax <- max(x,na.rm = T); Xmax <- which.max(x)
  #segments(x0=0,y0=Ymax,x1=Xmax,y1=Ymax,col=ColsD[disease],lty=2)
  
  ## Add K_min / K_max
  x_min <- c(1,K[disease] / diseaseTab_min[,disease])
  x_max <- c(1,K[disease] / diseaseTab_max[,disease])
  lines(0:nLociMax,x_min,col=ColsD[disease],lty=2)
  lines(0:nLociMax,x_max,col=ColsD[disease],lty=2)
  
}

legend(0,ymax,legend=c("Alzheimer's (5%)","MDD (15%)","Schizophrenia (1%)",
                     "Type 2 Diabetes (10%)","Coronary Artery Disease (6%)"),col=ColsD,pch=15:19,
       box.lty=0,cex=1.1,title="Disease (Prevalence)",title.adj = .1,title.col = 2)
#legend(7,ymax,legend="Higest Reduction",lty=2,box.lty=0)

## PLOT QUANTITATIVE TRAITS
qtTab <- matrix(NA,nrow=nLociMax,ncol=length(qtTraits))
colnames(qtTab) <- qtTraits

qtTab_se <- matrix(NA,nrow=nLociMax,ncol=length(qtTraits))
colnames(qtTab_se) <- qtTraits

for(trait in qtTraits){
  dt <- get(trait)
  qtTab[,trait]    <- cumsum(dt[,"Contribution"])
  qtTab_se[,trait] <- sqrt(cumsum(dt[,"q2"]))
}

ColsQt <- c(FG="dodgerblue",LDL="coral1",TG="goldenrod",SBP="darkred",DBP="lightgreen")
names(ColsQt) <- qtTraits

par(mar=c(5,7,3,1))
plot(c(0,nLociMax),c(15,0.0),type="n",cex.lab=1.3,
     xlab="Number of Edited Loci\n(Ranked from largest to smallest expected effect size)",
     ylab="Expected Reduction Population Mean (in trait SD)\n Among Edited Genomes",
     axes=FALSE)
axis(1,at=c(0,2,4,6,8,nLociMax),line = -0.5)
axis(2,at=seq(0,12,by=2),labels=c(0,seq(-2,-12,by=-2)))
#axis(2,at=,labels=c("0","-1","-2","-5"))
for(i in 1:ncol(qtTab)){
  x <- c(0,qtTab[,i])
  points(0:nLociMax,x,pch=14+i,col=ColsQt[qtTraits[i]],type="b")
  lines(0:nLociMax,x,col=ColsQt[qtTraits[i]],lwd=2)
  Ymax <- max(x,na.rm = T) 
  Xmax <- which.max(x)
  #segments(x0=0,y0=Ymax,x1=Xmax,y1=Ymax,col=ColsQt[qtTraits[i]],lty=2)
  ## Add SE's or SD of gain
  for(j in 1:nrow(qtTab)){
    segments(j,qtTab[j,i]-qtTab_se[j,i],j,qtTab[j,i]+qtTab_se[j,i],col=ColsQt[qtTraits[i]])
  }
  x_min <- c(0,qtTab[,i]-qtTab_se[,i])
  x_max <- c(0,qtTab[,i]+qtTab_se[,i])
  lines(0:nLociMax,x_min,col=ColsQt[qtTraits[i]],lty=2)
  lines(0:nLociMax,x_max,col=ColsQt[qtTraits[i]],lty=2)
}

legend(0,15,legend=qtTraits,col=ColsQt,pch=14+(1:nLociMax),box.lty=0,cex=1.1,lwd=2,lty=1)
#legend(7,15,legend="Highest Reduction",lty=2,box.lty=0)

par(op)
dev.off()

Dt <- NULL
for(trait in qtTraits){
  dt       <- get(trait)[,c("SNP","TIAF","BETA","Contribution","q2","OrderEditing")]
  dt$Trait <- trait
  dt$CumulativeMean <- cumsum(dt$Contribution)
  dt$CumulativeVarExp <- cumsum(dt$q2)
  dt$SD_gain <- sqrt(dt$CumulativeVarExp)
  dt   <- dt[,c("Trait","SNP","TIAF","BETA","Contribution","q2","OrderEditing",
                "CumulativeMean","CumulativeVarExp","SD_gain")]
  dt$BETA <- abs(dt$BETA)
  colnames(dt) <- c("Trait","SNP","Trait-increasing Allele Frequency","Trait-increasing Effect Size",
                    "Gain","Variance Explained","Order of Editing","Expected Cumulative Gain","Cumulative Variance Explained",
                    "SD of Gain")
  Dt         <- rbind(Dt,dt)
}

Ds <- NULL
for(disease in diseases){
  ds       <- get(disease)[,c("SNP","TIAF","Beta","Contribution","q2","OrderEditing")]
  ds$K     <- K[disease]
  ds$Trait <- disease
  ds$CumulativeMean <- cumsum(ds$Contribution)
  ds$CumulativeVarExp <- cumsum(ds$q2)
  ds$SD_gain <- sqrt(ds$CumulativeVarExp)
  Y  <- cumsum(ds[,"Contribution"])
  ds$Prevalence_Ratio <- K[disease]/pnorm(qnorm(1-K[disease]),mean=-ds$CumulativeMean,sd=1,lower.tail = FALSE)
  
  ## Added on 20/10/2023
  ds <- ds[,c("Trait","K","SNP","TIAF","Beta","Contribution","q2","OrderEditing","CumulativeMean","CumulativeVarExp","SD_gain","Prevalence_Ratio")]
  ds[,"Beta"] <- exp(abs(ds[,"Beta"]))
  colnames(ds) <- c("Disease","Assumed Prevalence","SNP","Risk Allele Frequency","Odds Ratio","Gain","Variance Explained (liability scaled)",
                    "Order of Editing","Expected Cumulative Gain","Cumulative Variance Explained (liability scale)",
                    "SD of Gain","Fold-reduction in Prevalence")
  Ds       <- rbind(Ds,ds)
}


write.table(Dt,"Example_of_10_loci_to_edit_for_5_qt_traits.txt",quote=F,row.names = F,sep="\t")
write.table(Ds,"Example_of_10_loci_to_edit_for_5_diseases.txt",quote=F,row.names = F,sep="\t")

