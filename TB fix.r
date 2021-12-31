# kanittaTB.Rcm

library(MASS)
library(rpart)
library(grid)
library(partykit)
library(epiDisplay)

setwd("C:/Users/USER/Dropbox/project/ICRIEMS") 
read.table("va.txt",h=T,as.is=T) -> m
str(m)

#data management 
m$age <- ifelse(is.na(m$age),99,m$age)          # replace 2 unknown ages by 99
m <- subset(m,m$age >= 5)                       # age at least 5 

m$VA <- ifelse(m$VA=="O758","N758",m$VA)# correct erroneous VA
m$VA <- ifelse(m$VA=="O758NA","N758NA",m$VA)# correct erroneous VA

cgs <- c("TB","se","HIV","oID",			# VR cause groups
	"liCa","luCa","oDC","oCa",
	"en","MN","IHD","str","oCV",
	"res","dig","gU","iDf","trA",
	"oIj","sui","oth")
cgA <- c("a","b","c","d","e","f",
	"g","h","i","j","k","l",
	"m","n","o","p","q","r","s","t","u")

vaIDs <- c(1:21)
vaGrps <- paste("0",vaIDs[1:9],":",cgA[1:9],sep="")  # put zero in front of single numbers (1-9)
vaGrps <- c(vaGrps,paste(vaIDs[10:21],":",cgA[10:21],sep=""))
Ch <- toupper(substr(m$VA,1,1))# Define VA cause groups
Int <- as.integer(substr(m$VA,2,3))

m$VAgrp <- ifelse(Ch=="A" & Int %in% c(15:19),vaGrps[1],
  ifelse(Ch=="A" & Int %in% c(40:41),vaGrps[2],
  ifelse(Ch=="B" & Int %in% c(20:24),vaGrps[3],
  ifelse(Ch %in% c("A","B"),vaGrps[4],
  ifelse(Ch=="C" & Int==22,vaGrps[5],
  ifelse(Ch=="C" & Int %in% c(30:39),vaGrps[6],
  ifelse(Ch=="C" & Int %in% c(15:21,23:26),vaGrps[7],
  ifelse(Ch=="C" | (Ch=="D" & Int<50),vaGrps[8], 
  ifelse(Ch=="E",vaGrps[9],ifelse(Ch %in% c("F","G"),vaGrps[10], 
  ifelse(Ch=="I" & Int %in% c(20:25),vaGrps[11],
  ifelse(Ch=="I" & Int %in% c(60:69),vaGrps[12],
  ifelse(Ch=="I",vaGrps[13],ifelse(Ch=="J",vaGrps[14],
  ifelse(Ch=="K",vaGrps[15],ifelse(Ch=="N",vaGrps[16],
  ifelse(Ch=="R",vaGrps[17],ifelse(Ch=="V",vaGrps[18],
  ifelse(Ch=="W" | (Ch=="X" & Int<60),vaGrps[19],
  ifelse(Ch=="X" & Int %in% c(60:84),vaGrps[20],
  vaGrps[21]))))) ))))) ))))) ))))) 

Ch <- toupper(substr(m$ncauseDR,1,1))# Define VR cause groups
Int <- as.integer(substr(m$ncauseDR,2,3))
m$DRgrp <- ifelse(Ch=="A" & Int %in% c(15:19),vaGrps[1],
  ifelse(Ch=="A" & Int %in% c(40:41),vaGrps[2],
  ifelse(Ch=="B" & Int %in% c(20:24),vaGrps[3],
  ifelse(Ch %in% c("A","B"),vaGrps[4],
  ifelse(Ch=="C" & Int==22,vaGrps[5],
  ifelse(Ch=="C" & Int %in% c(30:39),vaGrps[6],
  ifelse(Ch=="C" & Int %in% c(15:21,23:26),vaGrps[7],
  ifelse(Ch=="C" | (Ch=="D" & Int<50),vaGrps[8], 
  ifelse(Ch=="E",vaGrps[9],ifelse(Ch %in% c("F","G"),vaGrps[10], 
  ifelse(Ch=="I" & Int %in% c(20:25),vaGrps[11],
  ifelse(Ch=="I" & Int %in% c(60:69),vaGrps[12],
  ifelse(Ch=="I",vaGrps[13],ifelse(Ch=="J",vaGrps[14],
  ifelse(Ch=="K",vaGrps[15],ifelse(Ch=="N",vaGrps[16],
  ifelse(Ch=="R",vaGrps[17],ifelse(Ch=="V",vaGrps[18],
  ifelse(Ch=="W" | (Ch=="X" & Int<60),vaGrps[19],
  ifelse(Ch=="X" & Int %in% c(60:84),vaGrps[20],
  vaGrps[21]))))) ))))) ))))) ))))) 

addmargins(table(m$VAgrp,m$DRgrp))

nv <- ncol(m)
for (i in c(1:21)) {
 m[,nv+i] <- ifelse(as.integer(substr(m$DRgrp,1,2))==i,1,0)
 names(m)[nv+i] <- paste("DR",i,sep="")
}

#manage the province
nv <- ncol(m)
m$pro <- ifelse(m$pro==10,1,
	ifelse(m$pro==34,2,
	ifelse(m$pro==57,3,
	ifelse(m$pro==72,4,
 	ifelse(m$pro==90,5,
	ifelse(m$pro==26,6,
	ifelse(m$pro==42,7,
	ifelse(m$pro==56,8,
	ifelse(m$pro==86,9,0)))) )))) )

for (i in c(1:length(unique(m$pro)))) {
 m[,nv+i] <- ifelse(m$pro==i,1,0)
 names(m)[nv+i] <- paste("prov",i,sep="")
}

#tree based method
#m$tb <- as.factor(ifelse(m$VAgrp=="01:a","Tuberculosis","other"))

m$tb <- ifelse(m$VAgrp=="01:a",1,0)

names(m)
m.rp <- rpart(tb~.,m[,c(2:4,9:38,39)],cp=0.0001)
printcp(m.rp)

windows(10,6)
par(las=1,mar=c(4,4,4,2))
plotcp(m.rp)

windows(15,7)
#par(mar=c(1,0,1,0))
plot(m.rp,uniform=T); text(m.rp,use.n=T,cex=1,digits=3,font=2)
#m$gp0 <- predict(as.party(m.rp),type="node")	

#glm(family=binomial,data=m,tb~factor(gp0)) -> mod0
#summary(mod0)
#lroc(mod0)$auc

#pruning the tree
m.rp1 <- prune(m.rp,cp=0.0007)
windows(10,6)
par(las=1,mar=c(4,4,4,2))
plotcp(m.rp1)

windows(12,8)
par(mar=c(1,1,1,1))
plot(m.rp1,uniform=T); text(m.rp1,use.n=T,cex=1,font=2)

m.rp2 <- prune(m.rp,cp=0.0012)
windows(12,8)
par(mar=c(1,1,1,1))
plot(m.rp2,uniform=T); text(m.rp2,use.n=T,cex=1,font=2)

windows(10,6)
par(las=1,mar=c(4,4,4,2))
plotcp(m.rp2)
printcp(m.rp2)

m$gp <- predict(as.party(m.rp2),type="node")	

addmargins(table(m$gp))

glm(family=binomial,data=m,tb~factor(gp)) -> mod
summary(mod)
#lroc(mod)$auc

#-------------------------------------------------------------------------#
# plot confidence intervals 

source("C:/Users/USER/Dropbox/project/ICRIEMS/dcis.Rcm")

names(m)
glm.dci(m,yID=39,xIDs=40,delta=0.1) -> rez
#the model is not fitted well

xlab1 <- "Tree-based Group"
ylab <- "Tuberculosis Death Prevalence (%)"
titl <- paste(nrow(m)," Verbal Autopsies in Thailand: 2005",sep="")

windows(8,4)

par(oma=c(0,0,0,0),mar=c(2.5,2,2,1),las=1,mgp=c(1.1,0.1,0),tcl=0.2)

n1 <- length(unique(m$gp))

xCoord <- c(1:n1)

pc1 <- rez[c(1:n1),3]

yCoord <- 10*sqrt(pc1)
cilb1 <- rez[1:n1,4]
ciub1 <- rez[1:n1,5]

ymin <- -4
ymax <- 100

plot(1,type="n",xlim=c(1,max(xCoord)),ylim=c(ymin,ymax),
ylab="",xlab="",xaxt="n",yaxt="n")
abline(h=10*sqrt(c(1,5,15,20,25,40,50,70,80,90)),col=8)	#
meanPc <- 100*mean(m$tb)

yLab <- c(0,50,100)
yTick <- 10*sqrt(yLab)
axis(side=2,at=yTick,lab=yLab)

abline(h=10*sqrt(meanPc),col=2)
abline(h=yTick,col="grey")

for (i in c(1:n1)) {
points(xCoord[i]+c(0,0),10*sqrt(c(cilb1[i],ciub1[i])),type="l",lwd=2)
}
points(xCoord,yCoord,pch=20)


text(c(1:n1),10*sqrt(cilb1),adj=c(0.5,1),table(m$gp))

axis(side=1,at=c(1:n1),lab=c(1:n1))
axis(side=1,at=2*c(1:15),lab=2*c(1:15))

at1 <- (1+n1)/2
axis(side=1,at=at1,lab=xlab1,tcl=0,padj=1.4)

mtext(side=3,line=0.1,adj=-0.055,ylab)
mtext(side=3,line=0.1,adj=1,titl)


points(xCoord[c(3,7,8,13)],yCoord[c(3,7,8,13)],type="l",lwd=2,col=6)
points(xCoord[c(4,9,10,11,12,14,15,16)],yCoord[c(4,9,10,11,12,14,15,16)],type="l",lwd=2,col=3)
points(xCoord[c(17,18,19,20)],yCoord[c(17,18,19,20)],type="l",lwd=2,col=4)

# regroup by mustering

gp <- m$gp
m$gp1 <- ifelse(gp==4,1,ifelse(gp==7,2,	
	ifelse(gp %in% c(9,15,19,27),5,
	ifelse(gp==11,4,ifelse(gp==13,3,
	ifelse(gp %in% c(10,20,21,23,24,28,31,34),6,
	ifelse(gp %in% c(35,37,38,39),7,8))) ))))	

glm(family=binomial,data=m,tb~factor(gp1)) -> mod
summary(mod)
drop1(mod, type="Chisq")-> rez1
rez1

pval <- rez1$"Pr(>Chi)"[2:3]
pval1 <- ifelse(pval[1]<0.0001,"<0.0001",round(pval[1],4))

names(m)
glm.dci(m,yID=39,xIDs=41,delta=0.1) -> rez

xlab <- "Tree-based Group"
ylab <- "Tuberculosis Death Prevalence (%)"
titl <- paste(nrow(m)," Verbal Autopsies in Thailand: 2005",sep="")

windows(6,4)

par(oma=c(0,0,0,0),mar=c(2.5,2,2,1),las=1,mgp=c(1.1,0.1,0),tcl=0.2)

n <- length(unique(m$gp1))

xCoord <- c(1:n)

pc <- rez[c(1:n),3]

yCoord <- 10*sqrt(pc)
cilb <- rez[1:n,4]
ciub <- rez[1:n,5]

ymin <- -4
ymax <- 100

abline(h=10*sqrt(meanPc),col=2)

plot(1,type="n",xlim=c(0.5,max(xCoord)+0.5),ylim=c(ymin,ymax),
	ylab="",xlab="",xaxt="n",yaxt="n")
meanPc <- 100*mean(m$tb)

yLab <- c(0,2,10,30,60,100)
yTick <- 10*sqrt(yLab)
axis(side=2,at=yTick,lab=yLab)

abline(h=10*sqrt(meanPc),col=2)
abline(h=yTick,col="grey")

for (i in c(1:n)) {
 points(xCoord[i]+c(0,0),10*sqrt(c(cilb[i],ciub[i])),type="l",lwd=2)
}
points(xCoord,yCoord,pch=20)

text(c(1:n),10*sqrt(cilb),adj=c(0.5,1),table(m$gp1))

axis(side=1,at=c(1:n),lab=c(1:n))
axis(side=1,at=2*c(1:15),lab=2*c(1:15))

at1 <- (1+n)/2
axis(side=1,at=at1,lab=xlab,tcl=0,padj=1.4)

mtext(side=3,line=0.1,adj=-0.1,ylab)
mtext(side=3,line=0.1,adj=1,titl)

#--------------------------------------------
# crude and adjusted OR
# plot unadjusted means

pcCr <- 100*tapply(m$tb,m$gp1,mean)	#crude

points(xCoord-0.2,10*sqrt(pcCr),pch=21,bg=3,cex=0.8)

axis(side=1,at=c(1:n),lab=c(1:n))
axis(side=1,at=2*c(1:7),lab=2*c(1:7))

at <- (1+n1)/2
axis(side=1,at=at,lab=xlab,tcl=0,padj=1.4)

ptc <- c("1","2","3","4","5","6","7","8","9")
lg <- c("Crude","Adjusted")
legend("topleft",inset=c(0.34,0.08),leg=lg,pch=c(21,20),pt.bg=c(3,NA),pt.cex=c(0.8,1.2),bg="ivory",x.intersp=0.4)

#-----------------------------------------
# alternative method, simply using logistic regression
ag <- m$age; sx <- m$sex; gh <- m$ghos; tb <- m$DR1
rs <- m$DR14; id <- m$DR17
m$agSx <- ifelse(ag<40 & sx==1,1,ifelse(ag<40 & sx==2,2,
	ifelse(ag<50 & sx==1,3,ifelse(ag<50 & sx==2,4,
	ifelse(ag<60 & sx==1,5,ifelse(ag<60 & sx==2,6,
	ifelse(ag<70 & sx==1,7,ifelse(ag<70 & sx==2,8,
	ifelse(ag<80 & sx==1,9,ifelse(ag<80 & sx==2,10,
	ifelse(ag>=80 & sx==1,11,ifelse(ag>=80 & sx==2,12,13))))) ))))) ))
addmargins(table(m$agSx))

od13 <- c("b","c","d","e","f","g","h","i","j","k",
	"l","m","o","p","r","s","t","u")
od <- ifelse(substr(m$DRgrp,4,4) %in% od13,1,0)

m$drGh <- ifelse(tb==1 & gh==1,1,ifelse(tb==1 & gh==2,2,
	ifelse(rs==1 & gh==1,3,ifelse(rs==1 & gh==2,4,
	ifelse(id==1 & gh==1,5,ifelse(id==1 & gh==2,6,
	ifelse(od==1 & gh==1,7,ifelse(od==1 & gh==2,8,9)) ))) ))) 

addmargins(table(m$drGh))

glm(family=binomial,data=m,tb~factor(pro)+factor(agSx)+factor(drGh)) -> mod
summary(mod)
drop1(mod,test="Chisq")

#lroc(mod)$auc
# Democratic confidence intervals

options(scipen=8)

drop1(mod,test="Chisq") -> rez1
pval <- rez1$"Pr(>Chi)"[2:4]
pval1 <- ifelse(pval[1]<0.0001,"<0.0001",round(pval[1],4))
pval2 <- ifelse(pval[2]<0.0001,"<0.0001",round(pval[2],4))
pval3 <- ifelse(pval[3]<0.0001,"<0.0001",round(pval[3],4))

names(m)
glm.dci(m,yID=39,xIDs=c(1,42,43),reverse=F,delta=0.001) -> rez

xlab1 <- "Province"
xlab2 <- "Gender & Age Group"
xlab3 <- "Location & Reported Cause"
ylab <- "Tuberculosis Death Prevalence (%)"
titl <- paste(nrow(m)," Verbal Autopsies in Thailand: 2005",sep="")
 
windows(9,4)

par(oma=c(0,0,0,0),mar=c(2.5,2,2,1),las=1,mgp=c(1.1,0.1,0),tcl=0.2)

n1 <- length(unique(m$pro))
n2 <- length(unique(m$agSx))
n3 <- length(unique(m$drGh))

xCoord <- c((1:n1),(n1+2):(n1+n2+1),(n1+n2+3):(n1+n2+n3+2))

pc1 <- rez[c(1:n1),3]
pc2 <- rez[(n1+1):(n1+n2),3]
pc3 <- rez[(n1+n2+1):(n1+n2+n3),3]

yCoord <- 10*sqrt(c(pc1,pc2,pc3))
cilb1 <- rez[1:n1,4]
cilb2 <- rez[(n1+1):(n1+n2),4]
cilb3 <- rez[(n1+n2+1):(n1+n2+n3),4]

ciub1 <- rez[1:n1,5]
ciub2 <- rez[(n1+1):(n1+n2),5]
ciub3 <- rez[(n1+n2+1):(n1+n2+n3),5]

ymin <- 0
ymax <- 100

plot(1,type="n",xlim=c(0.5,max(xCoord)+0.5),ylim=c(ymin,ymax),
	ylab="",xlab="",xaxt="n",yaxt="n")
meanPc <- 100*mean(m$tb)

yLab <- c(0,2,10,30,60,100)
yTick <- 10*sqrt(yLab)
axis(side=2,at=yTick,lab=yLab)

abline(h=10*sqrt(meanPc),col=2)
abline(h=yTick,col="grey")
abline(v=c(n1+1,n1+n2+2),col="dimgrey")

for (i in c(1:n1)) {
 points(xCoord[i]+c(0,0),10*sqrt(c(cilb1[i],ciub1[i])),type="l",lwd=2)
}
for (i in c(1:n2)) {
 points(xCoord[n1+i]+c(0,0),10*sqrt(c(cilb2[i],ciub2[i])),type="l",lwd=2)
}
for (i in c(1:n3)) {
 points(xCoord[n1+n2+i]+c(0,0),10*sqrt(c(cilb3[i],ciub3[i])),type="l",lwd=2)
}
points(xCoord[(n1+n2+1):(n1+n2+n3)],yCoord[(n1+n2+1):(n1+n2+n3)],type="l")
points(xCoord[1:n1],yCoord[1:n1],pch=20)
points(xCoord[(n1+1):(n1+n2)],yCoord[(n1+1):(n1+n2)],pch=21,bg=c(4,2))
points(xCoord[(n1+n2+1):(n1+n2+n3)],yCoord[(n1+n2+1):(n1+n2+n3)],pch=21,bg=c("grey60","orange"))

# plot unadjusted means

pcCr1 <- 100*tapply(m$tb,m$pro,mean)	# crude percentages
pcCr2 <- 100*tapply(m$tb,m$agSx,mean)
pcCr3 <- 100*tapply(m$tb,m$drGh,mean)

points(xCoord[1:n1]-0.2,10*sqrt(pcCr1),pch=21,bg=3,cex=0.8)
points(xCoord[(n1+1):(n1+n2)]-0.2,10*sqrt(pcCr2),pch=21,bg=3,cex=0.8)
points(xCoord[(n1+n2+1):(n1+n2+n3)]-0.2,10*sqrt(pcCr3),pch=21,bg=3,cex=0.8)

axis(side=1,at=c(1:n1),lab=c(1:n1))
axis(side=1,at=2*c(1:4),lab=2*c(1:4))
axis(side=1,at=n1+0.5+2*c(1:(n2/2)),lab=c(1:(n2/2)))
axis(side=1,at=n1+n2+1.5+2*c(1:(n3/2)),lab=c(1:(n3/2)))

at1 <- (1+n1)/2
axis(side=1,at=at1,lab=xlab1,tcl=0,padj=1.4)
at2 <- (n1+2+n1+n2+1)/2
axis(side=1,at=at2,lab=xlab2,tcl=0,padj=1.4)
at3 <- (n1+n2+3+n1+n2+n3+2)/2
axis(side=1,at=at3,lab=xlab3,tcl=0,padj=1.4)

text(at1,ymin+2,adj=c(0.5,0),paste("p-value",pval1,sep=""))
text(at2,ymin+2,adj=c(0.5,0),paste("p-value",pval2,sep=""))
text(at3,ymin+2,adj=c(0.5,0),paste("p-value",pval3,sep=""))

mtext(side=3,line=0,adj=-0.05,ylab)
mtext(side=3,line=0.2,adj=1,titl)

ptc <- c("1","2","3","4","5","6","7","8","9")
lg1 <- c("Crude","Adjusted")
legend("topleft",inset=c(0.005,0.08),leg=lg1,pch=c(21,20),pt.bg=c(3,NA),pt.cex=c(0.8,1.2),bg="ivory",x.intersp=0.4)
lg2 <- c("Bangkok","Nak.Nayak","UbonRat.","Loei","Phayao","ChiangRai","SupanBuri","Chumpon","Songkla")
legend("topleft",inset=c(0.12,0.04),leg=lg2,pch=ptc,cex=0.9,bg="ivory",y.intersp=0.8)
lg3 <- c("Male","Female")
legend("topleft",inset=c(0.26,0.08),leg=lg3,pch=21,pt.bg=c(4,2),pt.cex=1,bg="ivory",x.intersp=0.8)
lg4 <- c(" 5-39","40-49","50-59","60-69","70-79"," 80+")
legend("topleft",inset=c(0.44,0.04),leg=lg4,pch=ptc[1:6],cex=0.9,bg="ivory",y.intersp=0.8)
lg5 <- c("In Hospital","Elsewhere")
legend("topright",inset=c(0.22,0.08),leg=lg5,pch=21,pt.bg=c("grey60","orange"),pt.cex=1,bg="ivory",x.intersp=0.8)
lg6 <- c("Tuberculosis","IllDefined","Respiratory","OtherDef.")
legend("topright",inset=c(0.005,0.04),leg=lg6,pch=ptc,cex=0.9,bg="ivory",y.intersp=0.8)

###########################

#COMPARING FOUR METHODS AND SAPARATE BETWEEN TRAINNING AND TESTING DATASET

###########################

#spliting data training and testing
library(caTools)
set.seed(101) 
sample = sample.split(m, SplitRatio = .70)
train = subset(m, sample == TRUE)
test  = subset(m, sample == FALSE)

# Logistic Regression model

mTr.lr <- glm(family=binomial,data=train,tb~factor(pro)+factor(agSx)+factor(drGh)) 
summary(mTr.lr)

prob.mTe.lr <- predict(mTr.lr,test,type="response")
pred.mTe.lr <- ifelse(prob.mTe.lr>0.5,1,0)
mTab <- table(pred.mTe.lr,test$tb)
addmargins(mTab)
pa.mTe.lr <- 100*sum(diag(mTab))/sum(mTab)
pa.mTe.lr

windows()
par(mfrow=c(2,1),oma=c(0,0,0,0),mar=c(3,3,3,4),
	las=1,mgp=c(1.1,0.1,0),tcl=0.2)
plot(prob.mTe.lr,test$tb,xlab="x coordinate",ylab="",ylim=c(-0.1,1.1),xlim=c(0,0.8))
mtext(side=3,adj=-0.1,line=0.2,"Tuberculosis")
mtext(side=3,adj=0.95,line=0.2,"9495 observation")
mtext(side=3,adj=1.07,line=0.2,"Total",font=3)
mtext(side=3,adj=0.5,line=1.4,"Logistic Regression Model",font=2)
testSum <- table(test$tb)
axis(side=4,at=c(0.2,0.8),lab=testSum)
abline(v=0.5,col=2,lwd=2)
text(c(0.2,0.8,0.2,0.8),c(0.2,0.2,0.8,0.8),mTab)
lg <- paste("Accuracy:",round(pa.mTe.lr,2),"%")
legend("right",inset=0.02,leg=lg,x.intersp=0.2,bg="ivory") 	#97.8

#################

# Rpart

#################

library(rpart)
library(partykit)
train$tb <- as.factor(train$tb)
mTr.rp <- rpart(tb~.,train[,c(2:4,9:38,39)],cp=0.0001)
mTr.rp2 <- prune(mTr.rp,cp=0.0012)

prob.mTe.rp <- predict(mTr.rp2,test,type="prob")
pred.mTe.rp <- ifelse(prob.mTe.rp[,2]>0.5,1,0)
mTab <- table(pred.mTe.rp,test$tb)
addmargins(mTab)
pa.mTe.rp <- 100*sum(diag(mTab))/sum(mTab)
pa.mTe.rp

xCoord <- prob.mTe.rp[,2]
yCoord <- as.integer(as.character(test$tb))

plot(xCoord,yCoord,xlab="x coordinate",ylab="",ylim=c(-0.1,1.1),xlim=c(0,0.8))
mtext(side=3,adj=-0.1,line=0.2,"Tuberculosis")
mtext(side=3,adj=0.95,line=0.2,"9495 observation")
mtext(side=3,adj=1.07,line=0.2,"Total",font=3)
mtext(side=3,adj=0.5,line=1.4,"Rpart",font=2)
testSum <- table(test$tb)
axis(side=4,at=c(0.2,0.8),lab=testSum)
abline(v=0.5,col=2,lwd=2)
text(c(0.2,0.8,0.2,0.8),c(0.2,0.2,0.8,0.8),mTab)
lg <- paste("Accuracy:",round(pa.mTe.rp,2),"%")
legend("right",inset=0.02,leg=lg,x.intersp=0.2,bg="ivory")
col <- as.integer(as.character(test$tb))				#97.88

#-----------------------------------------------
# ROC curves

glm(data=train,family=binomial,tb~factor(drGh)+factor(agSx)+factor(pro)) -> mod0
glm(data=train,family=binomial,tb~factor(drGh)) -> mod1
glm(data=train,family=binomial,tb~factor(agSx)) -> mod2
glm(data=train,family=binomial,tb~factor(pro)) -> mod3
glm(data=train,family=binomial,tb~factor(drGh)+factor(agSx)) -> mod4
glm(data=train,family=binomial,tb~factor(drGh)+factor(pro)) -> mod5
glm(data=train,family=binomial,tb~factor(agSx)+factor(pro)) -> mod6
glm(data=train,family=binomial,tb~1) -> mod7

glm(data=train,family=binomial,tb~factor(gp1)) -> mod8

trueYes <- train$tb

sens0 <- NULL
spec0 <- NULL
sens1 <- NULL
spec1 <- NULL
sens2 <- NULL
spec2 <- NULL
sens3 <- NULL
spec3 <- NULL
sens4 <- NULL
spec4 <- NULL
sens5 <- NULL
spec5 <- NULL
sens6 <- NULL
spec6 <- NULL
sens7 <- NULL
spec7 <- NULL
sens8 <- NULL
spec8 <- NULL

p10 <- NULL

for (i in c(1:10,10*(2:99)) ) {
  p1i <- i/1000
  predictYes7 <- ifelse(mod7$fit>=p1i,1,0)
  table(predictYes7,trueYes,useNA="always") -> tt
  sensi <- tt[2,2]/(tt[1,2]+tt[2,2])
  speci <- tt[1,1]/(tt[1,1]+tt[2,1])
  sens7 <- c(sens7,sensi)
  spec7 <- c(spec7,speci)
  predictYes0 <- ifelse(mod0$fit>=p1i,1,0)
  predictYes1 <- ifelse(mod1$fit>=p1i,1,0)
  predictYes2 <- ifelse(mod2$fit>=p1i,1,0)
  predictYes3 <- ifelse(mod3$fit>=p1i,1,0)
  predictYes4 <- ifelse(mod4$fit>=p1i,1,0)
  predictYes5 <- ifelse(mod5$fit>=p1i,1,0)
  predictYes6 <- ifelse(mod6$fit>=p1i,1,0)
  table(predictYes0,trueYes,useNA="always") -> tt
  sensi <- tt[2,2]/(tt[1,2]+tt[2,2])
  speci <- tt[1,1]/(tt[1,1]+tt[2,1])
  sens0 <- c(sens0,sensi)
  spec0 <- c(spec0,speci)
  table(predictYes1,trueYes,useNA="always") -> tt
  sensi <- tt[2,2]/(tt[1,2]+tt[2,2])
  speci <- tt[1,1]/(tt[1,1]+tt[2,1])
  sens1 <- c(sens1,sensi)
  spec1 <- c(spec1,speci)
  table(predictYes2,trueYes,useNA="always") -> tt
  sensi <- tt[2,2]/(tt[1,2]+tt[2,2])
  speci <- tt[1,1]/(tt[1,1]+tt[2,1])
  sens2 <- c(sens2,sensi)
  spec2 <- c(spec2,speci)
  table(predictYes3,trueYes,useNA="always") -> tt
  sensi <- tt[2,2]/(tt[1,2]+tt[2,2])
  speci <- tt[1,1]/(tt[1,1]+tt[2,1])
  sens3 <- c(sens3,sensi)
  spec3 <- c(spec3,speci)
  table(predictYes4,trueYes,useNA="always") -> tt
  sensi <- tt[2,2]/(tt[1,2]+tt[2,2])
  speci <- tt[1,1]/(tt[1,1]+tt[2,1])
  sens4 <- c(sens4,sensi)
  spec4 <- c(spec4,speci)
  table(predictYes5,trueYes,useNA="always") -> tt
  sensi <- tt[2,2]/(tt[1,2]+tt[2,2])
  speci <- tt[1,1]/(tt[1,1]+tt[2,1])
  sens5 <- c(sens5,sensi)
  spec5 <- c(spec5,speci)
  table(predictYes6,trueYes,useNA="always") -> tt
  sensi <- tt[2,2]/(tt[1,2]+tt[2,2])
  speci <- tt[1,1]/(tt[1,1]+tt[2,1])
  sens6 <- c(sens6,sensi)
  spec6 <- c(spec6,speci)
  predictYes8 <- ifelse(mod8$fit>=p1i,1,0)
  table(predictYes8,trueYes,useNA="always") -> tt
  sensi <- tt[2,2]/(tt[1,2]+tt[2,2])
  speci <- tt[1,1]/(tt[1,1]+tt[2,1])
  sens8 <- c(sens8,sensi)
  spec8 <- c(spec8,speci)
  p10 <- c(p10,p1i)
}
roc7 <- as.data.frame(cbind(p10,1-spec7,sens7))
roc7 <- subset(roc7,V2>0)
roc7 <- rbind(c(0,1,1),roc7,c(0,0,0))
roc1 <- as.data.frame(cbind(p10,1-spec1,sens1))
roc1 <- subset(roc1,V2>0)
roc1 <- rbind(c(0,1,1),roc1,c(0,0,0))
roc2 <- as.data.frame(cbind(p10,1-spec2,sens2))
roc2 <- subset(roc2,V2>0)
roc2 <- rbind(c(0,1,1),roc2,c(0,0,0))
roc3 <- as.data.frame(cbind(p10,1-spec3,sens3))
roc3 <- subset(roc3,V2>0)
roc3 <- rbind(c(0,1,1),roc3,c(0,0,0))
roc4 <- as.data.frame(cbind(p10,1-spec4,sens4))
roc4 <- subset(roc4,V2>0)
roc4 <- rbind(c(0,1,1),roc4,c(0,0,0))
roc5 <- as.data.frame(cbind(p10,1-spec5,sens5))
roc5 <- subset(roc5,V2>0)
roc5 <- rbind(c(0,1,1),roc5,c(0,0,0))
roc6 <- as.data.frame(cbind(p10,1-spec6,sens6))
roc6 <- subset(roc6,V2>0)
roc6 <- rbind(c(0,1,1),roc6,c(0,0,0))
roc8 <- as.data.frame(cbind(p10,1-spec8,sens8))
roc8 <- subset(roc8,V2>0)
roc8 <- rbind(c(0,1,1),roc8,c(0,0,0))
roc0 <- as.data.frame(cbind(p10,1-spec0,sens0))
roc0 <- subset(roc0,V2>0)
roc0 <- rbind(c(0,1,1),roc0,c(0,0,0))

windows(5,5)
par(mar=c(2.5,0.5,0,0.6),mgp=c(1.1,0.2,0),oma=c(0,1.8,2.5,0),las=1,tcl=-0.2)

plot(roc7[,2],roc7[,3],type="l",ylim=c(0,1.05),xlim=c(0,1),yaxs="i",
     xlab="False positive rate",ylab="",cex.lab=1.1,lwd=2,col=8)
npts <- nrow(roc8)
roc8r <- roc8[npts:1,]
polygon(c(roc0[,2],roc8r[,2]),c(roc0[,3],roc8r[,3]),col="grey")

points(c(0,1),c(0,1),type="l",col=8)

clrs <- c("pink",3,5,2,4,6,8)
points(roc1[,2],roc1[,3],type="l",col=clrs[1],lwd=2)
points(roc2[,2],roc2[,3],type="l",col=clrs[2],lwd=2)
points(roc3[,2],roc3[,3],type="l",col=clrs[3],lwd=2)
points(roc4[,2],roc4[,3],type="l",col=clrs[4],lwd=2)
points(roc5[,2],roc5[,3],type="l",col=clrs[5],lwd=2)
points(roc6[,2],roc6[,3],type="l",col=clrs[6],lwd=2)
points(roc8[,2],roc8[,3],type="l",col="brown",lwd=3)
points(roc0[,2],roc0[,3],type="l",lwd=2)
points(c(0,0,1),c(0,1,1),type="l")

mtext(side=3,adj=-0.08,line=0.2,"Sensitivity")
mtext(side=3,adj=0.5,line=1.3,"ROC Curves: Logistic Models")
mtext(side=3,adj=1,line=0.2,"TB Death Prevalence")

lg0 <- c("drGh+ageSex+Prov","drGh only","ageSex only","Province only",
         "drGh+ageSex","drGh+Prov","ageSex+Prov","Constant","Tree-based method")
legend("bottomright",lg0,inset=c(0.01,0.23),x.intersp=0.4,lwd=c(2,2,2,2,2,2,2,2,3),
       col=c(1,clrs,"brown"),bg="ivory")

Area <- function(X) {		# function to compute area of a polygon
  X <- rbind(X,X[1,])
  x <- X[,1]
  y <- X[,2]
  lx <- length(x)
  sum((x[2:lx]-x[1:lx-1])*(y[2:lx]+y[1:lx-1]))/2
}

xy0 <- rbind(roc0[,c(2,3)],c(1,0),c(1,1))
area0 <- 2*(-Area(xy0)-0.5)
xy8 <- rbind(roc8[,c(2,3)],c(1,0),c(1,1))
area8 <- 2*(-Area(xy8)-0.5)
titl <- "Area under curve above y=x"
lg1 <- paste("Logistic Regression: ",round(100*area0,2),"%",sep="")
lg2 <- paste("Tree-based Method: ",round(100*area8,2),"%",sep="")
lg <- c(lg1,lg2)
legend("bottomright",title=titl,inset=c(0.01,0.01),lg,lwd=2,
       col=c(1,"brown"),bg="ivory",x.intersp=0.1)

#-----------------------------------------end

