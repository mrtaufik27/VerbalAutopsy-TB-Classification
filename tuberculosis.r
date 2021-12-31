# kanittaTB.Rcm

library(MASS)
library(grid)
library(rpart)
library(partykit)
library(epiDisplay)

setwd("C:/Users/USER/Dropbox/project/research in practice and symposium") 
read.table("va.txt",h=T,as.is=T) -> m
str(m)

#data management 

m$age <- ifelse(is.na(m$age),99,m$age)          # replace 2 unknown ages by 99
m <- subset(m,m$age >= 5)                       # age at least 5 

m$VA <- ifelse(m$VA=="O758","N758",m$VA)        # correct erroneous VA
m$VA <- ifelse(m$VA=="O758NA","N758NA",m$VA)    # correct erroneous VA

cgs <- c("TB","se","HIV","oID",		              # VR cause groups
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

#m$tb <- as.factor(ifelse(m$VAgrp=="01:a","Tuberculosis","other"))

ag <- m$age; sx <- m$sex;
m$agSx <- ifelse(ag<40 & sx==1,1,ifelse(ag<40 & sx==2,2,
                                        ifelse(ag<50 & sx==1,3,ifelse(ag<50 & sx==2,4,
                                                                      ifelse(ag<60 & sx==1,5,ifelse(ag<60 & sx==2,6,
                                                                                                    ifelse(ag<70 & sx==1,7,ifelse(ag<70 & sx==2,8,
                                                                                                                                  ifelse(ag<80 & sx==1,9,ifelse(ag<80 & sx==2,10,
                                                                                                                                                                ifelse(ag>=80 & sx==1,11,ifelse(ag>=80 & sx==2,12,13))))) ))))) ))

m$tb <- ifelse(m$VAgrp=="01:a",1,0)
m$pro <- as.factor(as.character(m$pro))
m$ghos <- as.factor(m$ghos)
m$agSx <- as.factor(m$agSx)
m$DRgrp <- as.factor(m$DRgrp)

summary(m)
names(m)

#fix dataframe
dat <- m[c(1,2,8,9,10)]

#tree based method
names(m)
m.rp <- rpart(tb~.,m,cp=0.0001)
printcp(m.rp)

windows(10,6)
par(las=1,mar=c(4,4,4,2))
plotcp(m.rp)

windows(15,7)
#par(mar=c(1,0,1,0))
plot(m.rp,uniform=T); text(m.rp,use.n=T,cex=1,digits=3,font=5)
#m$gp0 <- predict(as.party(m.rp),type="node")	

#glm(family=binomial,data=m,tb~factor(gp0)) -> mod0
#summary(mod0)
#lroc(mod0)$auc

#pruning the tree
m.rp1 <- prune(m.rp,cp=0.0007)
windows(10,6)
par(las=1,mar=c(4,4,4,2))
plotcp(m.rp)

windows(12,8)
par(mar=c(1,1,1,1))
plot(m.rp1,uniform=T); text(m.rp1,use.n=T,cex=0.8,font=2)

m.rp2 <- prune(m.rp,cp=0.0012)
windows(12,8)
par(mar=c(1,1,1,1))
plot(m.rp2,uniform=T); text(m.rp2,use.n=T,cex=0.8,font=2)

windows(10,6)
par(las=1,mar=c(4,4,4,2))
plotcp(m.rp2)
printcp(m.rp2)

m$gp <- predict(as.party(m.rp2),type="node")	

addmargins(table(m$gp))

glm(family=binomial,data=m,tb~factor(gp)) -> mod
summary(mod)
lroc(mod)$auc

#-------------------------------------------------------------------------#
# plot confidence intervals 

source("C:/Users/ACER/Dropbox/project/research in practice/dcis.Rcm")

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
axis(side=1,at=at,lab=xlab,tcl=0,padj=1.4)

mtext(side=3,line=0.1,adj=-0.1,ylab)
mtext(side=3,line=0.1,adj=1,titl)

#--------------------------------------------
# crude and adjusted OR
# plot unadjusted means

pcCr <- 100*tapply(m$tb,m$gp1,mean)	#crude

points(xCoord-0.2,10*sqrt(pcCr),pch=21,bg=3,cex=0.8)

axis(side=1,at=c(1:n1),lab=c(1:n))
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

#COMPARING TWO METHODS AND SAPARATING BETWEEN TRAINNING AND TESTING DATASET

###########################

#spliting data training and testing
library(caTools)
set.seed(101) 
sample = sample.split(m, SplitRatio = .70)
train = subset(m, sample == TRUE)
test  = subset(m, sample == FALSE)

# Logistic Regression model

mTr.lr <- glm(family=binomial,data=mTr,tb~factor(pro)+factor(agSx)+factor(drGh)) 
summary(mTr.lr)

#mTr.lr <- glm(data=mTr[,c(2:4,9:38,39)],family=binomial,tb~.)
#summary(mTr.lr)

prob.mTe.lr <- predict(mTr.lr,mTe,type="response")
pred.mTe.lr <- ifelse(prob.mTe.lr>0.5,1,0)
mTab <- table(pred.mTe.lr,mTe$tb)
addmargins(mTab)
pa.mTe.lr <- 100*sum(diag(mTab))/sum(mTab)
pa.mTe.lr

windows()
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(3,3,3,4),
    las=1,mgp=c(1.1,0.1,0),tcl=0.2)
plot(prob.mTe.lr,mTe$tb,xlab="x coordinate",ylab="",ylim=c(-0.1,1.1),xlim=c(0,0.8))
mtext(side=3,adj=-0.1,line=0.2,"Tuberculosis")
mtext(side=3,adj=0.95,line=0.2,"9495 observation")
mtext(side=3,adj=1.07,line=0.2,"Total",font=3)
mtext(side=3,adj=0.5,line=1.4,"Logistic Regression Model",font=2)
testSum <- table(mTe$tb)
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
mTr$tb <- as.factor(mTr$tb)
mTr.rp <- rpart(tb~.,mTr[,c(2:4,9:38,39)],cp=0.0001)
mTr.rp2 <- prune(mTr.rp,cp=0.0012)

prob.mTe.rp <- predict(mTr.rp2,mTe,type="prob")
pred.mTe.rp <- ifelse(prob.mTe.rp[,2]>0.5,1,0)
mTab <- table(pred.mTe.rp,mTe$tb)
addmargins(mTab)
pa.mTe.rp <- 100*sum(diag(mTab))/sum(mTab)
pa.mTe.rp

xCoord <- prob.mTe.rp[,2]
yCoord <- as.integer(as.character(mTe$tb))

plot(xCoord,yCoord,xlab="x coordinate",ylab="",ylim=c(-0.1,1.1),xlim=c(0,0.8))
mtext(side=3,adj=-0.1,line=0.2,"Tuberculosis")
mtext(side=3,adj=0.95,line=0.2,"9495 observation")
mtext(side=3,adj=1.07,line=0.2,"Total",font=3)
mtext(side=3,adj=0.5,line=1.4,"Rpart",font=2)
testSum <- table(mTe$tb)
axis(side=4,at=c(0.2,0.8),lab=testSum)
abline(v=0.5,col=2,lwd=2)
text(c(0.2,0.8,0.2,0.8),c(0.2,0.2,0.8,0.8),mTab)
lg <- paste("Accuracy:",round(pa.mTe.rp,2),"%")
legend("right",inset=0.02,leg=lg,x.intersp=0.2,bg="ivory")
col <- as.integer(as.character(mTe$tb))				#97.88


# Logistic Regression model
mTr.lr <- glm(data=mTr[,c(2:4,9:38,39)],family=binomial,tb~.)
summary(mTr.lr)

prob.mTe.lr <- predict(mTr.lr,mTe,type="response")
pred.mTe.lr <- ifelse(prob.mTe.lr>0.5,1,0)
mTab <- table(pred.mTe.lr,mTe$tb)
addmargins(mTab)
pa.mTe.lr <- 100*sum(diag(mTab))/sum(mTab)
pa.mTe.lr

windows()
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(3,3,3,4),
	las=1,mgp=c(1.1,0.1,0),tcl=0.2)
plot(prob.mTe.lr,mTe$tb,ylab="",ylim=c(-0.1,1.1))
mtext(side=3,adj=-0.1,line=0.2,"Tuberculosis")
mtext(side=3,adj=0.95,line=0.2,"9495 observation")
mtext(side=3,adj=1.07,line=0.2,"Total",font=3)
mtext(side=3,adj=0.5,line=1.4,"Logistic Regression Model",font=2)
testSum <- table(mTe$tb)
axis(side=4,at=c(0.2,0.8),lab=testSum)
abline(v=0.5,col=2,lwd=2)
text(c(0.2,0.8,0.2,0.8),c(0.2,0.2,0.8,0.8),mTab)
lg <- paste("Accuracy:",round(pa.mTe.lr,2),"%")
legend("right",inset=0.02,leg=lg,x.intersp=0.2,bg="ivory") 	#97.69

#################

# Rpart

#################

library(rpart)
library(partykit)
mTr$tb <- as.factor(mTr$tb)
mTr.rp <- rpart(tb~.,mTr[,c(2:4,9:38,39)],cp=0.0001)
mTr.rp2 <- prune(mTr.rp,cp=0.0012)

prob.mTe.rp <- predict(mTr.rp2,mTe,type="prob")
pred.mTe.rp <- ifelse(prob.mTe.rp[,2]>0.5,1,0)
mTab <- table(pred.mTe.rp,mTe$tb)
addmargins(mTab)
pa.mTe.rp <- 100*sum(diag(mTab))/sum(mTab)
pa.mTe.rp

xCoord <- prob.mTe.rp[,2]
yCoord <- as.integer(as.character(mTe$tb))

plot(xCoord,yCoord,ylab="",ylim=c(-0.1,1.1))
mtext(side=3,adj=-0.1,line=0.2,"Tuberculosis")
mtext(side=3,adj=0.95,line=0.2,"9495 observation")
mtext(side=3,adj=1.07,line=0.2,"Total",font=3)
mtext(side=3,adj=0.5,line=1.4,"Rpart",font=2)
testSum <- table(mTe$tb)
axis(side=4,at=c(0.2,0.8),lab=testSum)
abline(v=0.5,col=2,lwd=2)
text(c(0.2,0.8,0.2,0.8),c(0.2,0.2,0.8,0.8),mTab)
lg <- paste("Accuracy:",round(pa.mTe.rp,2),"%")
legend("right",inset=0.02,leg=lg,x.intersp=0.2,bg="ivory")
col <- as.integer(as.character(mTe$tb))				#97.64

