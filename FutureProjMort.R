###Projection of mortality rates under future climate conditions

allplots <- c('YOHOPIPO','BBBPIPO','CCRPIPO','CRCRPIPO','FFS7CONTROL','FFS6BURN','FFS5BURN','SU','FRPIJE','FFS2BURN','LMCC','LOTHAR','LOGSEGI','UPTHAR','LOLOG','UPLOG','LOGPIJE','SFTRABMA','WTABMA','POFLABMA','PGABMA','EMSLOPE','EMRIDGE')

nplots <- length(allplots)

years <- seq(1999,2008,by=1)

#####Logit function
inv.logit <- function(x){
	a <- exp(x)/(1+exp(x))
}

######Long-term climate data
bframe6 <- data.frame(read.table("climate_data_long.txt",header=T))
basin <- c(27,0,1,3,10,9,8,21,11,7,13,17,15,25,16,24,14,20,26,19,18,6,4) 


##Average climate by plot
est.yrs <- seq(1976,2013,by=1)
for(p in 1:nplots){
	q <- which(bframe6[,3]==basin[p])
	if(length(q)>0){
	clim.temp <- bframe6[q,]
	Plt <- rep(p,length(est.yrs))
	precip.all <- numeric(0); JulMax.all <- numeric(0)
	JanMin.all <- numeric(0); AvTemp.all <- numeric(0)
	snow.all <- numeric(0); CWD.all <- numeric(0)
	for(t in est.yrs){
		precip.all <- c(precip.all,sum(clim.temp[which(clim.temp[,1]==t),4]))
		snow.all <- c(snow.all,sum(clim.temp[which(clim.temp[,1]==t),9]))
		CWD.all <- c(CWD.all,sum(clim.temp[which(clim.temp[,1]==t),15]))
		JulMax.all <- c(JulMax.all,clim.temp[which(clim.temp[,1]==t & clim.temp[,2]==10),6])
		JanMin.all <- c(JanMin.all,clim.temp[which(clim.temp[,1]==t & clim.temp[,2]==4),7])
		AvTemp.all <- c(AvTemp.all,mean(clim.temp[which(clim.temp[,1]==t),8]))
	}
	
	plot.clim <- cbind(precip.all,JulMax.all,JanMin.all,AvTemp.all,snow.all,CWD.all)
	mn.plot.clim <- apply(plot.clim,2,mean)
	
	if(p==1) {
		Climate.all <- cbind(Plt,est.yrs,plot.clim)
		Climate.mean <- c(p,mn.plot.clim)
		}
	if(p>1) {
		Climate.all <- rbind(Climate.all,cbind(Plt,est.yrs,plot.clim))
		Climate.mean <- rbind(Climate.mean,c(p,mn.plot.clim))
		}
	}
}

rownames(Climate.mean) <- allplots

###Species-level mean climatic variables
ABCO.JMx <- mean(Climate.mean[c(1:17),3])/10
ABMA.JMx <- mean(Climate.mean[c(8,10:11,13,15:21,23),3])/10
CADE.JMx <- mean(Climate.mean[c(1:10,15:16),3])/10
PILA.JMx <- mean(Climate.mean[c(1:17),3])/10
PIMO.JMx <- mean(Climate.mean[c(9,19,22:23),3])/10
PIPO.JMx <- mean(Climate.mean[c(1:4,6),3])/10
QU.JMx <- mean(Climate.mean[c(1:4,6,8:9,16),3])/10
SEGI.JMx <- mean(Climate.mean[c(11,13,15),3])/10

Spp.JMx <- c(ABCO.JMx,ABMA.JMx,CADE.JMx,PILA.JMx,PIMO.JMx,PIPO.JMx,QU.JMx,SEGI.JMx)

ABCO.P <- mean(Climate.mean[c(1:17),1])/1000
ABMA.P <- mean(Climate.mean[c(8,10:11,13,15:21,23),1])/1000
CADE.P <- mean(Climate.mean[c(1:10,15:16),1])/1000
PILA.P <- mean(Climate.mean[c(1:17),1])/1000
PIMO.P <- mean(Climate.mean[c(9,19,22:23),1])/1000
PIPO.P <- mean(Climate.mean[c(1:4,6),1])/1000
QU.P <- mean(Climate.mean[c(1:4,6,8:9,16),1])/1000
SEGI.P <- mean(Climate.mean[c(11,13,15),1])/1000

Spp.P <- c(ABCO.P,ABMA.P,CADE.P,PILA.P,PIMO.P,PIPO.P,QU.P,SEGI.P)

##CRCRPIPO
C.JMx.30 <- 2.83
C.JMx.S12 <- 3.13 #3 degrees hotter
C.JMx.S34 <- 3.43 #6 degrees hotter
C.JMxD.S12 <- C.JMx.S12-C.JMx.30
C.JMxD.S34 <- C.JMx.S34-C.JMx.30

C.P.30 <- 0.992
C.P.S13 <- 0.794 #20% lower precip
C.P.S24 <- 1.389 #40% higher precip
C.PD.S13 <- C.P.S13-C.P.30
C.PD.S24 <- C.P.S24-C.P.30

##EMRIDGE
E.JMx.30 <- 1.73
E.JMx.S12 <- 2.03 #3 degrees hotter
E.JMx.S34 <- 2.33 #6 degrees hotter
E.JMxD.S12 <- E.JMx.S12-E.JMx.30
E.JMxD.S34 <- E.JMx.S34-E.JMx.30

E.P.30 <- 1.356
E.P.S13 <- 1.085 #20% lower precip
E.P.S24 <- 1.898 #40% higher precip
E.PD.S13 <- E.P.S13-E.P.30
E.PD.S24 <- E.P.S24-E.P.30

###Species
Spp <- c("ABCO","ABMA","CADE","PILA","PIMO","PIPO","QU","SEGI")
nspp <- length(Spp)

#####Logit function
inv.logit <- function(x){
	a <- exp(x)/(1+exp(x))
}

#######Survival, small/lg seedlings, 0.2 m2 BA, no fire
IS.beta <- c(4.09,-2.3,-0.83,0.27) #intercept, size
Spp.beta <- c(0.19,-0.24,0.46,0.09,2.51,-1.08,-0.8,0.74) #ABCO,ABMA,CADE,PILA,PIMO,PIPO,QU,SEGI
fire.beta <- c(-1.15,-0.91,-0.05,-4.35,-0.98,-0.75,-5.7,-3.47,-0.16)
MnJmx.D.Beta<-c(-0.37,-0.12,-0.63,-0.34,-0.83,0.45,-0.42,-0.14)# ABCO,ABMA,CADE,PILA,PIMO,PIPO,QU,SEGI;
P.D.Beta<-c(0.13,0.03,0.13,0.15,0.37,-0.2,0.23,0.34) #P.D: ABCO,ABMA,CADE,PILA,PIMO,PIPO,QU,SEGI;
BA.Beta<- -0.03

############All species

### CRCRPIPO, Small seedlings, 30 yr av., 4 change scenarios
Surv.C.S <- matrix(0,5,nspp)
 for(s in 1:nspp){
  Surv.C.S[1,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+MnJmx.D.Beta[s]*(C.JMx.30-Spp.JMx[s])+P.D.Beta[s]*(C.P.30-Spp.P[s])+BA.Beta*0.2)
  Surv.C.S[2,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+MnJmx.D.Beta[s]*(C.JMx.S12-Spp.JMx[s])+P.D.Beta[s]*(C.P.S13-Spp.P[s])+BA.Beta*0.2)
  Surv.C.S[3,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+MnJmx.D.Beta[s]*(C.JMx.S12-Spp.JMx[s])+P.D.Beta[s]*(C.P.S24-Spp.P[s])+BA.Beta*0.2)
  Surv.C.S[4,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+MnJmx.D.Beta[s]*(C.JMx.S34-Spp.JMx[s])+P.D.Beta[s]*(C.P.S13-Spp.P[s])+BA.Beta*0.2)
  Surv.C.S[5,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+MnJmx.D.Beta[s]*(C.JMx.S34-Spp.JMx[s])+P.D.Beta[s]*(C.P.S24-Spp.P[s])+BA.Beta*0.2)  
 } 

rownames(Surv.C.S)<-c("base","3C.W","3C.D","6C.W","6C.D")
colnames(Surv.C.S)<-Spp

Mort.C.S <- 1-Surv.C.S

### CRCRPIPO, Large seedlings, 30 yr av., 4 change scenarios
Surv.C.L <- matrix(0,5,nspp)
 for(s in 1:nspp){
    Surv.C.L[1,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+MnJmx.D.Beta[s]*(C.JMx.30-Spp.JMx[s])+P.D.Beta[s]*(C.P.30-Spp.P[s])+BA.Beta*0.2)
  Surv.C.L[2,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+MnJmx.D.Beta[s]*(C.JMx.S12-Spp.JMx[s])+P.D.Beta[s]*(C.P.S13-Spp.P[s])+BA.Beta*0.2)
  Surv.C.L[3,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+MnJmx.D.Beta[s]*(C.JMx.S12-Spp.JMx[s])+P.D.Beta[s]*(C.P.S24-Spp.P[s])+BA.Beta*0.2)
  Surv.C.L[4,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+MnJmx.D.Beta[s]*(C.JMx.S34-Spp.JMx[s])+P.D.Beta[s]*(C.P.S13-Spp.P[s])+BA.Beta*0.2)
  Surv.C.L[5,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+MnJmx.D.Beta[s]*(C.JMx.S34-Spp.JMx[s])+P.D.Beta[s]*(C.P.S24-Spp.P[s])+BA.Beta*0.2) 
 } 

rownames(Surv.C.L)<-c("base","3C.W","3C.D","6C.W","6C.D")
colnames(Surv.C.L)<-Spp

Mort.C.L <- 1-Surv.C.L

### EMRIDGE, Small seedlings, 30 yr av., 4 change scenarios
Surv.E.S <- matrix(0,5,nspp)
 for(s in 1:nspp){
  Surv.E.S[1,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+MnJmx.D.Beta[s]*(E.JMx.30-Spp.JMx[s])+P.D.Beta[s]*(E.P.30-Spp.P[s])+BA.Beta*0.2)
  Surv.E.S[2,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+MnJmx.D.Beta[s]*(E.JMx.S12-Spp.JMx[s])+P.D.Beta[s]*(E.P.S13-Spp.P[s])+BA.Beta*0.2) 
  Surv.E.S[3,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+MnJmx.D.Beta[s]*(E.JMx.S12-Spp.JMx[s])+P.D.Beta[s]*(E.P.S24-Spp.P[s])+BA.Beta*0.2)
  Surv.E.S[4,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+MnJmx.D.Beta[s]*(E.JMx.S34-Spp.JMx[s])+P.D.Beta[s]*(E.P.S13-Spp.P[s])+BA.Beta*0.2)
  Surv.E.S[5,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+MnJmx.D.Beta[s]*(E.JMx.S34-Spp.JMx[s])+P.D.Beta[s]*(E.P.S24-Spp.P[s])+BA.Beta*0.2)  
 } 

rownames(Surv.E.S)<-c("base","3C.W","3C.D","6C.W","6C.D")
colnames(Surv.E.S)<-Spp

Mort.E.S <- 1-Surv.E.S

### EMRIDGE, Large seedlings, 30 yr av., 4 change scenarios
Surv.E.L <- matrix(0,5,nspp)
 for(s in 1:nspp){
    Surv.E.L[1,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+MnJmx.D.Beta[s]*(E.JMx.30-Spp.JMx[s])+P.D.Beta[s]*(E.P.30-Spp.P[s])+BA.Beta*0.2)
  Surv.E.L[2,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+MnJmx.D.Beta[s]*(E.JMx.S12-Spp.JMx[s])+P.D.Beta[s]*(E.P.S13-Spp.P[s])+BA.Beta*0.2)
  Surv.E.L[3,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+MnJmx.D.Beta[s]*(E.JMx.S12-Spp.JMx[s])+P.D.Beta[s]*(E.P.S24-Spp.P[s])+BA.Beta*0.2)
  Surv.E.L[4,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+MnJmx.D.Beta[s]*(E.JMx.S34-Spp.JMx[s])+P.D.Beta[s]*(E.P.S13-Spp.P[s])+BA.Beta*0.2)
  Surv.E.L[5,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+MnJmx.D.Beta[s]*(E.JMx.S34-Spp.JMx[s])+P.D.Beta[s]*(E.P.S24-Spp.P[s])+BA.Beta*0.2)
 } 

rownames(Surv.E.L)<-c("base","3C.W","3C.D","6C.W","6C.D")
colnames(Surv.E.L)<-Spp

Mort.E.L <- 1-Surv.E.L

Mort.C.S <- Mort.C.S[,c(1,3,4,6,7)]
Mort.C.L <- Mort.C.L[,c(1,3,4,6,7)]
Mort.E.S <- Mort.E.S[,c(2,5)]
Mort.E.L <- Mort.E.L[,c(2,5)]

tiff("Figure4.tiff", width = 6, height = 6, units = 'in', res = 300, compression = 'lzw')
par(mfrow=c(2,2),las=2)
barplot(Mort.C.S,cex.names=0.8,cex.axis=0.9,ylab="mortality probability",main="Small sdls, Low elevation site",ylim=c(0,0.45),beside=TRUE)

barplot(Mort.C.L,cex.names=0.8,cex.axis=0.9,ylab="mortality probability",main="Large sdls, Low elevation site",ylim=c(0,0.2),beside=TRUE,legend=rownames(Mort.C.L),args.legend=list(cex=0.8))

barplot(Mort.E.S,cex.names=0.8,cex.axis=0.9,ylab="mortality probability",main="Small sdls, High elevation site",ylim=c(0,0.45),beside=TRUE)

barplot(Mort.E.L,cex.names=0.8,cex.axis=0.9,ylab="mortality probability",main="Large sdls, High elevation site",ylim=c(0,0.2),beside=TRUE)
dev.off()
