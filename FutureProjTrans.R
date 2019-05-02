###Projection of transition rates under future climate conditions

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

##CRCRPIPO
C.JMx.30 <- 2.83
C.JMx.S12 <- 3.13 #3 degrees hotter
C.JMx.S34 <- 3.43 #6 degrees hotter

C.S.30 <- 0.28
C.S.S12 <- 0.17 #40% lower snow
C.S.S34 <- 0.34 #20% higher snow

##EMRIDGE
E.JMx.30 <- 1.73
E.JMx.S12 <- 2.03 #3 degrees hotter
E.JMx.S34 <- 2.33 #6 degrees hotter

E.S.30 <- 1.17
E.S.S12 <- 0.7 #40% lower snow
E.S.S34 <- 1.4 #20% higher snow
###Species
Spp <- c("ABCO","ABMA","CADE","PILA","PIMO","PIPO","QU","SEGI")
nspp <- length(Spp)

#####Logit function
inv.logit <- function(x){
	a <- exp(x)/(1+exp(x))
}

#######Transition, small/lg seedlings, 0.2 m2 BA, no fire
IS.beta <- c(-3.41,6.29,2.98,3.23) #intercept, size
Spp.beta <- c(2.66,-2.49,1.67,-1.8,0.35,-0.76,0.06,1.05) #ABCO,ABMA,CADE,PILA,PIMO,PIPO,QU,SEGI
fire.beta <- c(1.62,0.01,-1.02,0.52,1.07,0.05,-0.71,1.19,0.51)
Jmx.Beta<-c(-1.83,0.04,-1.35,-0.17,-0.66,-0.1,-1.13,-0.58)# ABCO,ABMA,CADE,PILA,PIMO,PIPO,QU,SEGI;
S.Beta<-c(-0.35,-0.19,0.66,0.3,-1.27,-1.09,1.1,-1.37)# ABCO,ABMA,CADE,PILA,PIMO,PIPO,QU,SEGI;
BA.Beta<-c(-0.04,0.04,-0.14,-0.05,-0.29,0.07,-0.02,-0.11) 
# BA:ABCO,ABMA,CADE,PILA,PIMO,PIPO,QU,SEGI

############Some species

### CRCRPIPO, Small seedlings, 30 yr av., 4 change scenarios
Trans.C.S <- matrix(0,5,nspp)
 for(s in 1:nspp){
 Trans.C.S[1,s]<- inv.logit(IS.beta[1]+IS.beta[2]+Spp.beta[s]+Jmx.Beta[s]*C.JMx.30+S.Beta[s]*C.S.30+BA.Beta[s]*0.2)
 Trans.C.S[2,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+Jmx.Beta[s]*C.JMx.S12 +S.Beta[s]*C.S.S12 +BA.Beta[s]*0.2) 
  Trans.C.S[3,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+Jmx.Beta[s]*C.JMx.S12+S.Beta[s]*C.S.S34+BA.Beta[s]*0.2)
  Trans.C.S[4,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+Jmx.Beta[s]*C.JMx.S34+S.Beta[s]*C.S.S34+BA.Beta[s]*0.2)
  Trans.C.S[5,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+Jmx.Beta[s]*C.JMx.S34+ S.Beta[s]*C.S.S12+ BA.Beta[s]*0.2)   
 } 

rownames(Trans.C.S)<-c("base","3C.LS","3C.MS","6C.LS","6C.MS")
colnames(Trans.C.S)<-Spp

### CRCRPIPO, Large seedlings, 30 yr av., 4 change scenarios
Trans.C.L <- matrix(0,5,nspp)
 for(s in 1:nspp){
 Trans.C.L[1,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+Jmx.Beta[s]*C.JMx.30+S.Beta[s]*C.S.30+BA.Beta[s]*0.2)
 Trans.C.L[2,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+Jmx.Beta[s]*C.JMx.S12+ S.Beta[s]*C.S.S12+ BA.Beta[s]*0.2) 
  Trans.C.L[3,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+Jmx.Beta[s]*C.JMx.S12+ S.Beta[s]*C.S.S34 +BA.Beta[s]*0.2)
  Trans.C.L[4,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+Jmx.Beta[s]*C.JMx.S34+ S.Beta[s]*C.S.S12 +BA.Beta[s]*0.2)
  Trans.C.L[5,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+Jmx.Beta[s]*C.JMx.S34+ S.Beta[s]*C.S.S34 +BA.Beta[s]*0.2)  
 } 

rownames(Trans.C.L)<-c("base","3C.LS","3C.MS","6C.LS","6C.MS")
colnames(Trans.C.L)<-Spp

### EMRIDGE, Small seedlings, 30 yr av., 4 change scenarios
Trans.E.S <- matrix(0,5,nspp)
 for(s in 1:nspp){
 Trans.E.S[1,s]<- inv.logit(IS.beta[1]+IS.beta[2]+Spp.beta[s]+Jmx.Beta[s]*E.JMx.30+S.Beta[s]*E.S.30+BA.Beta[s]*0.2)
 Trans.E.S[2,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+Jmx.Beta[s]*E.JMx.S12+S.Beta[s]*E.S.S12 +BA.Beta[s]*0.2) 
 Trans.E.S[3,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+Jmx.Beta[s]*E.JMx.S12+S.Beta[s]*E.S.S34 +BA.Beta[s]*0.2)
 Trans.E.S[4,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+Jmx.Beta[s]*E.JMx.S34+S.Beta[s]*E.S.S12 +BA.Beta[s]*0.2)
 Trans.E.S[5,s]<- inv.logit(IS.beta[1]+IS.beta[2] +Spp.beta[s]+Jmx.Beta[s]*E.JMx.S34+S.Beta[s]*E.S.S34 +BA.Beta[s]*0.2)   
 } 

rownames(Trans.E.S)<-c("base","3C.LS","3C.MS","6C.LS","6C.MS")
colnames(Trans.E.S)<-Spp

### EMRIDGE, Large seedlings, 30 yr av., 4 change scenarios
Trans.E.L <- matrix(0,5,nspp)
 for(s in 1:nspp){
  Trans.E.L[1,s]<- inv.logit(IS.beta[1]+IS.beta[4]+Spp.beta[s]+Jmx.Beta[s]*E.JMx.30+S.Beta[s]*E.S.30+BA.Beta[s]*0.2)
 Trans.E.L[2,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+Jmx.Beta[s]*E.JMx.S12+S.Beta[s]*E.S.S12 +BA.Beta[s]*0.2) 
  Trans.E.L[3,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+Jmx.Beta[s]*E.JMx.S12+S.Beta[s]*E.S.S34 +BA.Beta[s]*0.2)
  Trans.E.L[4,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+Jmx.Beta[s]*E.JMx.S34+S.Beta[s]*E.S.S12 +BA.Beta[s]*0.2)
  Trans.E.L[5,s]<- inv.logit(IS.beta[1]+IS.beta[4] +Spp.beta[s]+Jmx.Beta[s]*E.JMx.S34+S.Beta[s]*E.S.S34 +BA.Beta[s]*0.2)  
 } 

rownames(Trans.E.L)<-c("base","3C.LS","3C.MS","6C.LS","6C.MS")
colnames(Trans.E.L)<-Spp

Trans.C.S <- Trans.C.S[,c(1,3,4,6,7)]
Trans.C.L <- Trans.C.L[,c(1,3,4,6,7)]
Trans.E.S <- Trans.E.S[,c(2,5)]
Trans.E.L <- Trans.E.L[,c(2,5)]


tiff("Figure5.tiff", width = 6, height = 6, units = 'in', res = 300, compression = 'lzw')
par(mfrow=c(2,2),las=2,cex=0.7)
barplot(Trans.C.S,cex.names=0.9,cex.axis=0.9,ylab="transition probability",main="Small sdls, Low elevation site",ylim=c(0,0.9),beside=TRUE)

barplot(Trans.C.L,cex.names=0.9,cex.axis=0.9,ylab="transition probability",main="Large sdls, Low elevation site",ylim=c(0,0.5),beside=TRUE,legend=rownames(Trans.C.L))

barplot(Trans.E.S,cex.names=0.9,cex.axis=0.9,ylab="transition probability",main="Small sdls, High elevation site",ylim=c(0,0.9),beside=TRUE)

barplot(Trans.E.L,cex.names=0.9,cex.axis=0.9,ylab="transition probability",main="Large sdls, High elevation site",ylim=c(0,0.5),beside=TRUE)
dev.off()