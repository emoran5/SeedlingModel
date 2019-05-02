########## Gets Data Organized To Go Into Models

bframe1 <-data.frame(read.table("TaggedSdl.txt",header=T))

######################### Tree plots
#####Standard plots
allplots <- c('YOHOPIPO','CCRPIPO','CRCRPIPO','FFS7CONTROL','FFS6BURN','FFS5BURN','FFS2BURN','LOTHAR','UPTHAR','LOLOG','UPLOG','LOGPIJE','SFTRABMA','WTABMA','POFLABMA','PGABMA','BBBPIPO','FRPIJE','LMCC','LOGSEGI','EMSLOPE','EMRIDGE') 
nplots <- length(allplots)

####"Extra" plots (the 3 that are next to each other)
extplots <- c('SUPILA','SURIP','SUABCO')
neplots <- length(extplots)

plt.sdl <- numeric(0) 
for(p in 1:nplots){
	plt.sdl <- c(plt.sdl,which(bframe1[,1]==allplots[p]))
}
for(p in 1:neplots){
	plt.sdl <- c(plt.sdl,which(bframe1[,1]==extplots[p]))	
}

##########Only seedlings in selected plots
Seedlings <- bframe1[plt.sdl,]

#######Size classes = <10,10-25,25-50,50-75,75-100,<137 cm
SizeClass <- c(10,25,50,75,100,137)
N.Cls <- length(SizeClass)

Year <- seq(1999,2008,by=1)
NYr <- length(Year)

NSdl <- nrow(Seedlings)
max.col <- ncol(Seedlings)

Species <- c('ABCO','ABMA','ABPS','ABXX','CADE','PICO','PIJE','PILA','PIMO','PIPO','PIXX','PSME','QUCH', 'QUKE','SEGI')
NSpp <- length(Species)

#Make a survival matrix and size matrix. -33 indicates not yet tagged, -22 that seedling is dead
Surv.Mat <- matrix(-33,NSdl,NYr); Size.Mat <- matrix(-33,NSdl,NYr)
#Create species indicators
Spp.Ind <- matrix(0,NSdl,NSpp)
#Create plot indicators
Plt.Ind <- matrix(0,NSdl,nplots+neplots)

for(i in 1:NSdl){
	Spp.Ind[i,which(Species == Seedlings[i,'SPPCODE'])] <- 1
	Plt.Ind[i,which(c(allplots,extplots) == Seedlings[i,'PLOT_NAME'])] <- 1
	for(t in 1:NYr){
		if(Seedlings[i,7+t]>0){ #measured this year?
			Size.Mat[i,t] <- Seedlings[i,7+t]
			Surv.Mat[i,t] <- 1
		}
		if(Seedlings[i,7+t]==-77){ #transition to adult size
			Size.Mat[i,t] <- 150
			Surv.Mat[i,t] <- 1			
		}
		if(Seedlings[i,7+t]==-66){ #known dead
			Size.Mat[i,t] <- 0
			Surv.Mat[i,t] <- 0			
		}
		if(Seedlings[i,7+t]==-88){ #tag recovered, likely dead
			Size.Mat[i,t] <- 0
			Surv.Mat[i,t] <- 0			
		}	
		if(Seedlings[i,7+t]==-99){ #missing
			if(length(which(Seedlings[i,(7+t):max.col]>0))>=1){  #recorded in future
			if(length(which(Seedlings[i,(7+1):(7+t)]>0))>=1){  #recorded in past
				Surv.Mat[i,t] <- 1
				   aa <- which(Seedlings[i,8:(7+t)]>0)
				if(length(aa)>0) Size.Mat[i,t] <- Seedlings[i,7+max(aa)] #last recorded height
			}}
			else {												#not recorded in future
					Size.Mat[i,t] <- 0
					Surv.Mat[i,t] <- 0		
			}		
		}#end "missing" loop
		
	if(Seedlings[i,t+7]==-22 & Seedlings[i,(t+7)-1]>0) { #marked previously dead but clearly not
		Surv.Mat[i,t] <- 0
		Size.Mat[i,t] <- 0 
	}	
	if(Seedlings[i,t+7]==-22 & Seedlings[i,(t+7)-1]==-22) { #marked previously dead this year and last
		Surv.Mat[i,t] <- -22
		Size.Mat[i,t] <- -22 
	}
	if(Seedlings[i,t+7]==-22 & Seedlings[i,(t+7)-1]< -60) { #marked previously dead this year and dead or likely dead earlier
		Surv.Mat[i,t] <- -22
		Size.Mat[i,t] <- -22 
	}
			
	if(t>1){ 
		if(Surv.Mat[i,t]==0 & Surv.Mat[i,t-1]==0) {
		Surv.Mat[i,t] <- -22
		Size.Mat[i,t] <- -22 #previously dead
	}
		if(Surv.Mat[i,t]==0 & Surv.Mat[i,t-1]==-22) {
		Surv.Mat[i,t] <- -22
		Size.Mat[i,t] <- -22 #previously dead
	}}
 }#end year loop
 print(i)
}#end seedling loop


####Basal area and trees near seedling
bframe2 <-data.frame(read.table("PlotInfo2.txt",header=T))
bframe3 <- data.frame(read.table("treeyears.txt",header=T))
bframe4 <-read.csv("tree2.csv",header=T)
bframe5 <-read.csv("quadrat_precise.csv",header=T)

# to calculate distances
distmat <- function(x1,y1,x2,y2){
    xd <- outer(x1,x2,function(x1,x2) (x1 - x2)^2)
    yd <- outer(y1,y2,function(y1,y2) (y1 - y2)^2)
    d <- t(sqrt(xd + yd)) 
    return(d)
}

BA.Mat <- numeric(0); TNum.Mat <- numeric(0)

for (p in 1:nplots){
	Quad <-bframe5[which(bframe5[,1]==allplots[p]),]  #set of quadrat locations for plot
	nquad <- nrow(Quad)
		
	trees <- bframe4[bframe4$PLOT==allplots[p],]    #Trees in plot 
	ntree <- nrow(trees)
	treex <- trees$newX
	treey <- trees$newY
	
    #distances between trees and quadrats
	TQdist <- distmat(treex,treey,Quad[,6],Quad[,7]) 
	
	####calculate basal areas
	treeyears<-bframe3[which(bframe3[,'plot']==allplots[p]),2:8]
	#which census periods match sdl census?
  	rel.census <- which(treeyears>=1999 & treeyears<=2008) 
  	tree.sdl.yrs <- treeyears[rel.census]
  
	ba <- matrix(0,ntree,NYr)
	for(i in 1:ntree){
	for(t in 1:NYr){
		if(Year[t]<tree.sdl.yrs[1]) dbh <- trees[i,(7+rel.census[1])]
	  if(length(tree.sdl.yrs)>1){
		if(Year[t]>=tree.sdl.yrs[1] & Year[t]<tree.sdl.yrs[2]) dbh <- trees[i,
		(7+rel.census[1])]
	  if(length(tree.sdl.yrs)>2){	
		if(Year[t]>=tree.sdl.yrs[2] & Year[t]<tree.sdl.yrs[3]) dbh <- trees[i,
		(7+rel.census[2])]
		if(Year[t]>=tree.sdl.yrs[3]) dbh <- trees[i,(7+rel.census[3])]
	  }
	  else {if(Year[t]>=tree.sdl.yrs[2]) dbh <- trees[i,(7+rel.census[2])]}
	  }	
		
		if(length(which(dbh>=0))>0) ba[i,t]<-pi*((dbh/200)^2) #in meters, not cm
        if(length(which(trees[i,"MortalityYear"]>=0))>0 & trees[i,"MortalityYear"]<=Year[t])	
        ba[i,t]<-0
        if(length(which(trees[i,"IngrowthYear"]>=0))>0 & trees[i,"IngrowthYear"]>Year[t]) 
        ba[i,t]<-0
		}
	} 

##calculate tree number and basal area within 10 m of quadrat
	ba.quad <- matrix(0,nquad,NYr) 
	tn.quad <- matrix(0,nquad,NYr)  
	
	for (j in 1:nquad){
 		q <- which(TQdist[j,]<10)	
 		if(length(q)>0) {
 		tree.sub <- trees[q,]
 		ba.sub <- ba[q,]
 		for (t in 1:NYr){
 			if(length(q)>1){
 			  qq <- which(ba.sub[,t]>0)
 			  tn.quad[j,t] <- length(qq)
 			  ba.quad[j,t] <- sum(ba.sub[,t])
 			}
 			if(length(q)==1){
 			  if(ba.sub[t]>0) tn.quad[j,t] <- 1
 			  ba.quad[j,t] <- ba.sub[t] 
 			} 
 		} #end t loop
		}} #end if dist and nquad loops
	temp.BA <- cbind(rep(p,nquad),Quad[,2:3],Quad[,6:7],ba.quad)
	temp.Num <- cbind(rep(p,nquad),Quad[,2:3],Quad[,6:7],tn.quad)
	
	if(p==1) {
	BA.Mat <- temp.BA
	TNum.Mat <- temp.Num
	}
	if(p>1){
	BA.Mat <- rbind(BA.Mat,temp.BA)
	TNum.Mat <- rbind(TNum.Mat,temp.Num)
	}
	print(p)
} #end plot loop
  ##Add on additional plots

####SUPILA first
	Quad1 <-bframe5[which(bframe5[,1]==extplots[1]),]  #set of quadrat locations for plot
	nquad1 <- nrow(Quad1)  
 
 trees1 <- bframe4[bframe4$PLOT==extplots[1],]    #Trees in plot 
	ntree1 <- nrow(trees1)
	treex1 <- trees1$newX
	treey1 <- trees1$newY
	Qx1 <- Quad1[,6]
	Qy1 <- Quad1[,7]

######Add on SURIP
	Quad2 <-bframe5[which(bframe5[,1]==extplots[2]),]  #set of quadrat locations for plot
	nquad2 <- nrow(Quad2)
		
	trees2 <- bframe4[bframe4$PLOT==extplots[2],]
	trees2 <- trees2[c(1:290,292:length(treey2)),]   #Trees in plot 
	ntree2 <- nrow(trees2)
	treex2 <- trees2$newX
	treey2 <- trees2$newY-125
	Qx2 <- Quad2[,6]
	Qy2 <- Quad2[,7]-125
	
######Add on SUABCO 	
	Quad3 <-bframe5[which(bframe5[,1]==extplots[3]),]  #set of quadrat locations for plot
	nquad3 <- nrow(Quad3)
	
	trees3 <- bframe4[bframe4$PLOT==extplots[3],]
	  #Trees in plot 
	ntree3 <- nrow(trees3)
	treex3 <- trees3$newX+50
	treey3 <- trees3$newY-200
	Qx3 <- Quad3[,6]+50
	Qy3 <- Quad3[,7]-200

	
	treex <- c(treex1,treex2,treex3)
	treey <- c(treey1,treey2,treey3)
	Qx <- c(Qx1,Qx2,Qx3)
	Qy <- c(Qy1,Qy2,Qy3)
	
	ntree <- length(treex)
	trees <- rbind(trees1,trees2,trees3)
	nquad <- length(Qx)
 #distances between trees and quadrats
	TQdist <- distmat(treex,treey,Qx,Qy) 
	
	####calculate basal areas
	treeyears<-rbind(bframe3[which(bframe3[,'plot']==extplots[1]),2:8],bframe3[which(bframe3[,'plot']==extplots[2]),2:8],bframe3[which(bframe3[,'plot']==extplots[3]),2:8])
	#which census periods match sdl census?
  	rel.census <- which(treeyears[1,]>=1999 & treeyears[1,]<=2008) 
  	tree.sdl.yrs <- treeyears[1,rel.census]
  
	ba <- matrix(0,ntree,NYr)
	for(i in 1:ntree){
	for(t in 1:NYr){
		if(Year[t]<tree.sdl.yrs[1]) dbh <- trees[i,(7+rel.census[1])]
	  if(length(tree.sdl.yrs)>1){
		if(Year[t]>=tree.sdl.yrs[1] & Year[t]<tree.sdl.yrs[2]) dbh <- trees[i,
		(7+rel.census[1])]
	  if(length(tree.sdl.yrs)>2){	
		if(Year[t]>=tree.sdl.yrs[2] & Year[t]<tree.sdl.yrs[3]) dbh <- trees[i,
		(7+rel.census[2])]
		if(Year[t]>=tree.sdl.yrs[3]) dbh <- trees[i,(7+rel.census[3])]
	  }
	  else {if(Year[t]>=tree.sdl.yrs[2]) dbh <- trees[i,(7+rel.census[2])]}
	  }	
		if(length(which(dbh>=0))>0) ba[i,t]<-pi*((dbh/200)^2) #in meters, not cm
        if(length(which(trees[i,"MortalityYear"]>=0))>0 & trees[i,"MortalityYear"]<=Year[t])	
        ba[i,t]<-0
        if(length(which(trees[i,"IngrowthYear"]>=0))>0 & trees[i,"IngrowthYear"]>Year[t]) 
        ba[i,t]<-0
		}
	} 

##calculate tree number and basal area within 10 m of quadrat
	ba.quad <- matrix(0,nquad,NYr) 
	tn.quad <- matrix(0,nquad,NYr)  
	
	for (j in 1:nquad){
 		q <- which(TQdist[j,]<10)	
 		if(length(q)>0) {
 		tree.sub <- trees[q,]
 		ba.sub <- ba[q,]
 		for (t in 1:NYr){
 			if(length(q)>1){
 			  qq <- which(ba.sub[,t]>0)
 			  tn.quad[j,t] <- length(qq)
 			  ba.quad[j,t] <- sum(ba.sub[,t])
 			}
 			if(length(q)==1){
 			  if(ba.sub[t]>0) tn.quad[j,t] <- 1
 			  ba.quad[j,t] <- ba.sub[t] 
 			} 
 		} #end t loop
		}} #end if dist and nquad loops

plot.ind <- 	c(rep(23,nquad1),rep(24,nquad2),rep(25,nquad3))
Quad.ind <- rbind(Quad1[,2:3],Quad2[,2:3],Quad3[,2:3])	
		
BA.Temp <- cbind(plot.ind,Quad.ind,Qx,Qy,ba.quad)
TNum.Temp <- cbind(plot.ind,Quad.ind,Qx,Qy,tn.quad)		

colnames(BA.Mat) <- colnames(BA.Temp)
BA.Mat <- rbind(BA.Mat,BA.Temp)
colnames(TNum.Mat) <- colnames(TNum.Temp)
TNum.Mat <- rbind(TNum.Mat,TNum.Temp)

#### Climate
####Need to get new "basins" set up.
bframe6 <- data.frame(read.table("climate_data_long.txt",header=T))
basin <- c(27,1,3,10,9,8,7,17,25,16,24,14,20,26,19,18,0,11,13,15,6,4,21,21,21) 

for(p in 1:(nplots+neplots)){
	q <- which(bframe6[,3]==basin[p])
	if(length(q)>0){
	clim.temp <- bframe6[q,]
	Plt <- rep(p,length(Year))
	precip <- numeric(0); JulMax <- numeric(0)
	JanMin <- numeric(0); AvTemp <- numeric(0)
	snow <- numeric(0); CWD <- numeric(0)
	for(t in Year){
		precip <- c(precip,sum(clim.temp[which(clim.temp[,1]==t),4]))
		snow <- c(snow,sum(clim.temp[which(clim.temp[,1]==t),9]))
		CWD <- c(CWD,sum(clim.temp[which(clim.temp[,1]==t),15]))
		JulMax <- c(JulMax,clim.temp[which(clim.temp[,1]==t & clim.temp[,2]==10),6])
		JanMin <- c(JanMin,clim.temp[which(clim.temp[,1]==t & clim.temp[,2]==4),7])
		AvTemp <- c(AvTemp,mean(clim.temp[which(clim.temp[,1]==t),8]))
	}
	if(p==1) Climate <- cbind(Plt,Year,precip,JulMax,JanMin,AvTemp,snow,CWD)
	if(p>1) Climate <- rbind(Climate,cbind(Plt,Year,precip,JulMax,JanMin,AvTemp,snow,CWD))
	}
}

##Average climate by plot
est.yrs <- seq(1976,2013,by=1)
for(p in 1:(nplots+neplots)){
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


##### Setting up x's and y's
Y.surv <- numeric(0)  #survival vector
Y.grow <- numeric(0)  #growth vector
SPP <- numeric(0)  #Species indicator
Elev <- numeric(0) #elevation
PLT <- numeric(0)  #Plot indicator
YR <- numeric(0)  #Year indicator
SZ <- numeric(0)  #Size-class indicator
T.BA <- numeric(0) #tree BA
T.N <- numeric(0) #tree num
P.C <- numeric(0) #precip current yr
JMx.C <- numeric(0) #JulMax current yr
JMn.C <- numeric(0) #JanMin current yr
AT.C <- numeric(0)#AvTemp current yr
S.C <- numeric(0)#snow current yr
CWD.C <-numeric(0) #CWD current yr
P.P <- numeric(0)#precip prev yr
JMx.P <-numeric(0)#JulMax prev yr
JMn.P <-numeric(0)#JanMin prev yr
AT.P <-numeric(0)#AvTemp prev yr
S.P <-numeric(0)#snow prev yr
CWD.P <-numeric(0)#CWD prev yr

P.CD <- numeric(0) #precip current yr rel to mean
JMx.CD <- numeric(0) #JulMax current yr rel to mean
JMn.CD <- numeric(0) #JanMin current yr rel to mean
AT.CD <- numeric(0)#AvTemp current yr rel to mean
S.CD <- numeric(0)#snow current yr rel to mean
CWD.CD <-numeric(0) #CWD current yr rel to mean
P.PD <- numeric(0)#precip prev yr rel to mean
JMx.PD <-numeric(0)#JulMax prev yr rel to mean
JMn.PD <-numeric(0)#JanMin prev yr rel to mean
AT.PD <-numeric(0)#AvTemp prev yr rel to mean
S.PD <-numeric(0)#snow prev yr rel to mean
CWD.PD <-numeric(0)#CWD prev yr rel to mean

for(i in 1:4000){ #1:NSdl 
	lyr <- which(Surv.Mat[i,]>=0)
	for(t in lyr){
		if(t<10){
		if(Surv.Mat[i,t]==1){
			#The survival vector
			if(Surv.Mat[i,t+1]==1) Y.surv <- c(Y.surv,1)
			if(Surv.Mat[i,t+1]==0) Y.surv <- c(Y.surv,0)
			#The growth vector
			if(Size.Mat[i,t]>=Size.Mat[i,t+1]) Y.grow <- c(Y.grow,0)
			if(Size.Mat[i,t]<Size.Mat[i,t+1])  Y.grow <- c(Y.grow,1)
			
			#The X's
				SPP <- rbind(SPP,Spp.Ind[i,])				#species
				PLT <- rbind(PLT,Plt.Ind[i,])				#plot
					yr <- rep(0,NYr); yr[t+1]<- 1
				Elev <- c(Elev,bframe2[which(bframe2[,     1]==allplots[which(Plt.Ind[i,]==1)]),2])                   #elevation
				YR <- rbind(YR,yr) 							#year
					a <- which(SizeClass==Size.Mat[i,t])
					sz<- rep(0,N.Cls); sz[a] <- 1
				SZ<- rbind(SZ,sz)   						#size class
					qd <- which(BA.Mat[,1]==which(Plt.Ind[i,]==1) & BA.Mat[,
				       2]==Seedlings[i,2] & BA.Mat[,3]==Seedlings[i,3])
				T.BA <- c(T.BA,BA.Mat[qd,5+t])				#local basal area
				T.N <- c(T.N,TNum.Mat[qd,5+t])				#local tree number
				
				    clm.c <- which(Climate[,1]==which(Plt.Ind[i,]==1) & 
				    Climate[,2]==Year[t+1])
				    clm.p <- which(Climate[,1]==which(Plt.Ind[i,]==1) & 
				    Climate[,2]==Year[t])	
				    				
			  P.C <- c(P.C,Climate[clm.c,3]); P.P <- c(P.P,Climate[clm.p,3])
			JMx.C <- c(JMx.C,Climate[clm.c,4]); JMx.P <- c(JMx.P,Climate[clm.p,4])
		JMn.C <- c(JMn.C,Climate[clm.c,5]); JMn.P <- c(JMn.P,Climate[clm.p,5])
		AT.C <- c(AT.C,Climate[clm.c,6]); AT.P <- c(AT.P,Climate[clm.p,6])
			  S.C <- c(S.C,Climate[clm.c,7]); S.P <- c(S.P,Climate[clm.p,7])
			  CWD.C <- c(CWD.C,Climate[clm.c,8]); CWD.P <- c(CWD.P,Climate[clm.p,8])
			  
			 P.CD <- c(P.CD,Climate[clm.c,3]-Climate.mean[which(Plt.Ind[i,]==1),2])
			 P.PD <-c(P.PD,Climate[clm.p,3]-Climate.mean[which(Plt.Ind[i,]==1),2])
		JMx.CD <- c(JMx.CD,Climate[clm.c,4]-Climate.mean[which(Plt.Ind[i,]==1),3]) 
		JMx.PD <- c(JMx.PD,Climate[clm.p,4]-Climate.mean[which(Plt.Ind[i,]==1),3])
		JMn.CD <- c(JMn.CD,Climate[clm.c,5]-Climate.mean[which(Plt.Ind[i,]==1),4]) 
		JMn.PD <- c(JMn.PD,Climate[clm.p,5]-Climate.mean[which(Plt.Ind[i,]==1),4])
			AT.CD <- c(AT.CD,Climate[clm.c,6]-Climate.mean[which(Plt.Ind[i,]==1),5]) 
			AT.PD <- c(AT.PD,Climate[clm.p,6]-Climate.mean[which(Plt.Ind[i,]==1),5])
			S.CD <- c(S.CD,Climate[clm.c,7]-Climate.mean[which(Plt.Ind[i,]==1),6]) 
			S.PD <- c(S.PD,Climate[clm.p,7]-Climate.mean[which(Plt.Ind[i,]==1),6])
		CWD.CD <- c(CWD.CD,Climate[clm.c,8]-Climate.mean[which(Plt.Ind[i,]==1),7]) 
		CWD.PD <- c(CWD.PD,Climate[clm.p,8]-Climate.mean[which(Plt.Ind[i,]==1),7])
			  
			  }# end if sdl. alive  loop
		} #end "if t < 10"
	} #end t loop
print(i)
if(length(Y.surv)!=length(T.BA)) break
} #end i loop

#Length Y.surv = 63904

#Fire parameter
N1 <- length(Y.surv) 

Fire <- matrix(0,N1,3)#disturbed this year, -1 yr, fire 2-7 years ago
for(j in 1:N1){
	if(PLT[j,1]==1){ #YOHOPIPO fire in 2007 
		if(which(YR[j,]==1)==9) Fire[j,1]<-1
		if(which(YR[j,]==1)==10) Fire[j,2]<-1
	}
	if(PLT[j,5]==1){ #FFS6Burn fire in 2001 fall (so 2002 "current")
		if(which(YR[j,]==1)==4) Fire[j,1]<-1
		if(which(YR[j,]==1)==5) Fire[j,2]<-1
		if(which(YR[j,]==1)==6) Fire[j,3]<-1
		if(which(YR[j,]==1)==7) Fire[j,3]<-1
		if(which(YR[j,]==1)==8) Fire[j,3]<-1
		if(which(YR[j,]==1)==9) Fire[j,3]<-1
		if(which(YR[j,]==1)==10) Fire[j,3]<-1
	}	
		if(PLT[j,6]==1){ #FFS5Burn fire in 2001 fall (so 2002 "current")
		if(which(YR[j,]==1)==4) Fire[j,1]<-1
		if(which(YR[j,]==1)==5) Fire[j,2]<-1
		if(which(YR[j,]==1)==6) Fire[j,3]<-1
		if(which(YR[j,]==1)==7) Fire[j,3]<-1
		if(which(YR[j,]==1)==8) Fire[j,3]<-1
		if(which(YR[j,]==1)==9) Fire[j,3]<-1
		if(which(YR[j,]==1)==10) Fire[j,3]<-1
	}	
		if(PLT[j,7]==1){ #FFS2Burn fire in 2001 fall (so 2002 "current")
		if(which(YR[j,]==1)==4) Fire[j,1]<-1
		if(which(YR[j,]==1)==5) Fire[j,2]<-1
		if(which(YR[j,]==1)==6) Fire[j,3]<-1
		if(which(YR[j,]==1)==7) Fire[j,3]<-1
		if(which(YR[j,]==1)==8) Fire[j,3]<-1
		if(which(YR[j,]==1)==9) Fire[j,3]<-1
		if(which(YR[j,]==1)==10) Fire[j,3]<-1
	}
		if(PLT[j,8]==1){ #LOTHAR fire in 2004 (so 2005 "current")
		if(which(YR[j,]==1)==7) Fire[j,1]<-1
		if(which(YR[j,]==1)==8) Fire[j,2]<-1
		if(which(YR[j,]==1)==9) Fire[j,3]<-1
		if(which(YR[j,]==1)==10) Fire[j,3]<-1
	}
		if(PLT[j,9]==1){ #UPTHAR fire in 2004 (so 2005 "current")
		if(which(YR[j,]==1)==7) Fire[j,1]<-1
		if(which(YR[j,]==1)==8) Fire[j,2]<-1
		if(which(YR[j,]==1)==9) Fire[j,3]<-1
		if(which(YR[j,]==1)==10) Fire[j,3]<-1
	}
	print(j)
}

save.image(file="TagSdl_Data5.RData")