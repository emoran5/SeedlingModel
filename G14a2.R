#Growth model
#G14 - Size+ Species + Fire X Size + Clim1 X Spp + Clim2 X Spp + BA X SPP 

load("TagSdl_Data5.RData")

library(mvtnorm)

###functions
b.update <- function(LM,X1,sinv,vinvert,bprior){  ##Update betas 
    sx <- crossprod(X1)*sinv    
    sy <- crossprod(X1,LM)*sinv
  bigv <- solve(sx+vinvert)  #inverse of (n/sig + 1/priorvar)
  smallv <- sy+vinvert%*%bprior  #(n*meany/sig + priormean/priorvar)
  b <- t(rmvnorm(1,bigv%*%smallv,bigv))
  return(b)
}

v.update <- function(Lm,X,beta,N){  #update variance
   sx <- crossprod((Lm-X%*%beta))
    u1 <- S1 + 0.5*N +1
    u2 <- S2 + 0.5*sx          
  return(1/rgamma(1,u1,u2))
}

lm.calc <- function(X,beta,sig,N){
	q <- X%*%beta
	e <- rnorm(N,0,sqrt(sig))
	lm <- q+e
	return(lm)
}

l1func <- function (Lm,Y,N){   #Calculates likelihood of Y given current parameters
	theta <- inv.logit(Lm)
	l1 <- rep(1,N1)
	for (i in 1:N1){
		l1[i] <- (theta[i]^Y[i])*((1-theta[i])^(1-Y[i]))
	}
	return(l1)
}

logit <- function(x){
	a <- log(x/(1-x))
	return(a)
}
inv.logit <- function(x){
	a <- exp(x)/(1+exp(x))
}

###########################Simple data generation
###Pull out growth data for seedlings that are alive
Live <- which(Y.surv==1)
Gro <- Y.grow[Live] 
N1 <- length(Gro)

#note: divide spring snow and precip by 1000 (so in m units, not mm)
S.C <- S.C[Live]/1000; P.C <- P.C[Live]/1000
S.P <- S.P[Live]/1000; P.P <- P.P[Live]/1000
#divide CWD by 100 (so in 10-cm units, not mm)
CWD.C <- CWD.C[Live]/100; CWD.P <- CWD.P[Live]/100;
#divide temp by 10 (so in 1/10 degrees):
JMx.C <- JMx.C[Live]/10; JMn.C <- JMn.C[Live]/10; AT.C <- AT.C[Live]/10
JMx.P <- JMx.P[Live]/10; JMn.P <- JMn.P[Live]/10; AT.P <- AT.P[Live]/10

JMx.CD <- JMx.CD[Live];JMn.CD <- JMn.CD[Live]
JMx.PD <- JMx.PD[Live];JMn.PD <- JMn.PD[Live]

AT.CD <- AT.CD[Live];AT.PD <- AT.PD[Live]
S.CD <- S.CD[Live]/100; P.CD <- P.CD[Live]/100
S.PD <- S.PD[Live]/100; P.PD <- P.PD[Live]/100
CWD.CD <- CWD.CD[Live]/100; CWD.PD <- CWD.PD[Live]/100

#Reduced size indicator; 10, 25-50, 75-135
s.10 <- SZ[Live,1]; s.25_50 <- apply(SZ[Live,2:3],1,sum); s.75 <- apply(SZ[Live,4:6],1,sum)
SZ.red <- cbind(s.10,s.25_50,s.75)

#Reduced species: ABCO, ABMA,"maybe abies",CADE,PILA,PIMO,PIPO,"Other pines",PSME,"Oaks", SEGI
AB <- apply(SPP[Live,3:4],1,sum)
PI <- apply(SPP[Live,c(6:7,11)],1,sum)
QU <- apply(SPP[Live,13:14],1,sum)
SPP.red <- cbind(SPP[Live,1:2],AB,SPP[Live,c(5,8,9,10)],PI,SPP[Live,12],QU,SPP[Live,15])

######Select Climate Variables
#Clim1 <- JMx.CD; Clim2 <- S.CD #A1
#Clim1 <- JMx.CD; Clim2 <- P.CD #B1
#Clim1 <- JMx.PD; Clim2 <- P.CD #H1
#Clim1 <- JMx.CD; Clim2 <- P.PD #K1
#Clim1 <- JMx.CD; Clim2 <- (S.CD+S.PD)/2 #N1
#Clim1 <- (JMx.CD+JMx.PD)/2; Clim2 <- (S.CD+S.PD)/2 #O1

Clim1 <- JMx.C; Clim2 <- S.C #A2
#Clim1 <- JMx.C; Clim2 <- P.C #B2
#Clim1 <- JMx.P; Clim2 <- P.C #H2
#Clim1 <- JMx.C; Clim2 <- P.P #K2
#Clim1 <- JMx.C; Clim2 <- (S.C+S.P)/2 #N2
#Clim1 <- (JMx.C+JMx.P)/2; Clim2 <- (S.C+S.P)/2 #O2

T.BA <- T.BA[Live]

Fire <- Fire[Live,]

##Size X Fire indicator
SZ.Fire <- matrix(0,N1,9)
for(i in 1:N1){
	a <- which(SZ.red[i,]==1); b <- which(Fire[i,]==1)
	if(length(a)>0 & length(b)>0){
	if(a==1 & b==1) SZ.Fire[i,1]<- 1
	if(a==1 & b==2) SZ.Fire[i,2]<- 1
	if(a==1 & b==3) SZ.Fire[i,3]<- 1
	if(a==2 & b==1) SZ.Fire[i,4]<- 1
	if(a==2 & b==2) SZ.Fire[i,5]<- 1
	if(a==2 & b==3) SZ.Fire[i,6]<- 1
	if(a==3 & b==1) SZ.Fire[i,7]<- 1
	if(a==3 & b==2) SZ.Fire[i,8]<- 1
	if(a==3 & b==3) SZ.Fire[i,9]<- 1
	}
}

### Climate X Species matrix
SPP.Clim1 <- matrix(0,N1,11)
  for(i in 1:N1){
  	SPP.Clim1[i,which(SPP.red[i,]==1)]<-Clim1[i]
  	}
SPP.Clim2 <- matrix(0,N1,11)
  for(i in 1:N1){
  	SPP.Clim2[i,which(SPP.red[i,]==1)]<-Clim2[i]
  	}  	

### Basal area X Species matrix
SPP.red.BA <- matrix(0,N1,11)
  for(i in 1:N1){
  	SPP.red.BA[i,which(SPP.red[i,]==1)]<-T.BA[i]
  	}

####### Setup design matrix
X1 <- cbind(rep(1,N1),SZ.red,SPP.red,SZ.Fire,SPP.Clim1,SPP.Clim2,SPP.red.BA) 

I <- ncol(X1)  #columns of design matrix

###### Priors
Bpm <- c(-2.6,rep(0,3),rep(0,11),rep(0,9),rep(0,11),rep(0,11),rep(0,11)) #Intercept,size (red), BA
Bps <- c(3,rep(3,3),rep(3,11),rep(3,9),rep(3,11),rep(3,11),rep(3,11))   #variances
bprior <- matrix(Bpm,I,1)  #beta mean matrix
vinvert <- solve(diag(Bps))  # 1/beta variances

 ###IG prior for error parameter
#S1 <- 2.2 #must be >2.  Find appropriate as follows...#2.01
#  IGM <- 0.2  #desired variance mean (square of sd) # - 0.16
#S2 <- (S1-1)*IGM 
# IGV <- (S2*S2)/((S1-1)^2)/(S1-2); IGV  #check if has high enough variance - sd 1.6

S1 <-2.1
S2 <-0.35


###Starting values
Beta <- c(runif(1,-3,0),runif(3,-1,1),runif(11,-1,1),runif(9,-1,1),runif(11,-0.5,0.5),runif(11,-0.5,0.5),runif(11,-0.5,0.5))
Sig <- runif(1,0.05,0.5)  

LM <- lm.calc(X1,Beta,Sig,N1)

ncyc <- 25000   #Number of Gibbs Steps
bgibbs <- matrix(NA,I,ncyc+1)
sgibbs <- rep(NA,ncyc+1)


bgibbs[,1] <- Beta
sgibbs[1] <- Sig

thin <- seq(15000,(ncyc+1),by=20) #Burn-in and thinning

for (g in 1:ncyc){
	#update beta and sigma
	sinv <- 1/Sig
	Beta <- b.update(LM,X1,sinv,vinvert,bprior)
	Sig <- v.update(LM,X1,Beta,N1)
	bgibbs[,g+1] <- Beta
	sgibbs[g+1] <- Sig
	
	Lm.new <- lm.calc(X1,Beta,Sig,N1)
	
	l1.now <- l1func(LM,Gro,N1)
	l1.new <- l1func(Lm.new,Gro,N1)	
	
	pnow <- log(l1.now); pnew <- log(l1.new)
	r <- exp(pnew-pnow)
	
	z <- runif(N1,0,1)
	
	for(i in 1:N1){
		if (r[i]>z[i]) LM[i] <- Lm.new[i]
	}
	
}

beta.thin <- bgibbs[,thin]
sig.thin <- sgibbs[thin]

#### Posterior parameter output
outpar <- matrix(NA,I+1,4)
outpar[,1]<- c(apply(beta.thin,1,mean),mean(sig.thin))
outpar[,2]<- c(apply(beta.thin,1,sd),sd(sig.thin))
outpar[,3]<- c(apply(beta.thin,1, quantile, probs= 0.025),quantile(sig.thin,0.025))
outpar[,4]<- c(apply(beta.thin,1, quantile, probs= 0.975),quantile(sig.thin,0.975))
colnames(outpar) <- c('estimate','se','.025','.975')
outpar    

######### Predictive Loss Calculation
ng <- length(thin)
ypred <- matrix(0,ng,N1)
for (j in 1:ng){
	bg <- beta.thin[,j]; sd <- sqrt(sig.thin[j])
	
	q <- X1%*%bg
	e <- rnorm(N1,0,sd)
	lam <- q+e
	theta <- inv.logit(lam)
	for(i in 1:N1){
	ypred[j,i] <- rbinom(1,1,theta[i])		
	}	
}

ym <- apply(ypred,2,mean)
yv <- apply(ypred,2,var)

gm <- sum((ym-Gro)^2)
pm <- sum(yv)
Dm <- gm +pm
Dm  
