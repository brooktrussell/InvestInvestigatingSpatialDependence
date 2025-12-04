#Code to replicate the simulation study in the
#paper "Investigating Spatial Dependence in the Degree 
#of Asymptotic Dependence between a Satellite 
#Precipitation Product and Station Data in the Northern 
#US Rocky Mountains"

#Define asymptotic dependence estimators:
#
#function to calculate gamma hat
GammaCalc<-function(x,y,qtile=.97){
  #remove observations with missing values
  new.dat<-na.omit(cbind(x,y))
  
  #uses 'unitFr' function
  unitFr<-function(dat.vec){n<-length(na.omit(dat.vec))
  trash.rank<-rank(dat.vec)
  U<-ifelse(is.na(dat.vec)==TRUE,NA,trash.rank/(n+1))
  F<- -1/log(U);return(F)}
  #uses 'gam.hat' function
  gam.hat<-function(x.vec,y.vec){sum(abs(x.vec-y.vec)/(x.vec+y.vec))/length(x.vec)}
  
  #transform x and y vectors to unit Frechet
  #using a rank transform
  x.fr<-unitFr(new.dat[,1])
  y.fr<-unitFr(new.dat[,2])
  rad.fr<-x.fr+y.fr
  mat<-cbind(x.fr,y.fr)
  
  cutoff <- quantile(rad.fr,qtile)
  mat2<-cbind(mat,rad.fr)
  mat3<-subset(mat2,mat2[,3]>=cutoff)
  gamma.hat<-gam.hat(mat3[,1],mat3[,2])
  
  return(list(gamma.hat=gamma.hat))#return gamma hat

}#end GammaCalc function


#function to calculate chi hat
ChiCalc<-function(x,y,qtile=.97){#P(y>u|x>u)
  #remove observations with missing values
  new.dat<-na.omit(cbind(x,y))
  
  #uses 'unitFr' function
  unitFr<-function(dat.vec){n<-length(na.omit(dat.vec))
  trash.rank<-rank(dat.vec)
  U<-ifelse(is.na(dat.vec)==TRUE,NA,trash.rank/(n+1))
  F<- -1/log(U);return(F)}
  
  x.u <- quantile(new.dat[,1],qtile)
  y.u <- quantile(new.dat[,2],qtile)
  
  mat2<-subset(new.dat,new.dat[,1]>x.u)
  mat3<-subset(mat2,mat2[,2]>y.u)
  chi.hat<-nrow(mat3)/nrow(mat2)
  
  return(list(chi.hat=chi.hat))#return chi hat
}#end chicalc function


#begin simulations
#for alpha=.35,.85
#for alpha=.35,.85
nsim <- 10000
nobs_vec <- seq(500,7500,by=500)
library(evd)

set.seed(1)
gammaBar_hat_distMat <- matrix(NA,nsim,length(nobs_vec))
chi_hat_distMat <- matrix(NA,nsim,length(nobs_vec))
#for alpha=.35
for (i in 1:nsim){
  for (j in 1:length(nobs_vec)){
    trashdat <- rbvevd(nobs_vec[j],dep=.35,model="log",mar1 = c(1,1,1),mar2=c(1,1,1))
    gammaBar_hat_distMat[i,j] <- ChiCalc(trashdat[,1],trashdat[,2],qtile=.95)$chi.hat
    chi_hat_distMat[i,j] <- 1 - GammaCalc(trashdat[,1],trashdat[,2],qtile=.95)$gamma.hat
  }
}

pdf(file="SimStudyBoxPlots_HighDep_12042025.pdf",w=8.5,h=9)
par(mfrow=c(2,1))
boxplot(c(gammaBar_hat_distMat)~rep(1:length(nobs_vec),each=nsim),ylim=c(.44,.96),xaxt="n",xlab="Sample Size",ylab=expression(hat(bar(gamma))))
axis(1,at=1:length(nobs_vec),labels=nobs_vec)
boxplot(c(chi_hat_distMat)~rep(1:length(nobs_vec),each=nsim),ylim=c(.44,.96),xaxt="n",xlab="Sample Size",ylab=expression(hat(chi)))
axis(1,at=1:length(nobs_vec),labels=nobs_vec)
dev.off()

gammaBar_hat_distMat2 <- matrix(NA,nsim,length(nobs_vec))
chi_hat_distMat2 <- matrix(NA,nsim,length(nobs_vec))
#for alpha=.85
for (i in 1:nsim){
  for (j in 1:length(nobs_vec)){
    trashdat <- rbvevd(nobs_vec[j],dep=.85,model="log",mar1 = c(1,1,1),mar2=c(1,1,1))
    gammaBar_hat_distMat2[i,j] <- ChiCalc(trashdat[,1],trashdat[,2],qtile=.95)$chi.hat
    chi_hat_distMat2[i,j] <- 1 - GammaCalc(trashdat[,1],trashdat[,2],qtile=.95)$gamma.hat
  }
}

pdf(file="SimStudyBoxPlots_LowDep_12042025.pdf",w=8.5,h=9)
par(mfrow=c(2,1))
boxplot(c(gammaBar_hat_distMat2)~rep(1:length(nobs_vec),each=nsim),ylim=c(0,.6),xaxt="n",xlab="Sample Size",ylab=expression(hat(bar(gamma))))
axis(1,at=1:length(nobs_vec),labels=nobs_vec)
boxplot(c(chi_hat_distMat2)~rep(1:length(nobs_vec),each=nsim),ylim=c(0,.6),xaxt="n",xlab="Sample Size",ylab=expression(hat(chi)))
axis(1,at=1:length(nobs_vec),labels=nobs_vec)
dev.off()
