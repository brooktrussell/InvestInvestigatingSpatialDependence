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
nobs_vec <- c(1250,2500,5000,10000)#seq(500,7500,by=500)
library(evd)

set.seed(1)
gammaBar_hat_distMat <- matrix(NA,nsim,length(nobs_vec))
chi_hat_distMat <- matrix(NA,nsim,length(nobs_vec))
#for alpha=.35
for (i in 1:nsim){
  for (j in 1:length(nobs_vec)){
    trashdat <- rbvevd(nobs_vec[j],dep=.35,model="log",mar1 = c(1,1,1),mar2=c(1,1,1))
    gammaBar_hat_distMat[i,j] <- ChiCalc(trashdat[,1],trashdat[,2],qtile=.99)$chi.hat
    chi_hat_distMat[i,j] <- 1 - GammaCalc(trashdat[,1],trashdat[,2],qtile=.99)$gamma.hat
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
    gammaBar_hat_distMat2[i,j] <- ChiCalc(trashdat[,1],trashdat[,2],qtile=.99)$chi.hat
    chi_hat_distMat2[i,j] <- 1 - GammaCalc(trashdat[,1],trashdat[,2],qtile=.99)$gamma.hat
  }
}

pdf(file="SimStudyBoxPlots_LowDep_12042025.pdf",w=8.5,h=9)
par(mfrow=c(2,1))
boxplot(c(gammaBar_hat_distMat2)~rep(1:length(nobs_vec),each=nsim),ylim=c(0,.6),xaxt="n",xlab="Sample Size",ylab=expression(hat(bar(gamma))))
axis(1,at=1:length(nobs_vec),labels=nobs_vec)
boxplot(c(chi_hat_distMat2)~rep(1:length(nobs_vec),each=nsim),ylim=c(0,.6),xaxt="n",xlab="Sample Size",ylab=expression(hat(chi)))
axis(1,at=1:length(nobs_vec),labels=nobs_vec)
dev.off()

angdist_log <- function(w,a){.5*(1/a - 1)*(w*(1-w))^(-1 - 1/a)*(w^(-1/a) + (1-w)^(-1/a))^(a-2)}
angdist_log <- Vectorize(angdist_log)

pdf(file="SimStudyBoxPlots.pdf",w=8.5,h=9)
par(mfrow=c(2,2))
boxplot(c(gammaBar_hat_distMat)~rep(1:length(nobs_vec),each=nsim),ylim=c(0,1),xaxt="n",xlab="Sample Size",ylab=expression(hat(bar(gamma))))
axis(1,at=1:length(nobs_vec),labels=nobs_vec)
integrand <- function(x) {abs(2*x-1)*angdist_log(x,a=0.35)}
abline(h=1 - integrate(integrand,lower=0,upper=1)$value)
boxplot(c(chi_hat_distMat)~rep(1:length(nobs_vec),each=nsim),ylim=c(0,1),xaxt="n",xlab="Sample Size",ylab=expression(hat(chi)))
axis(1,at=1:length(nobs_vec),labels=nobs_vec)
abline(h=2-2^.35)
#
boxplot(c(gammaBar_hat_distMat2)~rep(1:length(nobs_vec),each=nsim),ylim=c(0,1),xaxt="n",xlab="Sample Size",ylab=expression(hat(bar(gamma))))
axis(1,at=1:length(nobs_vec),labels=nobs_vec)
integrand <- function(x) {abs(2*x-1)*angdist_log(x,a=0.85)}
abline(h=1 - integrate(integrand,lower=0,upper=1)$value)
boxplot(c(chi_hat_distMat2)~rep(1:length(nobs_vec),each=nsim),ylim=c(0,1),xaxt="n",xlab="Sample Size",ylab=expression(hat(chi)))
axis(1,at=1:length(nobs_vec),labels=nobs_vec)
abline(h=2-2^.85)
dev.off()



#################################
#Simulate two data sets, one with high AD and one with low AD

library(evd)

unitFr<-function(dat.vec){n<-length(na.omit(dat.vec))
trash.rank<-rank(dat.vec)
U<-ifelse(is.na(dat.vec)==TRUE,NA,trash.rank/(n+1))
F<- -1/log(U);return(F)}

angdist_log <- function(w,a){.5*(1/a - 1)*(w*(1-w))^(-1 - 1/a)*(w^(-1/a) + (1-w)^(-1/a))^(a-2)}
angdist_log <- Vectorize(angdist_log)

wseq <- seq(.01,.99,length=100)

sim1_gum <- rbvevd(1000,dep=.35,model="log")
sim1_frech <- apply(sim1_gum,2,FUN=unitFr)
hist_dat <- sim1_frech[which(rowSums(sim1_frech) > quantile(rowSums(sim1_frech),.95)),]
pdf(file="~/Downloads/Sim1_plots.pdf",w=8.5,h=3.3)
par(mfrow=c(1,3),mar = c(5.1, 4.1, 4.1, 2.1))
plot(sim1_gum,xlab=expression(X[1]),ylab=expression(X[2]))
plot(sim1_frech,xlab=expression(Y[1]),ylab=expression(Y[2]))
hist(hist_dat[,1]/rowSums(hist_dat),main="",ylim=c(0,2.35),freq = FALSE,xlab="Angular Component",breaks=seq(0,1,by=.2))
lines(wseq,angdist_log(wseq,a=.35))
dev.off()
1-c(unlist(GammaCalc(sim1_frech,.95)))
ChiCalc(sim1_frech,.95)
2 - 2^.35


sim1_gum <- rbvevd(1000,dep=.85,model="log")
sim1_frech <- apply(sim1_gum,2,FUN=unitFr)
hist_dat <- sim1_frech[which(rowSums(sim1_frech) > quantile(rowSums(sim1_frech),.95)),]
pdf(file="~/Downloads/Sim2_plots.pdf",w=8.5,h=3.3)
par(mfrow=c(1,3),mar = c(5.1, 4.1, 4.1, 2.1))
plot(sim1_gum,xlab=expression(X[1]),ylab=expression(X[2]))
plot(sim1_frech,xlab=expression(Y[1]),ylab=expression(Y[2]))
hist(hist_dat[,1]/rowSums(hist_dat),main="",ylim=c(0,2.25),freq = FALSE,xlab="Angular Component",breaks=seq(0,1,by=.2))
lines(wseq,angdist_log(wseq,a=.85))
dev.off()
1-c(unlist(GammaCalc(sim1_frech,.95)))
ChiCalc(sim1_frech,.95)
2 - 2^.85


nsims <- 50
useq <- seq(.9,.999,by=.001)

gamma_norm_mat <- matrix(NA,nsims,length(useq))
chi_norm_mat <- matrix(NA,nsims,length(useq))

library(mvtnorm)

set.seed(1)
for (i in 1:nsims){
  for (j in 1:length(useq)){
    bvn_dat <- rmvnorm(50000,sigma = matrix(c(1,.6,.6,1),2,2))
    gamma_norm_mat[i,j] <- GammaCalc(bvn_dat[,1],bvn_dat[,2],qtile=useq[j])$gamma.hat
    chi_norm_mat[i,j] <- ChiCalc(bvn_dat[,1],bvn_dat[,2],qtile=useq[j])$chi.hat
  }
}


pdf(file="~/Downloads/NormalAsyDep_122025.pdf",8.5,4.5)
par(mfrow=c(1,2))
matplot(y=1-t(gamma_norm_mat),x=useq,type="l",lty=1,col="gray",xlab="Quantile",ylab=expression(hat(bar(gamma))))
matplot(y=t(chi_norm_mat),x=useq,type="l",lty=1,col="gray",xlab="Quantile",ylab=expression(hat(chi)))
dev.off()
