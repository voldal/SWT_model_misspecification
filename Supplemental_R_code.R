#This document provides R code to implement the methods described in the paper 'Model misspecification in stepped wedge trials: Random effects for time or treatment'
#Emily C. Voldal, Fan Xia, Avi Kenny, Patrick J. Heagerty, and James P. Hughes

##############################################################################
#
#Outline
#
##############################################################################

#1. Introduction/notes
#2. Functions required for both cases
#3. Time-fitted random treatment case (Case 1)
#4. Treatment-fitted random time case (Case 2)
#5. Use and examples

##############################################################################
#
#1. Introduction/notes
#
##############################################################################


#This code was written in R, version 3.6.1

#This code calculates validity and efficiency (or the associated variances) based on a specific design and true variance components
#See main paper for examples of application

#If you're not interested in details, just run sections 1-4 and scroll down to the 'use and examples' section to get started.

#As written, this code is for Root 1 and assumes the fixed effect for time is linear
#These can be easily modified; see comments for guidance
#Notation is generally consistent with the paper, unless noted otherwise in comments
#For readability, we'll use false, true, and optimal variances as shorthand for var_{m}(\hat\theta), var_{s}(\hat\theta), and var_{m}(\hat\theta_t)

#For updates and any corrections, see https://github.com/voldal/SWT_model_misspecification

#Packages used:
library(swCRTdesign)#Version 3.1
#(See package documentation and paper https://doi.org/10.1016/j.cmpb.2020.105514 for details on use)
library(nleqslv)#Version 3.3.2
#Only needed for Case 2
#Can be easily replaced with your favorite numerical root-finding package/method

#This code works for any classic design, e.g.:
swDsn(c(1,1))
swDsn(c(1,1,1,1,1))
#It also works for SWT designs with extra time before, after, or between cross-overs, e.g.:
swDsn(c(1,1),all.ctl.time0 = FALSE)
swDsn(c(1,1),extra.time=2)
swDsn(c(0,0,1,1))
swDsn(c(0,1,0,1,0))
#It is not currently adapted to work if there are wash-out periods, different number of subjects between time points or clusters, etc.
#Any design built using swDsn where the 'clusters' argument is a vector of 0's and 1's (and tx.effect.frac is left at its default) should be acceptable

##############################################################################
#
#2. Functions required for both cases
#
##############################################################################

#These 'secondary' functions are called by functions for calculations in both cases
#These functions return the covariance matrix for the fixed effects in each model
#But they are not scaled by N (that is, this would be the covariance if there was one cluster per sequence)
#Variance component arguments (sigmasq, tausq, etc.) should be the true values if the model is correctly specified, and the roots misspecified parameters converge to if the model is misspecified

#As written, functions in this file include three fixed effects: an intercept, a time slope, and a treatment effect.
#Change the design matrix Xmat.i (and update the position of the treatment indicator) if you want to model time in a different way (splines, indicator variables, etc.)

#####
#Function: calculate covariance for a model with random treatment effects
#####

#For a model where we specified random treatment (in addition to the random intercept effect), what is the covariance matrix of the fixed effects?

GetCov_thisdesign_trt <- function(design,K,sigmasq,tausq,etasq){
  
  #Info from design
  design.mat <- design$swDsn#Design matrix
  M <- design$n.waves#Total number of waves
  J <- design$total.time#Total number of time periods
  
  RunningSum <- 0#Add up pieces here
  
  for(i in 1:M){#Sum over sequences
    
    #Design matrix X (fixed effects) For K=1:
    Xmat.i <- matrix(c(rep(1,length=J),0:(J-1),design.mat[i,]),nrow=J,ncol=3)
    #Expand X for the correct K:
    Xmat.i <- Xmat.i[rep(1:J,each=K),]#Keep it ordered by time
    
    #Design matrix Z for random treatment effect
    Zmat.trt.i <- as.matrix(cbind(rep(1,length=J),design.mat[i,]))
    Zmat.trt.i <- Zmat.trt.i[rep(1:J,each=K),]
    
    #Matrix G (covariance matrix for random effects)
    G <- diag(c(tausq,etasq))
    
    #Covariance matrix
    SigmaMat.i <- Zmat.trt.i %*% G %*% t(Zmat.trt.i)+sigmasq*diag(J*K)  
    #Element inside the sum
    Piece.i <- t(Xmat.i) %*% solve(SigmaMat.i) %*% Xmat.i#Sum these up over all sequences, then invert and take [3,3]
    RunningSum <- RunningSum+Piece.i
    
  }
  

  return(solve(RunningSum))#This returns A^-1 (with N already taken out)
  #For model-based variance of treatment estimate, take [3,3] of this.  For Sandwich, use this for bread.
}

#####
#Function: calculate covariance for a model with random time effects
#####

#For a model where we specified random time (in addition to the random intercept effect), what is the covariance matrix of the fixed effects?

GetCov_thisdesign_time <- function(design,K,sigmasq,tausq,gammasq){
  
  #Info from design
  design.mat <- design$swDsn#Design matrix
  M <- design$n.waves#Total number of waves
  J <- design$total.time#Total number of time periods
  
  RunningSum <- 0#Add up pieces here
  
  for(i in 1:M){#Sum over sequences
    
    #Design matrix X (fixed effects) For K=1:
    Xmat.i <- matrix(c(rep(1,length=J),0:(J-1),design.mat[i,]),nrow=J,ncol=3)
    #Expand X for the correct K:
    Xmat.i <- Xmat.i[rep(1:J,each=K),]#Keep it ordered by time
    
    #Design matrix Z for random treatment effect
    Zmat.time.i <- as.matrix(cbind(rep(1,length=J),diag(J)))
    Zmat.time.i <- Zmat.time.i[rep(1:J,each=K),]
    
    #Matrix G (covariance matrix for random effects)
    G <- diag(c(tausq,rep(gammasq,length=J)))
    
    #Covariance matrix
    SigmaMat.i <- Zmat.time.i %*% G %*% t(Zmat.time.i)+sigmasq*diag(J*K) 
    #Element inside the sum
    Piece.i <- t(Xmat.i) %*% solve(SigmaMat.i) %*% Xmat.i#Sum these up over all sequences, then invert and take [3,3]
    RunningSum <- RunningSum+Piece.i
    
  }
  
  return(solve(RunningSum))#This returns A^-1 (with N already taken out)
  #For model-based variance of treatment estimate, take [3,3] of this.  For Sandwich, use this for bread.
  
}



##############################################################################
#
#3. Time-fitted random treatment case (Case 1)
#
##############################################################################


#####
#Function: calculate 'meat' for a sandwich estimate, for a misspecified model with time random effects
#####

#This 'secondary' function is called by GetVars_thisdesign_Case1_Root1
#This returns the 'meat' (or B) for the Sandwich-form estimate of the true variance of a misspecified model 
#Requires both true variances and roots of variances from the misspecified model

#Consistent with the 'bread' functions, there is no N here.

GetMeat_thisdesign_falsetime <- function(design,K,sigmasq,tausq,gammasq,sigmasq_t,tausq_t,etasq_t){
  
  #Info from design
  design.mat <- design$swDsn#Design matrix
  M <- design$n.waves#Total number of waves
  J <- design$total.time#Total number of time periods
  
  RunningSum <- 0#Add up pieces here
  
  for(i in 1:M){#Sum over sequences
    
    #Design matrix X (fixed effects) For K=1:
    Xmat.i <- matrix(c(rep(1,length=J),0:(J-1),design.mat[i,]),nrow=J,ncol=3)
    #Expand X for the correct K:
    Xmat.i <- Xmat.i[rep(1:J,each=K),]#Keep it ordered by time
    
    #Design matrix Z for random time effect
    Zmat.time.i <- as.matrix(cbind(rep(1,length=J),diag(J)))
    Zmat.time.i <- Zmat.time.i[rep(1:J,each=K),]
    
    #Design matrix Z for random treatment effect
    Zmat.trt.i <- as.matrix(cbind(rep(1,length=J),design.mat[i,]))
    Zmat.trt.i <- Zmat.trt.i[rep(1:J,each=K),]
    
    #Matrix G (covariance matrix for random effects)
    G.time <- diag(c(tausq,rep(gammasq,length=J)))#Misspecified
    G.trt <- diag(c(tausq_t,etasq_t))#Correctly specified
    
    #Covariance matrix
    SigmaMat.time.i <- Zmat.time.i %*% G.time %*% t(Zmat.time.i)+sigmasq*diag(J*K)#Misspecified
    SigmaMat.trt.i <- Zmat.trt.i %*% G.trt %*% t(Zmat.trt.i)+sigmasq_t*diag(J*K)#Correctly specified
    #Element inside the sum
    Piece.i <- t(Xmat.i) %*% solve(SigmaMat.time.i) %*% SigmaMat.trt.i %*% solve(SigmaMat.time.i) %*% Xmat.i
    RunningSum <- RunningSum+Piece.i
    
  }
  
  return(RunningSum)#B (with N taken out already)
}


#####
#Function: calculate three different variances (false, true, and optimal) for a certain design and setting
#####

#This is the 'primary' function to implement these methods for case 1.  
#Arguments:
#design - design of the SWT, specified with swDsn with only one cluster per sequence.
#K - Number of observations per cluster per time period
#sigmasqTrue - the true value of sigma squared (residual variance from the correctly specified model)
#tausqTrue - the true value of tau squared (random intercept variance from the correctly specified model)
#etasqTrue - the true value of eta squared (random treatment variance from the correctly specified model)

#Function: for any design, take true parameter values and return three variances for Root 1 (if there was only one cluster per sequence)
GetVars_thisdesign_Case1_Root1 <- function(design,K,sigmasqTrue,tausqTrue,etasqTrue){
  M <- design$n.waves#Total number of sequences
  J <- design$total.time#Total number of time periods
  R <- sum(design$swDsn)#Total number of cluster x time periods on treatment (Sum of T's)
  S <- sum(rowSums(design$swDsn)^2)#Total number of squared cluster x time periods on treatment (Sum of T^2's)
  
  #Variance of treatment effect estimate from correctly specified model
  TrueModelBased <- GetCov_thisdesign_trt(design=design,K=K,sigmasq=sigmasqTrue,tausq=tausqTrue,etasq=etasqTrue)[3,3]
  #Covariance matrix for fixed effects from misspecified model, using Root 1 for the values of the random effect variances
  AInv <- GetCov_thisdesign_time(design=design,K=K,sigmasq=sigmasqTrue,tausq=etasqTrue*(S-R)/((J^2-J)*M)+tausqTrue,gammasq = etasqTrue*(J*R-S)/((J^2-J)*M))
  #Model-based variance of treatment effect estimate from misspecified model
  FalseModelBased <- AInv[3,3]
  #Meat for sandwich estimator from the misspecified model; requires both true variances and roots of variances from the misspecified model
  Meat <- GetMeat_thisdesign_falsetime(design=design,K=K,sigmasq=sigmasqTrue,tausq=etasqTrue*(S-R)/((J^2-J)*M)+tausqTrue,gammasq=etasqTrue*(J*R-S)/((J^2-J)*M),sigmasq_t=sigmasqTrue,tausq_t=tausqTrue,etasq_t=etasqTrue)
  #True variance of treatment effect from misspecified model
  SandwichBased <- (AInv %*% Meat %*% AInv)[3,3]
  
  #Returns false, true, and optimal variances
  VectorOfVars <- c(FalseModelBased,SandwichBased,TrueModelBased)
  names(VectorOfVars) <- c("False var","True var","Optimal var")
  return(VectorOfVars)
}



##############################################################################
#
#4. Treatment-fitted random time case (Case 2)
#
##############################################################################

#Note that these extra functions are necessary because Root 1 doesn't have an easy closed-form solution
#For roots that have closed-form solutions (Roots 3 and 4), you can skip all the numerical stuff and just plug in the formulas for the roots.
#Note: instead of T (the number of time points on treatment), we're using 'Tm' to represent T for sequence m

#####
#Function: Equation from system corresponding to taking derivatives with respect to sigma for a single cluster
#####

SigmaEq <- function(J,K,Tm,sigma_t,tau_t,gamma_t,sigma,tau,eta){
  value <- -(J*K)/(2)*(1)/(sigma^2)-0.5*(((J*K)/(sigma^2)+(1)/(tau^2))*((Tm*K)/(sigma^2)+(1)/(eta^2))-((Tm*K)/(sigma^2))^2)^(-1)*((-2*J*Tm*K^2)/((sigma^2)^3)-(J*K)/((sigma^2)^2*eta^2)-(K*Tm)/((sigma^2)^2*tau^2)+(2*K^2*Tm^2)/((sigma^2)^3))+ (1)/(2 *(sigma^2)^2)* J*K*(tau_t^2+gamma_t^2+sigma_t^2)-((1)/(sigma^2))^3*((1)/((((K*Tm)/(sigma^2)+(1)/(eta^2))*((J*K)/(sigma^2)+(1)/(tau^2))-((K*Tm)/(sigma^2))^2))*(((Tm*K)/(sigma^2)+(1)/(eta^2))*(J^2*K^2*tau_t^2+J*K^2* gamma_t^2)+2*((-Tm*K)/(sigma^2))*(J*Tm*K^2*tau_t^2+ Tm*K^2*gamma_t^2)+((J*K)/(sigma^2)+(1)/(tau^2))*(Tm^2*K^2*tau_t^2+Tm*K^2*gamma_t^2))+sigma_t^2* (1)/((((K*Tm)/(sigma^2)+(1)/(eta^2))*((J*K)/(sigma^2)+(1)/(tau^2))-((K*Tm)/(sigma^2))^2))*(2*(J-Tm)*K*(Tm*K)/(sigma^2)+J*K*(1)/(eta^2)+Tm*K*(1)/(tau^2)))+0.5* ((1)/(sigma^2))^4*((1)/((((K*Tm)/(sigma^2)+(1)/(eta^2))*((J*K)/(sigma^2)+(1)/(tau^2))-((K*Tm)/(sigma^2))^2)^2)*(((J-Tm)*(Tm^2*K^3)/((sigma^2)^2)+2*(J-Tm)*K*(Tm*K)/(sigma^2)*(1)/(eta^2)+J*K*(1)/((eta^2)^2))*(J^2*K^2*tau_t^2+J*K^2* gamma_t^2)+2*(Tm*K*(1)/(tau^2)*(1)/(eta^2)-(J-Tm)*(Tm^2*K^3)/((sigma^2)^2))*(J*Tm*K^2*tau_t^2+ Tm*K^2*gamma_t^2)+((J-Tm)*(J*Tm*K^3)/((sigma^2)^2)+2*(J-Tm)*K*(Tm*K)/(sigma^2)*(1)/(tau^2)+Tm*K*(1)/((tau^2)^2))*(Tm^2*K^2*tau_t^2+Tm*K^2*gamma_t^2))+sigma_t^2*(1)/((((K*Tm)/(sigma^2)+(1)/(eta^2))*((J*K)/(sigma^2)+(1)/(tau^2))-((K*Tm)/(sigma^2))^2)^2)*(2*(J-Tm)^2*(Tm^2*K^4)/((sigma^2)^2)+2*(J-Tm)*J*K^2*(Tm*K)/(sigma^2)*(1)/(eta^2)+J^2*K^2*(1)/((eta^2)^2) +2*Tm^2*K^2*(1)/(tau^2)*(1)/(eta^2)+2*(J-Tm)*K*(Tm^2*K^2)/(sigma^2)*(1)/(tau^2)+Tm^2*K^2*(1)/((tau^2)^2)))
  return(value*sigma)#Chain rule: derivative above is with respect to sigma^2 (since ==0, dropping constant)
}

#####
#Function: Equation from system corresponding to taking derivatives with respect to tau for a single cluster
#####

TauEq <- function(J,K,Tm,sigma_t,tau_t,gamma_t,sigma,tau,eta){
  value <- -0.5*(1)/(tau^2)+0.5*(((J*K)/(sigma^2)+(1)/(tau^2))*((Tm*K)/(sigma^2)+(1)/(eta^2))-((Tm*K)/(sigma^2))^2)^(-1)*((1)/(tau^2))^2*((K*Tm)/(sigma^2)+(1)/(eta^2))+0.5* ((1)/(sigma^2))^2*((1)/(tau^2))^2*(1)/((((K*Tm)/(sigma^2)+(1)/(eta^2))*((J*K)/(sigma^2)+(1)/(tau^2))-((K*Tm)/(sigma^2))^2)^2)*(((Tm*K)/(sigma^2)+(1)/(eta^2))^2*(J^2*K^2*tau_t^2+J*K^2* gamma_t^2)-2*((Tm*K)/(sigma^2)+(1)/(eta^2))*(Tm*K)/(sigma^2)*(J*Tm*K^2*tau_t^2+ Tm*K^2*gamma_t^2)+((Tm*K)/(sigma^2))^2*(Tm^2*K^2*tau_t^2+Tm*K^2*gamma_t^2)+sigma_t^2*(J*K*((Tm*K)/(sigma^2)+(1)/(eta^2))^2-((Tm*K)/(sigma^2)+(2)/(eta^2))*(Tm^2*K^2)/(sigma^2)))
  return(value*tau)#Chain rule: derivative above is with respect to tau^2 (since ==0, dropping constant)
}

#####
#Function: Equation from system corresponding to taking derivatives with respect to eta for a single cluster
#####

EtaEq <- function(J,K,Tm,sigma_t,tau_t,gamma_t,sigma,tau,eta){
  value <- -0.5/(eta^2)+0.5*(((J*K)/(sigma^2)+(1)/(tau^2))*((Tm*K)/(sigma^2)+(1)/(eta^2))-((Tm*K)/(sigma^2))^2)^(-1)*((1)/(eta^2))^2*((J*K)/(sigma^2)+(1)/(tau^2))+0.5* ((1)/(sigma^2))^2*((1)/(eta^2))^2*((1)/((((K*Tm)/(sigma^2)+(1)/(eta^2))*((J*K)/(sigma^2)+(1)/(tau^2))-((K*Tm)/(sigma^2))^2)^2))*(((Tm*K)/(sigma^2))^2*(J^2*K^2*tau_t^2+J*K^2* gamma_t^2)-2*(Tm*K)/(sigma^2)*((J*K)/(sigma^2)+(1)/(tau^2))*(J*Tm*K^2*tau_t^2+ Tm*K^2*gamma_t^2)+((J*K)/(sigma^2)+(1)/(tau^2))^2*(Tm^2*K^2*tau_t^2+Tm*K^2*gamma_t^2)+sigma_t^2*((-Tm^2*K^2)/(sigma^2)*((J*K)/(sigma^2)+(2)/(tau^2)) +  Tm*K*((J*K)/(sigma^2)+(1)/(tau^2))^2))
  return(value*eta)#Chain rule: derivative above is with respect to eta^2 (since ==0, dropping constant)
}

#####
#Function: Return values from entire system of equations
#####

#Arguments:
#design - design of the SWT, specified with swDsn with only one cluster per sequence.
#K - Number of observations per cluster per time period
#sigmasqTrue - the true value of sigma squared (residual variance from the correctly specified model)
#tausqTrue - the true value of tau squared (random intercept variance from the correctly specified model)
#gammasqTrue - the true value of gamma squared (random time variance from the correctly specified model)
#vector.of.mis.var.vals - the vector of values (sigmasqMis,tausqMis,etasqMis) that we are checking as candidates for the roots of the parameters from the misspecified model
#Returns:
#vector of the values from the three equations in the system of equations


full_system_of_equations <- function(design,K,sigmasqTrue,tausqTrue,gammasqTrue,vector.of.mis.var.vals){#vector.of.mis.var.vals=(sigmasqMis,tausqMis,etasqMis)
  J <- design$total.time#Total number of time periods
  T_vector <- rowSums(design$swDsn)#Number of cluster x time periods on treatment (T's) for each sequence
  #Converting from variances to SD's, since the system of equations is written in terms of SD's
  sigma_t <- sigmasqTrue^.5
  tau_t <- tausqTrue^.5
  gamma_t <- gammasqTrue^.5
  x <- vector.of.mis.var.vals^.5
  
  #Initialize three equation sums
  F1 <- 0
  F2 <- 0
  F3 <- 0
  
  for(i in 1:length(T_vector)){#Sum over all the sequences
    #Equation based on sigma
    F1 <- F1+SigmaEq(J=J,K=K,Tm=T_vector[i],sigma_t=sigma_t,tau_t=tau_t,gamma_t=gamma_t,sigma=x[1],tau=x[2],eta=x[3])
    #Equation based on tau
    F2 <- F2+TauEq(J=J,K=K,Tm=T_vector[i],sigma_t=sigma_t,tau_t=tau_t,gamma_t=gamma_t,sigma=x[1],tau=x[2],eta=x[3])
    #Equation based on eta
    F3 <- F3+EtaEq(J=J,K=K,Tm=T_vector[i],sigma_t=sigma_t,tau_t=tau_t,gamma_t=gamma_t,sigma=x[1],tau=x[2],eta=x[3])
  }
  
  #Return list of three values; the right vector.of.mis.var.vals will make this c(0,0,0) 
  list <- c(F1=F1,F2=F2,F3=F3)
  return(list)
}



#####
#Function: calculate 'meat' for a sandwich estimate, for a misspecified model with treatment random effects
#####

#This 'secondary' function is called by GetVars_thisdesign_Case2_Root1
#This returns the 'meat' (or B) for the Sandwich-form estimate of the true variance of a misspecified model 
#Requires both true variances and roots of variances from the misspecified model

#Consistent with the 'bread' functions, there is no N here.

GetMeat_thisdesign_falsetrt <- function(design,K,sigmasq,tausq,etasq,sigmasq_t,tausq_t,gammasq_t){
  
  #Info from design
  design.mat <- design$swDsn#Design matrix
  M <- design$n.waves#Total number of waves
  J <- design$total.time#Total number of time periods
  
  RunningSum <- 0#Add up pieces here
  
  for(i in 1:M){#Sum over sequences
    
    #Design matrix X (fixed effects):
    Xmat.i <- matrix(c(rep(1,length=J),0:(J-1),design.mat[i,]),nrow=J,ncol=3)
    Xmat.i <- Xmat.i[rep(1:J,each=K),]#Keep it ordered by time
    
    #Design matrix Z for random time effect
    Zmat.time.i <- as.matrix(cbind(rep(1,length=J),diag(J)))
    Zmat.time.i <- Zmat.time.i[rep(1:J,each=K),]
    
    #Design matrix Z for random treatment effect
    Zmat.trt.i <- as.matrix(cbind(rep(1,length=J),design.mat[i,]))
    Zmat.trt.i <- Zmat.trt.i[rep(1:J,each=K),]
    
    #Matrix G (covariance matrix for random effects)
    G.time <- diag(c(tausq_t,rep(gammasq_t,length=J)))#Correctly specified
    G.trt <- diag(c(tausq,etasq))#Misspecified
    
    #Covariance matrix
    SigmaMat.time.i <- Zmat.time.i %*% G.time %*% t(Zmat.time.i)+sigmasq_t*diag(J*K)#Correctly specified
    SigmaMat.trt.i <- Zmat.trt.i %*% G.trt %*% t(Zmat.trt.i)+sigmasq*diag(J*K)#Misspecified
    #Element inside the sum
    Piece.i <- t(Xmat.i) %*% solve(SigmaMat.trt.i) %*% SigmaMat.time.i %*% solve(SigmaMat.trt.i) %*% Xmat.i#Sum these up over all sequences, then invert and take [3,3]
    RunningSum <- RunningSum+Piece.i
    
  }
  
  return(RunningSum)#B (with N taken out already)
}

#####
#Function: calculate three different variances (false, true, and optimal) for a certain design and setting
#####

#This is the 'primary' function to implement these methods for case 2.  
#Arguments:
#design - design of the SWT, specified with swDsn with only one cluster per sequence.
#K - Number of observations per cluster per time period
#sigmasqTrue - the true value of sigma squared (residual variance from the correctly specified model)
#tausqTrue - the true value of tau squared (random intercept variance from the correctly specified model)
#gammasqTrue - the true value of gamma squared (random time variance from the correctly specified model)
#sigmasqMis - the root for sigma squared, found via numerical solution to system of equations (residual variance from the misspecified model)
#tausqMis - the root for tau squared, found via numerical solution to system of equations (random intercept variance from the misspecified model)
#etasqMis - the root for eta squared, found via numerical solution to system of equations (random treatment variance from the misspecified model)

#Function: for any design, take true parameter values and values of roots from misspecified model and return three variances for Root 1 (if there was only one cluster per sequence)

GetVars_thisdesign_Case2_Root1 <- function(design,K,sigmasqTrue,tausqTrue,gammasqTrue,sigmasqMis,tausqMis,etasqMis){
  M <- design$n.waves#Total number of sequences
  J <- design$total.time#Total number of time periods
  R <- sum(design$swDsn)#Total number of cluster x time periods on treatment (Sum of T's)
  S <- sum(rowSums(design$swDsn)^2)#Total number of squared cluster x time periods on treatment (Sum of T^2's)
  
  #Variance of treatment effect estimate from correctly specified model
  TrueModelBased <- GetCov_thisdesign_time(design=design,K=K,sigmasq=sigmasqTrue,tausq=tausqTrue,gammasq=gammasqTrue)[3,3]
  #Covariance matrix for fixed effects from misspecified model, using Root 1 for the values of the random effect variances
  AInv <- GetCov_thisdesign_trt(design=design,K=K,sigmasq=sigmasqMis,tausq=tausqMis,etasq = etasqMis)
  #Model-based variance of treatment effect estimate from misspecified model
  FalseModelBased <- AInv[3,3]
  #Meat for sandwich estimator from the misspecified model; requires both true variances and roots of variances from the misspecified model
  Meat <- GetMeat_thisdesign_falsetrt(design=design,K=K,sigmasq=sigmasqMis,tausq=tausqMis,etasq=etasqMis,sigmasq_t=sigmasqTrue,tausq_t=tausqTrue,gammasq_t=gammasqTrue)
  #True variance of treatment effect from misspecified model
  SandwichBased <- (AInv %*% Meat %*% AInv)[3,3]
  
  #Returns false, true, and optimal variances
  VectorOfVars <- c(FalseModelBased,SandwichBased,TrueModelBased)
  names(VectorOfVars) <- c("False var","True var","Optimal var")
  return(VectorOfVars)
}

##############################################################################
##############################################################################

##############################################################################
#
#5. Use and examples
#
##############################################################################

#####
#Time-fitted random treatment case (Case 1)
#####

#GetVars_thisdesign_Case1_Root1 is the 'primary' function to implement these methods for case 1.  
#Arguments:
#design - design of the SWT, specified with swDsn with only one cluster per sequence; that is, the 'clusters' argument should be a vector of just 0's and 1's.
#K - Number of observations per cluster per time period
#sigmasqTrue - the true value of sigma squared (residual variance from the correctly specified model)
#tausqTrue - the true value of tau squared (random intercept variance from the correctly specified model)
#etasqTrue - the true value of eta squared (random treatment variance from the correctly specified model)
#Returns:
#vector of false, true, and optimal variances of the treatment effect estimate, i.e. c(var_{m}(\hat\theta), var_{s}(\hat\theta), var_{m}(\hat\theta_t))

my.vars <- GetVars_thisdesign_Case1_Root1(design=swDsn(c(1,1,1)),K=2,sigmasqTrue=2,tausqTrue=.5,etasqTrue=.25)
my.vars#Recall that these are not scaled by N!  
#If you want to use these variances directly, you need to scale them by the number of clusters per sequence (see paper for formulas)
#But recall that this work is based on asymptotics, so applications should be considered carefully

my.vars[[1]]/my.vars[[2]]#Validity
my.vars[[3]]/my.vars[[2]]#Efficiency

#####
#Treatment-fitted random time case (Case 2)
#####

#Because we need to find the roots numerically, the code for Case 2 is a little more involved (for Root 1 - if you're interested in Roots 3 or 4, you can just use the closed-form roots).
#The steps to follow are outlined below:
#1. Find numerical roots for the parameters in the misspecified model
#2. Plug in your settings and the numerical roots to get the validity and efficiency of the treatment effect estimator


#####
#Step 1: Find numerical roots.
#I am presenting one way of doing this here, but there are many others (in R or other programs) which may be superior
#We are going to use a grid of lots of different starting values to try to find Root 1
#Take care: you could run into other roots (where one or more components are zero) here.  
#It may be especially hard to distinguish the roots when the true random effect variances are small, and more sophisticated techniques may be necessary
#I'm basing the rough range of this grid on the true variance settings
#If the roots you're finding look unreasonable, try a larger and/or finer grid of starting values!

#Build the grid of values to start on
grid_for_sigmasqMis <- seq(0.5,2,length=5)#I picked a range from half to double my true sigma squared
grid_for_tausqMis <- seq(0.0000001,1,length=10)#I picked a range from 0 to my true sigma squared
grid_for_etasqMis <- seq(0.0000001,1,length=10)#I picked a range from 0 to my true sigma squared
xstart <- as.matrix(expand.grid(grid_for_sigmasqMis,grid_for_tausqMis,grid_for_etasqMis),ncol=3)

#For convenience, list my settings of interest here:
my.design <- swDsn(c(1,1,1,1,1,1))
my.K <- 10
my.sigmasqTrue <- 1
my.tausqTrue <- 0.1
my.gammasqTrue <- 0.01


#Now, we look for roots
#Depending on the size of your grid, this may take a little time.  There are definitely faster ways to do this, especially if you have some decent starting values in mind.

#Arguments for full_system_of_equations:
#design - design of the SWT, specified with swDsn with only one cluster per sequence; that is, the 'clusters' argument should be a vector of just 0's and 1's.
#K - Number of observations per cluster per time period
#sigmasqTrue - the true value of sigma squared (residual variance from the correctly specified model)
#tausqTrue - the true value of tau squared (random intercept variance from the correctly specified model)
#gammasqTrue - the true value of gamma squared (random time variance from the correctly specified model)
#vector.of.mis.var.vals - the vector of values (sigmasqMis,tausqMis,etasqMis) that we are checking as candidates for the roots of the parameters from the misspecified model
#Returns:
#vector of the values from the three equations in the system of equations

roots <- searchZeros(xstart,function(x) full_system_of_equations(design=my.design,K=my.K,sigmasqTrue=my.sigmasqTrue,tausqTrue=my.tausqTrue,gammasqTrue=my.gammasqTrue,vector.of.mis.var.vals=x))$x
#Check to make sure this looks reasonable, and that it 'looks like' Root 1.
roots
#In this case, we found exactly one solution and it looks reasonable.
#If it had hit a different root, we would probably expect one of the components to be closer to zero (e.g. e-10)
#When in doubt, try expanding the range of the grid.  


#####
#Step 2: Plug everything in.

#Arguments for GetVars_thisdesign_Case2_Root1:
#design - design of the SWT, specified with swDsn with only one cluster per sequence; that is, the 'clusters' argument should be a vector of just 0's and 1's.
#K - Number of observations per cluster per time period
#sigmasqTrue - the true value of sigma squared (residual variance from the correctly specified model)
#tausqTrue - the true value of tau squared (random intercept variance from the correctly specified model)
#gammasqTrue - the true value of gamma squared (random time variance from the correctly specified model)
#sigmasqMis - the root for sigma squared, found via numerical solution to system of equations (residual variance from the misspecified model)
#tausqMis - the root for tau squared, found via numerical solution to system of equations (random intercept variance from the misspecified model)
#etasqMis - the root for eta squared, found via numerical solution to system of equations (random treatment variance from the misspecified model)
#Returns:
#vector of false, true, and optimal variances of the treatment effect estimate, i.e. c(var_{m}(\hat\theta), var_{s}(\hat\theta), var_{m}(\hat\theta_t))


#Recall that this design needs to match the design used to find the roots.
my.vars <- GetVars_thisdesign_Case2_Root1(design=my.design,K=my.K,sigmasqTrue=my.sigmasqTrue,tausqTrue=my.tausqTrue,gammasqTrue=my.gammasqTrue,sigmasqMis=roots[1],tausqMis=roots[2],etasqMis=roots[3])
my.vars#Recall that these are not scaled by N!  
#If you want to use these variances directly, you need to scale them by the number of clusters per sequence (see paper for formulas)
#But recall that this work is based on asymptotics, so applications should be considered carefully

my.vars[[1]]/my.vars[[2]]#Validity
my.vars[[3]]/my.vars[[2]]#Efficiency


