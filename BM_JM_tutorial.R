## BM_JM_Tutorial.R
##########################################################################################
##########################################################################################
## Description: A step-by-step implementation of BM-JM and the associated estimation    ##
## and dynamic prediction procedures described in "Bayesian multivariate joint modeling ##
## of longitudinal, recurrent, and competing risk terminal events in patients with      ##
## chronic kidney disease".                                                             ##
## Note that this tutorial is set for multivariate datasets with longitudinal, recurrent##
## and two competing risk terminal processes, but the codes can be easily modified to   ##
## accommodate more than two competing risk terminal processes.                         ##
##########################################################################################
##########################################################################################
## Functions implemented: 
## BM_JM_prepare_jags.R: R codes to prepare the multivariate data to be passed onto JAGS;
## BM_JM_jags.R: R file containing the JAGS model of the proposed BM-JM;
## BM_JM_prepare_prediction.R: R file containing functions to be used in dynamic predictions;
## BM_JM_measure_prediction_accuracy.R: R functions and codes to calculate dynamic AUC and BS 
##                                      values from the dynamic prediction output.
#############################################################################
## Tutorial Outline:
## 1. Load a multivariate training data set
## 2. Perform BM-JM estimation and inference, and visualize the estimation output 
 #    (BM_JM_prepare_jags.R, BM_JM_jags.R)
## 3. Load an independent multivariate testing data set
## 4. Perform BM-JM dynamic prediction algorithm on the testing data set  
 #    (BM_JM_prepare_prediction.R)
## 4. Measure the predictive accuracy of BM-JM  
 #    (BM_JM_measure_prediction_accuracy.R) 
## 5. Visualization of BM-JM dynamic prediction results
#############################################################################

################################################
## PACKAGES
################################################
# Install missing packages
list.of.packages <- c("rjags", "parallel", "nimble", "R2WinBUGS", "mvtnorm", "Matrix",
                      "MASS", "doParallel", "basicMCMCplots", "tidyverse",
                      "latex2exp", "timeROC", "survival", "pec")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) 

# Load libraries
library(rjags)
library(R2WinBUGS)
library(parallel)
require(R2WinBUGS) 
library(nimble)
library(mvtnorm)
library(Matrix)
library(MASS)
library(doParallel)
library(basicMCMCplots)
library(tidyverse)
library(latex2exp)
library(timeROC)
library(survival)
library(pec)

################################################
# Set working direction
################################################
setwd("~/Desktop/BM-JM codes")


################################################
# SPECIFY: the number of baseline (X) and event-varying (Z) covariates, and the number of maximum recurrent events per subject
################################################

ncX <- 1 ## number of baseline covariates (X)
ncZ <- 1 ## number of covariates for event-varying effects (Z)
max.nr <- 3 ## maximum number of recurrent events across all subjects

################################################
################################################
## 1. Load a multivariate training data set
################################################ 
################################################

# Load longitudinal outcome
load("dat.long.RData")
# Load recurrent outcomes
load("dat.recurrent.RData")
# Load competing-risk terminal outcomes
load("dat.terminal.RData")

################################################ 
################################################
## 2. Perform BM-JM estimation and inference, and visualize estimation output
################################################ 
################################################

################################################
# Prepare data for the following JAGS procedure
# Outputs are a list of `Data` that will be passed onto JAGS
################################################ 
source("BM_JM_prepare_jags.R")
  
################################################
# Create `text` file containing the JAGS model 
################################################
# Load `glm` mode
load.module("glm")
# Source JAGS model for BM-JM
source("BM_JM_jags.R")
# Create text file
filename <- file.path("BM_JM_jags.txt")
write.model(model, filename)
  
################################################
## Set parameters that will be tracked in the following JAGS 
## Notation follows the paper, only sigmasq.long stands for sigmasq.epsilon
################################################
params <- c("beta.l", "gamma", "phi.l", "beta.r", "eta.r0", "eta.r1", 
            "alpha1", "alpha",  "phi.r", "phi.r1", "beta.t.1", "beta.t.2", 
            "eta.t.1", "eta.t.2", "zeta.1", "zeta.2", "phi.t.1", "phi.t.2", 
            "sigmasq.long", "sigmasq.nu", "sigmasq.b0", "sigmasq.b1", "rho.01",
            "Bs.gammas.r", "Bs.gammas.t.1", "Bs.gammas.t.2")
  
###############################################
# Run JAGS model with MCMC chains running in parallel
###############################################
## 
# Number of MCMC chains 
n.chains <- 3

# Parallel chains
###############################################
# Create functions to run JAGS chain in parallel
## {{ 
## input: j representing a random scalar;
## To obtain reproducible simulation output, you must set the seed of the RJAGS random number generator. 
## This works differently than in base R. Instead of using set.seed(), you will specify a starting seed using inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = j)
## Note you can modify the numbers of adaption, update, running iteration, and thinning here.
## }}
## {{
## output: This produces a wrapper function that will be run in each core of your computer.
## }}
coda.samples.wrapper <- function(j)
  {
    model.fit <- jags.model(file = "BM_JM_jags.txt",
                            inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = j), ## to make sure each chain starts from a different set of initial values
                            data = Data, n.chains = 1, n.adapt = 2000)
    update(model.fit, 3000)
    
    coda.samples(model.fit, params,  n.iter = 10000, thin = 5)
  }
  
# Run JAGS model in parallel
post.samples <- mclapply(1:n.chains, coda.samples.wrapper,  mc.cores = n.chains)

# Transform BM-JM estimation output into a mcmc list
for(ii in 1:length(post.samples)){ 
  post.samples[[ii]] <- post.samples[[ii]][[1]]
}
class(post.samples) <- "mcmc.list"

# Visualize BM-JM estimation results
summary(post.samples)

###########################################################################
###########################################################################
# 3. Load an independent testing data set, which also includes longitudinal, 
# recurrent, and two competing-risk terminal outcomes
###########################################################################
###########################################################################
load("pred.subject.RData")

#################################################################### 
####################################################################
# 4. Perform BM-JM dynamic prediction procedure on the testing data set
#################################################################### 
####################################################################

####################################################################
# Source functions which are needed in the dynamic prediction procedure
#################################################################### 
source("BM_JM_prepare_prediction_3.R")


################################################################# 
# Perform BM-JM dynamic prediction procedure on the testing data set
################################################################# 

# Number of subjects in the testing data set
n.sample <- nrow(pred.subject$dat.terminal)
 
# Check patient IDs in the testing data set
nd.id <- pred.subject$dat.terminal$subj.id

# Create a list that stores dynamic prediction outputs
PredOut <- list()

# Do dynamic predictions
for (k in 1:n.sample){
  PredOut[[k]] <- joint.predict(T.start = c(0, 0.2), # Prediction time points (t in the paper)
                                tt = c(0.2, 0.4, 0.6, 0.8, 1), # Prediction time windows (\Delta t in the paper)
                                pred.subject = pred.subject,# Key information from the testing data set
                                output = post.samples, # BM-JM estimation output
                                M = 200, # Number of Monte Carlo iterations in dynamic predictions (H in the paper)
                                sub.id = nd.id[k] # Subject to be predicted 
                               )
  }


################################################
################################################
# 5. Measure the accuracy of BM-JM prediction 
################################################ 
################################################
source("BM_JM_measure_prediction_accuracy_2.R")

# Dynamic AUC and BS values in different predition time t and prediction time window \Delta t
AUCst.1
AUCst.2
BrierS.s.1
BrierS.s.2


#################################################################
#################################################################
## 6. Visualize BM-JM prediction output for subjects of interest
################################################################# 
#################################################################
T.start = c(0, 0.2)
tt = c(0.2, 0.4, 0.6, 0.8, 1)

# For example, we want to see the CIF for the 6th subject in the testing dataset
# For T.start = 0
d1 <- PredOut[[6]]$cif.out[[1]]$cond.cif.t1
d1[which(d1 >=1,arr.ind=TRUE)] <- 1
d1[which(d1 <=0,arr.ind=TRUE)] <- 0
d2 <- PredOut[[6]]$cif.out[[1]]$cond.cif.t2
d2[which(d2 >=1,arr.ind=TRUE)] <- 1
d2[which(d2 <=0,arr.ind=TRUE)] <- 0
# Medians
Tstart0.mean_1 <- apply(d1, 2, median)
Tstart0.mean_2 <- apply(d2, 2, median)
# Percentiles
Tstart0.q2.5_1 <- apply(d1,2, quantile, probs = 0.025 , na.rm = TRUE)
Tstart0.q2.5_2 <- apply(d2,2, quantile, probs = 0.025 , na.rm = TRUE)
Tstart0.q97.5_1 <-apply(d1,2, quantile, probs = 0.975 , na.rm = TRUE)
Tstart0.q97.5_2 <-apply(d2,2, quantile, probs = 0.975 , na.rm = TRUE)


# Plot the estimated CIF and 95% credibel intervals for the selected subject
par(mfrow = c(1,2))
xtime = c(0.2, 0.4, 0.6, 0.8, 1)
plot(x = xtime, y = Tstart0.mean_1, type="b", pch=19,  xlim = c(0, 1),
     main = "(A) 1st Competing Risk: t = 0",lwd = 2, cex.main=1.3, 
     xlab = "", ylab = "", cex.axis=1.1, col = "skyblue4", ylim = c(0,1))
title(xlab = "Time",ylab=TeX(r'($\pi^{(1)}_{p}(t,s)$)', bold=F), line=2, cex.lab=1.2)
segments(xtime, Tstart0.q2.5_1, xtime, Tstart0.q97.5_1, col = "skyblue4", lwd = 2)
segments(xtime-0.02, Tstart0.q2.5_1, xtime+0.02, Tstart0.q2.5_1, col = "skyblue4", lwd = 2)
segments(xtime-0.02, Tstart0.q97.5_1, xtime+0.02, Tstart0.q97.5_1, col = "skyblue4", lwd = 2)

xtime = c(0.2, 0.4, 0.6, 0.8, 1)
plot(x = xtime, y = Tstart0.mean_2, type="b", pch=19,  xlim = c(0, 1),
     main = "(B) 2nd Competing Risk: t = 0",lwd = 2, cex.main=1.3, 
     xlab = "", ylab = "", cex.axis=1.1, col = "indianred4", ylim = c(0,1))
title(xlab = "Time",ylab=TeX(r'($\pi^{(2)}_{p}(t,s)$)', bold=F), line=2, cex.lab=1.2)
segments(xtime, Tstart0.q2.5_2, xtime, Tstart0.q97.5_2, col = "indianred4", lwd = 2)
segments(xtime-0.02, Tstart0.q2.5_2, xtime+0.02, Tstart0.q2.5_2, col = "indianred4", lwd = 2)
segments(xtime-0.02, Tstart0.q97.5_2, xtime+0.02, Tstart0.q97.5_2, col = "indianred4", lwd = 2)
