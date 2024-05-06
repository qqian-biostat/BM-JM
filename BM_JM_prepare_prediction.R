##########################################################################################
##########################################################################################
## This R file contains all the functions needed in BM-JM dynamic prediction procedures ##
##########################################################################################
##########################################################################################


######################################################
### Function to generate baseline splines 
######################################################

## {{input: x: time points, vector;
##          sk: nodes for Gauss quadrature, vector. }}
## 
## {{output: a list of Lam: B spline matrix for baseline;
##                     Lam.s: B spline matrix for baseline in 15-point quadrature.}}
Bspline.tool <- function(x, sk){
  
  tlo <- 0
  nknots = 10 ## 10 knots are used in each baseline function fit
  tlo = 0
  nseg = nknots-2
  bdeg = 2
  
  P <- x/2
  st <- outer(P, sk + 1) 
  
  thi = floor(max(x)) + 1
  Lam = bbase1(x, tlo, thi, nseg, bdeg)  #for baseline
  Lam.s = bbase1(c(t(st)), tlo, thi, nseg, bdeg)  #for baseline in 15-point quadrature
  DDal <- diag(ncol(Lam)) 
  return(list(Lam=Lam, Lam.s=Lam.s))
}



######################################################################
# Function to give the density of a multivariate t distribution
# Defined by Dr. Dimitris Rizopoulos (taken from JMbayes package)
# Related paper: "Description Shared parameter models for the joint modeling of longitudinal and time-to-event
# data using MCMC; Dimitris Rizopoulos (2016) <doi:10.18637/jss.v072.i07>".
# https://github.com/drizopoulos/JMbayes/blob/master/R
######################################################################

## {{input: x: a numeric matrix or vector;
##          mu: mean values of x, a numeric vector or matrix representing the mean vector of the multivariate t distribution;
##          Sigma or invSigma: positive-definite covariance matrix of multivariate t distribution, must provide one of them;
##          df: a scalar, degree of freedom in multivariate t distribution;
##          log: logical, TRUE or FALSE. If log=TRUE, then the logarithm of the density is returned;
##          prop: TRUE or FALSE}}
## 
## {{output: density of a multivariate t distribution.}}

dmvt <-
  function (x, mu, Sigma = NULL, invSigma = NULL, df, log = FALSE, prop = TRUE) {
    if (!is.numeric(x)) 
      stop("'x' must be a numeric matrix or vector")
    if (!is.matrix(x)) 
      x <- rbind(x)
    p <- length(mu)
    if (is.null(Sigma) && is.null(invSigma))
      stop("'Sigma' or 'invSigma' must be given.")
    if (!is.null(Sigma)) {
      if (is.list(Sigma)) {
        ev <- Sigma$values
        evec <- Sigma$vectors
      } else {
        ed <- eigen(Sigma, symmetric = TRUE)
        ev <- ed$values
        evec <- ed$vectors
      }
      if (!all(ev >= -1e-06 * abs(ev[1]))) 
        stop("'Sigma' is not positive definite")
      invSigma <- evec %*% (t(evec)/ev)
      if (!prop)
        logdetSigma <- sum(log(ev))
    } else {
      if (!prop)
        logdetSigma <- c(-determinant(invSigma)$modulus)
    }
    ss <- x - rep(mu, each = nrow(x))
    quad <- rowSums((ss %*% invSigma) * ss)/df
    if (!prop)
      fact <- lgamma((df + p)/2) - lgamma(df/2) - 0.5 * (p * (log(pi) + log(df)) + logdetSigma)
    if (log) {
      if (!prop) as.vector(fact - 0.5 * (df + p) * log(1 + quad)) else as.vector(- 0.5 * (df + p) * log(1 + quad))
    } else {
      if (!prop) as.vector(exp(fact) * ((1 + quad)^(-(df + p)/2))) else as.vector(((1 + quad)^(-(df + p)/2)))
    }
  }

######################################################################
# Function to generate random deviates from multivariate t distribution
# Defined by Dr. Dimitris Rizopoulos (taken from JMbayes package)
# Related paper: "Description Shared parameter models for the joint modeling of longitudinal and time-to-event
# data using MCMC; Dimitris Rizopoulos (2016) <doi:10.18637/jss.v072.i07>".
# https://github.com/drizopoulos/JMbayes/blob/master/R
######################################################################

## {{input: n: a numeric scalar denoting the number of random draws;
##          mu: a numeric vector or matrix representing the mean vector of the multivariate t distribution;
##          Sigma: positive-definite covariance matrix of multivariate t distribution;
##          df: a scalar, degree of freedom in multivariate t distribution;
##          }}
## 
## {{output: random deviates from multivariate t distribution.}}

rmvt <- function (n, mu, Sigma, df) {
  p <- length(mu)
  if (is.list(Sigma)) {
    ev <- Sigma$values
    evec <- Sigma$vectors
  } else {
    ed <- eigen(Sigma, symmetric = TRUE)
    ev <- ed$values
    evec <- ed$vectors
  }
  X <- drop(mu) + tcrossprod(evec * rep(sqrt(pmax(ev, 0)), each = p), 
                             matrix(rnorm(n * p), n)) / rep(sqrt(rchisq(n, df)/df), each = p)
  if (n == 1L) drop(X) else t.default(X)
}


############################################################
# Function to prepare a multivariate data set into a list that 
# includes key informations from the data set.
############################################################
## {{input: multivariate data, including a longitudinal outcome (dat.long), 
##          recurrent events outcome (dat.recurrent), and competing-risk
##          terminal events outcome (dat.terminal).}}
## {{output: a list of `Data` that includes key information, including:
##           wk, ## weights for Gauss quadrature to be used in approximating survival function
##           sk, ## nodes for Gauss quadrature to be used in approximating survival function
##           y, ## longitudinal outcome
##           times, ## time points for longitudinal outcome 
##           x1, ## secondary covariates for new patient
##           z1, ## primary covariates for new patient
##           n.r, ## number of recurrent events for new patient
##           event.r, ## recurrent event status for new patient
##           T1.r, ## recurrent gap times for new patient
##           T1.t, # terminal event times for new patient
##           Lam.sr, ## B spline matrix in 15-point quadrature for baseline function in recurrent submodel
##           Lam.r, ## B spline matrix for baseline function in recurrent submodel
##           Lam.st, ## B spline matrix in 15-point quadrature for baseline function in terminal submodel
##           Lam.t ## B spline matrix for baseline function in terminal submodel. }}
prepare.data <- function(dat.long, dat.recurrent, dat.terminal){
  ## First, gather information on the recurrent events
  # Time.r0 <- dat.recurrent$T.rec0 ## these are the event times for all subjects at all recurrent events (stacked)
  Time.r <- dat.recurrent$T.gap ## these are the gap times for all subjects at all recurrent events (stacked)
  
  event.r <- dat.recurrent$event.r ## these are the event status for all subjects at all recurrent event times (stacked) 
  
  subjs <- dat.recurrent$subj.id ## subjects with recurrent data
  
  n <- length(unique(subjs)) ## n: number of subjects 
  
  N.r <- dim(dat.recurrent)[1] ## N.r: dimension of the recurrent event data
  
  ## data.id has one row per subject
  data.id <- as.data.frame(dat.recurrent[!duplicated(subjs), ])
  
  x1 <- data.id[, 5] ## secondary covariates
  
  z1 <- data.id[,6]  ## primary covariates 
  
  ## Number of observations (rows) for each subject in the recurrent events data, denoted as n.r
  nts.r <- as.vector(table(x = subjs))
  
  n.r <- c(0, nts.r) 
  
  #############################################
  ## Second, gather information on the terminal events
  
  Time.t <- dat.terminal$T.terminal ## these are the terminal event times for all subjects (stacked)
 
  #############################################
  ## Weights and nodes for Gauss quadrature to be used in approximation of survival function -- both in terminal and recurrent events
  
  wk <- c(0.0630920926299786, 0.140653259715526, 0.190350578064785,
          0.209482141084728, 0.190350578064785, 0.140653259715526,
          0.0630920926299786, 0.0229353220105292, 0.10479001032225,
          0.169004726639268, 0.204432940075299, 0.204432940075299,
          0.169004726639268, 0.10479001032225, 0.0229353220105292)
  sk <- c(-0.949107912342758, -0.741531185599394, -0.405845151377397,
          0, 0.405845151377397, 0.741531185599394, 0.949107912342758,
          -0.991455371120813, -0.864864423359769, -0.586087235467691,
          -0.207784955007898, 0.207784955007898, 0.586087235467691,
          0.864864423359769, 0.991455371120813)
  
  ordsk <- order(sk)
  sk <- sk[ordsk]
  wk <- wk[ordsk]
  
  K <- length(sk) ## number of quadrature points
  
  ## For recurrent
  P.r <- Time.r/2 # Use the recurrent gap times
  st.r <- outer(P.r, sk + 1) ## dimension: nxlength(sk)
  
  ## For terminal
  P.t <- Time.t/2
  st.t <- outer(P.t, sk + 1) ## dimension:nxlength(sk)
  
  
  ############################################
  ## Define parameters for P-splines used in fitting the baseline hazard functions
  
  nknots = 10 ## 10 knots are used in each baseline function fit
  tlo = 0
  nseg = nknots-2
  bdeg = 2
  
  ############################################ 
  ## P-spline terms to estimate the baseline hazard function in the recurrent events model 
  
  thi = floor(max(Time.r)) + 1    # Use the recurrent gap times
  
  B0 = bbase1(Time.r, tlo, thi, nseg, bdeg) # Use the recurrent gap times
  Lam.r <- B0 ## dimension: nxnknots, one row for each subject -- this is for the hazard function
  
  
  B0 = bbase1(c(t(st.r)), tlo, thi, nseg, bdeg)
  Lam.sr <- B0 
  
  ############################################
  ## P-spline terms to estimate the baseline hazard function in the terminal event model 
  
  thi = floor(max(Time.t)) + 1
  
  B0 = bbase1(Time.t, tlo, thi, nseg, bdeg)
  Lam.t <- B0 
  
  
  B0 = bbase1(c(t(st.t)), tlo, thi, nseg, bdeg)
  Lam.st <- B0 
  
  ############################################
  ## Lastly, gather information on the longitudinal outcome

  times <- dat.long$t ## repeated measurement times
  
  ############################################
  # Data list that include key information from the multivariate data set
  
  Data <- list(wk = wk, ## weights used in approximating survival function
               sk = sk,
               y = dat.long$Y, ## longitudinal outcome
               times = times, 
               x1 = as.matrix(x1),
               z1 = as.matrix(z1),
               n.r = n.r, 
               event.r = event.r,
               T1.r = Time.r, ## recurrent gap times
               T1.t = Time.t, # terminal event time in terminal submodel
               Lam.sr = Lam.sr,
               Lam.r = Lam.r,
               Lam.st = Lam.st,
               Lam.t = Lam.t
  )
  
  return(Data)
}


#########################################################################   
# Function to define the log-posterior distribution of the random effects
#########################################################################  
## {{input: u, ## vector of values of random effects in calculating log posterior of;
##          new.data, ## a list from the output of `prepare.data` function containing information for testing subjects;
##          output, ## an MCMC list from BM-JM estimation JAGS output;
##          theta, ## a list of values of random effect for the current iteration;
##          nd.id ## subject id to be predicted in the current iteration, scalar.}}
## {{output: scalar representing the value of log posterior distribution}}

log.posterior.b <- function(u, # value of u to calculate log posterior of
                            new.data, # testing subjects data
                            output, # BM-JM estimation output 
                            theta, # value of theta for this iteration
                            nd.id){ # which subject id in the data
  # Longitudinal portion
  mu.y <- c(new.data$x1 * theta$beta.l) + c(new.data$z1 * theta$phi.l) + 
    theta$gamma * new.data$times + u[1] + u[2] * new.data$times
  log.y <- dnorm(new.data$y, mean = mu.y, sd = sqrt(theta$sigmasq.long), log = TRUE)
  log.p.y <- sum(log.y)
  sigma.b <- matrix(0, nrow = 2, ncol = 2)
  sigma.b[1,1] <- theta$sigmasq.b0
  sigma.b[2,2] <- theta$sigmasq.b1
  sigma.b[1,2] <- sigma.b[2,1] <- theta$rho.01
  log.p.b <- dmvnorm(matrix(c(u[1], u[2]), ncol = 2), mean=c(0,0), sigma = sigma.b, log=TRUE)
  log.p.nu <- dnorm(u[3], mean = 0, sd = sqrt(theta$sigmasq.nu), log = TRUE)
  
  # Recurrent portion
  surv.r <- matrix(nrow = sum(new.data$n.r), ncol = K)
  lin.pred.rs <- matrix(nrow = sum(new.data$n.r), ncol = K)
  log.h0.rs <- matrix(nrow = sum(new.data$n.r), ncol = K)
  surv.r <- matrix(nrow = sum(new.data$n.r), ncol = K)
  log.surv.r <- log.h0.r <- haz.r <- log.rec <- vector(length = sum(new.data$n.r))
  
  for (j in 1:sum(new.data$n.r)){
    
    log.h0.r[j] <- new.data$Lam.r[j, ]%*%theta$Bs.gammas.r
    
    haz.r[j] <- exp(log.h0.r[j] + new.data$x1*theta$beta.r +  theta$phi.r1[j]* new.data$z1 + 
                      theta$eta.r0*u[1] + theta$eta.r1*u[2] + u[3] + sum(theta$alpha1[1:j]))
    
    for(k in 1:K)
    {
      log.h0.rs[j, k] <- new.data$Lam.sr[K * (j - 1) + k, ]%*%theta$Bs.gammas.r
      lin.pred.rs[j, k] <- exp(theta$beta.r*new.data$x1 + theta$phi.r1[j]*new.data$z1 +
                                 sum(theta$alpha1[1:j]) + theta$eta.r0*u[1] + 
                                 theta$eta.r1*u[2] + u[3])
      surv.r[j, k] <- exp(log.h0.rs[j, k]) * lin.pred.rs[j, k]
     }
    log.surv.r[j] <-  (-new.data$T1.r[j]/2*new.data$wk%*%surv.r[j,])
    log.rec[j] <- new.data$event.r[j] * log(haz.r[j]) + log.surv.r[j]
    }
    log.surv.r = sum(log.rec)
  
  
  # Survival portion
   surv.t.1 <- surv.t.2 <- vector(length = K)
   lin.pred.ts.1 <- lin.pred.ts.2 <- vector(length = K)
   log.h0.ts.1 <- log.h0.ts.2 <- vector(length = K)
   mu.ts <- vector(length = K)
   surv.t.1 <- surv.t.2 <- vector(length = K)
  
   for(k in 1:K)
   {
    log.h0.ts.1[k] <- new.data$Lam.st[k, ]%*%theta$Bs.gammas.t.1
    log.h0.ts.2[k] <- new.data$Lam.st[k, ]%*%theta$Bs.gammas.t.2
    mu.ts[k] <- new.data$x1*theta$beta.l+new.data$z1*theta$phi.l + 
                theta$gamma*(new.data$T1.t/2*(sk[k]+1))+ u[1] + 
                u[2]*(new.data$T1.t/2*(sk[k]+1)) 
    lin.pred.ts.1[k] <- exp(theta$beta.t.1*new.data$x1 + theta$phi.t.1*new.data$z1 +
                              theta$eta.t.1 * mu.ts[k] + theta$zeta.1 * u[3])
    lin.pred.ts.2[k] <- exp(theta$beta.t.2*new.data$x1 + theta$phi.t.2*new.data$z1 +
                              theta$eta.t.2 * mu.ts[k] + theta$zeta.2 * u[3])
    surv.t.1[k] <- exp(log.h0.ts.1[k]) * lin.pred.ts.1[k]
    surv.t.2[k] <- exp(log.h0.ts.2[k]) * lin.pred.ts.2[k]
   }
   log.surv.t.1 <-  (-new.data$T1.t/2*new.data$wk%*%surv.t.1)
   log.surv.t.2 <-  (-new.data$T1.t/2*new.data$wk%*%surv.t.2)
   log.surv.t <- log.surv.t.1 + log.surv.t.2
  
  # Return posterior distribution of the random effects
   return(sum(c(log.p.y,log.surv.t, log.surv.r,log.p.b, log.p.nu)))
 }

################################################################################################## 
# Function to dynamically predict cumulative incidence function for subjects in the testing data set
################################################################################################## 
## {{input: T.start, ## a vector of prediction times, t in the paper;
##          tt, ## a vector of prediction time windows, \Delta t in the paper;
##          pred.subject ## a list from the output of `prepare.data` function containing information for testing subjects;  
##          output, ## an MCMC list from BM-JM estimation JAGS output;
##          M, ## number of Monte Carlo iterations used in dynamic prediction, scalar, H in the paper;
##          sub.id, ## subject id to be predicted, scalar.
## }}
## {{output: list of lists, with the first level of list is different prediction times, 
##           and the second level of list contains two sublist:
##           cond.cif.t1 and cond.cif.t2;
##           Within each sublist is a matrix (#prediction time windows*M):
##                       with columns representing different prediction time windows, 
##                       and rows representing different MC iterations. }}

joint.predict <- function(T.start, # prediction time, t in the paper
                          tt, # prediction time window, \Delta t in the paper
                          pred.subject, # key information from the testing data set
                          output, # BM-JM estimation output
                          M = 200, # Number of Monte Carlo iterations in dynamic prediction, H in the paper
                          sub.id){ # subject id to be predicted
  
  # Clean raw data for use in JAGS, sample random thetas for use in prediction.
  dat.long.run <- pred.subject$dat.long[pred.subject$dat.long$subj_id == sub.id,]
  dat.recurrent.run <- pred.subject$dat.recurrent[pred.subject$dat.recurrent$subj.id == sub.id,]
  dat.terminal.run <- pred.subject$dat.terminal[pred.subject$dat.terminal$subj.id == sub.id,]
  new.data <- prepare.data(dat.long.run, dat.recurrent.run, dat.terminal.run)
  thetas <- mcmc(do.call(rbind, output))
  thetas <- thetas[sample(nrow(thetas), nrow(thetas)),]
  n.thetas <- min(nrow(thetas), M) 
 
  # Prepare empirical Bayes estimates and associated variances
  modes.u <- matrix(0, nrow = n.thetas, ncol = 3) # point estimates go here
  u.old <- u.new <- modes.u
  invVars.u <- Vars.u <- vector("list", n.thetas) # observed information goes here
  temp.theta  <- vector("list", n.thetas)
  
  #Calculate Empirical Bayes Estimates
  for (i in seq_len(n.thetas)) {
    
    #Aggregate the parameters for use in the prediction
    temp.theta[[i]]$Bs.gammas.r <- thetas[i,][startsWith(names(thetas[i,]),"Bs.gammas.r")]
    temp.theta[[i]]$Bs.gammas.t.1 <- thetas[i,][startsWith(names(thetas[i,]),"Bs.gammas.t.1")]
    temp.theta[[i]]$Bs.gammas.t.2 <- thetas[i,][startsWith(names(thetas[i,]),"Bs.gammas.t.2")]
    temp.theta[[i]]$alpha <- thetas[i,][startsWith(names(thetas[i,]),"alpha")]
    temp.theta[[i]]$alpha1 <- thetas[i,][startsWith(names(thetas[i,]),"alpha1")]
    temp.theta[[i]]$beta.l <- thetas[i,][startsWith(names(thetas[i,]),"beta.l")]
    temp.theta[[i]]$beta.r <- thetas[i,][startsWith(names(thetas[i,]),"beta.r")]
    temp.theta[[i]]$beta.t.1 <- thetas[i,][startsWith(names(thetas[i,]),"beta.t.1")]
    temp.theta[[i]]$beta.t.2 <- thetas[i,][startsWith(names(thetas[i,]),"beta.t.2")]
    temp.theta[[i]]$eta.r0 <- thetas[i,][startsWith(names(thetas[i,]),"eta.r0")]
    temp.theta[[i]]$eta.r1 <- thetas[i,][startsWith(names(thetas[i,]),"eta.r1")]
    temp.theta[[i]]$eta.t.1  <- thetas[i,][startsWith(names(thetas[i,]),"eta.t.1")]
    temp.theta[[i]]$eta.t.2  <- thetas[i,][startsWith(names(thetas[i,]),"eta.t.2")]
    temp.theta[[i]]$gamma  <- thetas[i,][startsWith(names(thetas[i,]),"gamma")]
    temp.theta[[i]]$phi.l <- thetas[i,][startsWith(names(thetas[i,]),"phi.l")]
    temp.theta[[i]]$phi.r <- thetas[i,][startsWith(names(thetas[i,]),"phi.r")]
    temp.theta[[i]]$phi.r1 <- thetas[i,][startsWith(names(thetas[i,]),"phi.r1")]
    temp.theta[[i]]$phi.t.1 <- thetas[i,][startsWith(names(thetas[i,]),"phi.t.1")]
    temp.theta[[i]]$phi.t.2 <- thetas[i,][startsWith(names(thetas[i,]),"phi.t.2")]
    temp.theta[[i]]$rho.01  <- thetas[i,][startsWith(names(thetas[i,]),"rho.01")]
    temp.theta[[i]]$sigmasq.b0 <- thetas[i,][startsWith(names(thetas[i,]),"sigmasq.b0")]
    temp.theta[[i]]$sigmasq.b1 <- thetas[i,][startsWith(names(thetas[i,]),"sigmasq.b1")]
    temp.theta[[i]]$sigmasq.long <- thetas[i,][startsWith(names(thetas[i,]),"sigmasq.long")]
    temp.theta[[i]]$sigmasq.nu <- thetas[i,][startsWith(names(thetas[i,]),"sigmasq.nu")]
    temp.theta[[i]]$zeta.1 <- thetas[i,][startsWith(names(thetas[i,]),"zeta.1")]
    temp.theta[[i]]$zeta.2 <- thetas[i,][startsWith(names(thetas[i,]),"zeta.2")]
    
    # Define the function to minimize as the negative log posterior distribution
    ff <- function(u=u,new.data,output,theta,nd.id = sub.id){
      -log.posterior.b(u=u,new.data=new.data,output=output,theta=temp.theta[[i]],nd.id=sub.id)
    }
    
    # Find our estimated u value (and associated Hessian) via BFGS
    opt <- optim(c(0,0,0),fn=ff,new.data=new.data,output=output,theta=temp.theta[[i]], hessian=TRUE, method="BFGS")
    
    # Point estimate
    modes.u[i, ] <- opt$par
    
    # Get estimated information matrix (scale to improve acceptance rate)
    invVars.u[[i]] <- opt$hessian/1.6 # 1.6 is just a scale parameter here
    Vars.u[[i]] <- 1.6*solve(opt$hessian)   
  }
  
  # Generate the full set of independent proposals from multivariate-t for MH
  success.rate <- matrix(FALSE, M, n.thetas)
  u.old <- u.new <- modes.u
  
  proposed.u <- mapply(rmvt, mu = split(modes.u, 1:nrow(modes.u)), Sigma = Vars.u, 
                       MoreArgs = list(n = M, df = 4), SIMPLIFY = FALSE)
  dmvt.proposed <- mapply(dmvt, x = proposed.u, mu = split(modes.u, 1:nrow(modes.u)),
                          Sigma = Vars.u, MoreArgs = list(df = 4, log = TRUE), 
                          SIMPLIFY = FALSE)
  
  #For given parameter values theta
  for (m in 1:M){
    #Draw new random effect values from posterior using MH
    for (i in seq(n.thetas)) {
      
      p.u <- proposed.u[[i]][m,]
      dmvt.old <- dmvt(x = u.old[i,], mu = modes.u[i,], Sigma = Vars.u[[i]], df = 4, log = TRUE)
      dmvt.prop <- dmvt.proposed[[i]][m]
      a <- min(exp(log.posterior.b(u=p.u, new.data=new.data, output=output,
                                   theta=temp.theta[[i]],nd.id=sub.id) + 
                     dmvt.old - 
                     log.posterior.b(u=u.old[i,], new.data=new.data, output=output,
                                     theta=temp.theta[[i]],nd.id=sub.id) - 
                     dmvt.prop), 1)
      ind <- runif(1) <= a
      success.rate[m, i] <- ind
      if (!is.na(ind) && ind)
        u.new[i,] <- p.u
    }
    u.old <- u.new
  }
  
  # Estimate survival at given times for parameter values theta
  # and simulated random effects u using the 15-point quadratures
  surv.function <- function(times, u, temp.theta){
    surv.probs.1 <- vector(length = length(times))
    surv.probs.2 <- vector(length = length(times)) 
    x1.new <- new.data$x1
    z1.new <- new.data$z1
    sk <- new.data$sk
    
    for (m in times){
      # for (j in seq(B)){
      
      log.h0.ts.1 <- vector(length=K)
      log.h0.ts.2 <- vector(length=K)
      surv.1 <- vector(length=K)
      surv.2 <- vector(length=K)
      mu.surv <- vector(length=K)
      lin.pred.1 <- vector(length=K)
      lin.pred.2 <- vector(length=K)
      
      ## survival function
      for(k in 1:K)
      {
        log.h0.ts.1[k] <- Bspline.tool(m, sk)$Lam.s[k, ]%*%temp.theta$Bs.gammas.t.1[]
        log.h0.ts.2[k] <- Bspline.tool(m, sk)$Lam.s[k, ]%*%temp.theta$Bs.gammas.t.2[]
        mu.surv[k] <-  x1.new*temp.theta$beta.l + z1.new*temp.theta$phi.l +
          temp.theta$gamma*(m/2*(sk[k]+1))+ u[1]+ u[2]*(m/2*(sk[k]+1))
        
        lin.pred.1[k] <- exp(x1.new*temp.theta$beta.t.1+z1.new*temp.theta$phi.t.1 +
                               temp.theta$eta.t.1*mu.surv[k]+
                               temp.theta$zeta.1*u[3])
        lin.pred.2[k] <- exp(x1.new*temp.theta$beta.t.2+z1.new*temp.theta$phi.t.2 +
                               temp.theta$eta.t.2*mu.surv[k]+
                               temp.theta$zeta.2*u[3])
        surv.1[k] <- exp(log.h0.ts.1[k])*lin.pred.1[k]
        surv.2[k] <- exp(log.h0.ts.2[k])*lin.pred.2[k]
      }
      surv.probs.1[which(times==m)] <-  (exp(-m/2*(wk%*%surv.1)))
      surv.probs.2[which(times==m)] <-  (exp(-m/2*(wk%*%surv.2)))
    }
    # }
    surva <- surv.probs.1 * surv.probs.2
    return(list(surva = surva, surv.probs.1 = surv.probs.1, surv.probs.2 =surv.probs.2))  
  }
  
  # Calculate CIF 
  cif.function <- function(T.start = T.start, # prediction starting time point, t in paper; it's a point, not vector!
                           tt = tt, # prediction end time points, s in paper
                           u = u, 
                           temp.theta = temp.theta)
  {
    B <- nrow(u)
    Time.n <- length(tt)
    surv.probs.t1.1 <- surv.probs.t1.2 <- vector(length=B)
    surv.probs.m.1 <- surv.probs.m.2 <- matrix(nrow=B, ncol=Time.n)
    x1.new <- new.data$x1
    z1.new <- new.data$z1
    sk <- new.data$sk 
    t1 <- T.start # prediction start point
    t2 <- t1 + tt # prediction end point
    cif.t1.1 <- cif.t1.2 <- cif.m.1 <- cif.m.2 <- cif.t1 <- cif.t2 <- matrix(nrow = B, ncol = Time.n)
    cond.cif.t1 <- cond.cif.t2 <- cond.cif.t1.corrected <- cond.cif.t2.corrected <- matrix(nrow = B, ncol = Time.n)
    colnames(cif.t1) <- round(tt,2)
    colnames(cif.t2) <- round(tt,2)
    ## survival function
    ## t1 part
    for (j in seq(B)){
      surv.t1.1 <-surv.t1.2 <- matrix(NA, nrow = B, ncol=K)
      log.h0.t1s.1 <- log.h0.t1s.2 <- matrix(NA, nrow = B, ncol=K)
      mu.surv.t1 <- matrix(NA, nrow = B, ncol=K)
      lin.pred.t1.1 <- lin.pred.t1.2 <- matrix(NA, nrow = B, ncol=K)
      for (k in 1:K)
      {## Hazard portion
        log.h0.t1s.1[j,k] <- Bspline.tool(t1, sk)$Lam.s[k, ]%*%temp.theta[[j]]$Bs.gammas.t.1
        log.h0.t1s.2[j,k] <- Bspline.tool(t1, sk)$Lam.s[k, ]%*%temp.theta[[j]]$Bs.gammas.t.2
        mu.surv.t1[j,k] <-  x1.new*temp.theta[[j]]$beta.l + z1.new*temp.theta[[j]]$phi.l +
          temp.theta[[j]]$gamma*(t1/2*(sk[k]+1))+ u[j,1]+ u[j,2]*(t1/2*(sk[k]+1))
        
        lin.pred.t1.1[j,k] <- exp(x1.new*temp.theta[[j]]$beta.t.1+z1.new*temp.theta[[j]]$phi.t.1 +
                                    temp.theta[[j]]$eta.t.1*mu.surv.t1[j,k]+
                                    temp.theta[[j]]$zeta.1*u[j,3])
        lin.pred.t1.2[j,k] <- exp(x1.new*temp.theta[[j]]$beta.t.2+z1.new*temp.theta[[j]]$phi.t.2 +
                                    temp.theta[[j]]$eta.t.2*mu.surv.t1[j,k]+
                                    temp.theta[[j]]$zeta.2*u[j,3])
        surv.t1.1[j,k] <- exp(log.h0.t1s.1[j,k])*lin.pred.t1.1[j,k]*surv.function((t1/2*(sk[k]+1)), u[j,], temp.theta = temp.theta[[j]])$surva
        surv.t1.2[j,k] <- exp(log.h0.t1s.2[j,k])*lin.pred.t1.2[j,k]*surv.function((t1/2*(sk[k]+1)), u[j,], temp.theta = temp.theta[[j]])$surva
      }
      
      surv.probs.t1.1[j] <-  t1/2*(wk%*%surv.t1.1[j,]) 
      surv.probs.t1.2[j] <-  t1/2*(wk%*%surv.t1.2[j,])
    }  
    
    ## time m part   
    for (m in t2){
      for (j in seq(B)){
        
        log.h0.tms.1 <- log.h0.tms.2 <- matrix(NA, nrow = B, ncol=K)
        surv.m.1 <-surv.m.2 <- matrix(NA, nrow = B, ncol=K)
        mu.surv.m <- matrix(NA, nrow = B, ncol=K)
        lin.pred.m.1 <- lin.pred.m.2 <- matrix(NA, nrow = B, ncol=K)
        
        for (k in 1:K)
        {## Hazard portion
          
          log.h0.tms.1[j,k] <- Bspline.tool(m, sk)$Lam.s[k, ]%*%temp.theta[[j]]$Bs.gammas.t.1
          log.h0.tms.2[j,k] <- Bspline.tool(m, sk)$Lam.s[k, ]%*%temp.theta[[j]]$Bs.gammas.t.2
          
          mu.surv.m[j,k] <-  x1.new*temp.theta[[j]]$beta.l + z1.new*temp.theta[[j]]$phi.l +
            temp.theta[[j]]$gamma*(m/2*(sk[k]+1))+ u[j,1]+ u[j,2]*(m/2*(sk[k]+1))
          
          lin.pred.m.1[j,k] <- exp(x1.new*temp.theta[[j]]$beta.t.1+z1.new*temp.theta[[j]]$phi.t.1 +
                                     temp.theta[[j]]$eta.t.1*mu.surv.m[j,k]+
                                     temp.theta[[j]]$zeta.1*u[j,3])
          lin.pred.m.2[j,k] <- exp(x1.new*temp.theta[[j]]$beta.t.2+z1.new*temp.theta[[j]]$phi.t.2 +
                                     temp.theta[[j]]$eta.t.2*mu.surv.m[j,k]+
                                     temp.theta[[j]]$zeta.2*u[j,3])
          
          surv.m.1[j,k] <- exp(log.h0.tms.1[j,k])*lin.pred.m.1[j,k]*surv.function((m/2*(sk[k]+1)), u[j,], temp.theta = temp.theta[[j]])$surva
          surv.m.2[j,k] <- exp(log.h0.tms.2[j,k])*lin.pred.m.2[j,k]*surv.function((m/2*(sk[k]+1)), u[j,], temp.theta = temp.theta[[j]])$surva
        }
        
        surv.probs.m.1[j, which(t2==m)] <-  m/2*(wk%*%surv.m.1[j,])
        surv.probs.m.2[j, which(t2==m)] <-  m/2*(wk%*%surv.m.2[j,])
        
        ## calculate the CIF/S
        cif.t1[j, which(t2 ==m)] <- surv.probs.m.1[j, which(t2 ==m)] - surv.probs.t1.1[j]
        cif.t2[j, which(t2 ==m)] <- surv.probs.m.2[j, which(t2 ==m)] - surv.probs.t1.2[j]
        
        
        cond.cif.t1[j, which(t2 ==m)] <- cif.t1[j, which(t2 ==m)]/surv.function(t1, u[j,], temp.theta[[j]])$surva
        cond.cif.t2[j, which(t2 ==m)] <- cif.t2[j, which(t2 ==m)]/surv.function(t1, u[j,], temp.theta[[j]])$surva
        colnames(cond.cif.t1) <- colnames(cond.cif.t2) <- round(tt, 2)
        
      }
    } 
    
    return(list(cond.cif.t1 = cond.cif.t1, cond.cif.t2 = cond.cif.t2))
  }    
  
  cif.out <- c()
  for (i in 1:length(T.start)){
    cif.out[[i]] <- cif.function(T.start = T.start[i], tt = tt, u = u.new, temp.theta = temp.theta)
  }
  
  return(list(cif.out = cif.out)) 
}
