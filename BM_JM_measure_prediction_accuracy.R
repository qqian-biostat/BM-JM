#####################################################################################################
#####################################################################################################
## This R file contains functions and codes to calculate dynamic AUC and BS from BM-JM predictions ##
#####################################################################################################
#####################################################################################################

######################################################################
## Define related functions (from Dr. Paul Blanche's timeROC package) 
## Related paper: Blanche P, Dartigues J, Jacqmin-Gadda H (2013). “Estimating and Comparing time-dependent areas 
## under receiver operating characteristic curves for censored event times with competing risks.” Statistics in Medicine, 32(30), 
## 5381–5397. http://onlinelibrary.wiley.com/doi/10.1002/sim.5958/full.
## More details can be found in:
## https://cran.r-project.org/web/packages/timeROC/index.html
######################################################################

##########################################################################################
# Function to compute iid decomposition of KM estimator of the censoring survival distribution
##########################################################################################
## {{ Input
## times: a vector of observed (censored) event times
## status: a vector of event indicators, 1 if non censored, 0 if censored
## }}

## {{ Output
## iid.mat: a matrix of the iid representation of KM for all times in vector times
## }}
Compute.iid.KM <- function(times, status){
  times <- times[order(times)]
  status <- status[order(times)] 
  n <- length(times)
  mat.data<-cbind(times,as.numeric(status==0))
  colnames(mat.data)<-c("T","indic.Cens")
  # Compute the empirical survival function corresponding to the counting process 1(\tilde{eta}=0, \tilde{T}<=t)
  hatSdeltaCensTc<-1-cumsum(mat.data[,c("indic.Cens")])/n  
  # Build the matrix required for computing  dM_C(u) for all time u (all observed times \tilde{T}_i)
  temp1 <- cbind(mat.data[,c("T","indic.Cens")],1-(1:n)/n,hatSdeltaCensTc)
  temp1 <- rbind(c(0,0,1,1),temp1) # Add the first row corresponding to time t=0
  colnames(temp1)<-c("T","indic.Cens","hatSTc","hatSdeltaCensTc")
  # compute hazard function of the censoring
  lambdaC<-(temp1[-1,"indic.Cens"])/(n:1)  
  # Add the column of the hazard function of the censoring (equal to 0 at time t=0)
  temp1<-cbind(temp1,c(0,lambdaC))
  colnames(temp1)[ncol(temp1)]<-"lambdaC"
  # Cumulative hazard of censoring
  LambdaC<-cumsum(lambdaC)         
  # Add the column of the cumulative hazard function of the censoring (equal to 0 at time t=0)
  temp1 <- cbind(temp1,c(0,LambdaC))
  colnames(temp1)[ncol(temp1)]<-"LambdaC"
  temp2<-temp1[-1,]
  # Compute  martingale of censoring \hat{M}_{C_i}(u) for all time u (all observed times \tilde{T}_i) using previous matrix
  # We obtain a matrix. Each column contains the vector of M_{C_i}(\tilde{T}_j) for  all j.
  hatMC<-matrix(NA,n,n)
  for (i in 1:n){
    hatMC[,i] <-temp2[i,2]*as.numeric(temp2[i,1]<=temp2[,"T"])- c(temp2[0:i,"LambdaC"], rep(temp2[i,6],(n-i)))
  }  
  # Compute d \hat{M}_{C_i} (u) for all time u (all observed times \tilde{T}_i)
  dhatMC<-rbind(hatMC[1,],hatMC[-1,]-hatMC[-nrow(hatMC),])
  # Compute d \hat{M}_{C_i} (u)/(S_{\tilde{T}}(u)) for all time u (all observed times \tilde{T}_i)
  # We need this for integrals in the martingale representation of the Kaplan-Meier estimator of the censoring survival function
  # function to divide d \hat{M}_{C_i} (u) by (S_{\tilde{T}}(u))
  MulhatSTc<-function(v){
    n <- length(v)
    v/c(1,1-(1:(n-1))/n)      # c(1,1-(1:(n-1))/n) is the at risk probability (S_{\tilde{T}}(u))
  }
  # Apply the function for each column (corresponding to the
  # vector M_{C_i}(u)  for all time u (all observed times \tilde{T}_i), 
  # time \tilde{T}_i corresponds to the i-th row of the matrix)
  dhatMCdivST<-apply(dhatMC,2,MulhatSTc)
  # Compute \int_0^{\tilde{T}_j} d{ \hat{M}_{C_l} (u) } / (S_{\tilde{T}}(u)) for each subject l, we compute for all time \tilde{T}_j.
  # l=column, j=row
  MatInt0TcidhatMCksurEff<-apply(dhatMCdivST,2,cumsum)  # (Remark : on of the row corresponds to the previous step...) 
  colnames(MatInt0TcidhatMCksurEff)<-paste("M_{C_",1:length(times),"}",sep="")
  rownames(MatInt0TcidhatMCksurEff)<-times  
  return(MatInt0TcidhatMCksurEff)  
}

#####################################################
## Function to calculate dynamic AUC values
#####################################################
## {{ Input
## T: a vector of (censored) event-times
## delta: a vector of event indicators, 1 if non censored, 0 if censored
## marker: a vector of the prediction values for which we want to compute the time-dependen ROC curves.
## other_markers: a matrix that contains values of other predictions that we want to take into account for computing the inverse probability of censoring weights. 
##                Default value is other_markers=NULL.
## cause: a scalar that represents the event of interest for which we aim to compute the time-dependent ROC curve, 1 or 2 if competing risk exists. Without competing risks,
##        it must be the value that indicates a non-censored observation (usually 1).
## weighting: the method used to compute the weights. weighting="marginal" uses the Kaplan-Meier estimator of the censoring distribution. weighting="cox" and
##            weighting="aalen" model the censoring by the Cox model and the additive Aalen model respectively. 
##            Default value is weighting="marginal".
## times: a vector of prediction time points "t" at which we want to compute the time-dependent ROC curve. If vector times contains only a single value, then value zero is added.
## ROC: a logical value that indicates if we want to save the estimates of sensitivities and specificties. Default value is ROC = TRUE.
## iid: a logical value that indicates if we want to compute the iid-representation of the area under time-dependent ROC curve estimator. 
##      iid = TRUE is required for computation of all inference procedures 
## }}

## {{ Output
## iid.mat: a matrix of the iid representation of KM for all times in vector times
## }}

timeROC <- function(T,delta,marker,other_markers=NULL,cause,weighting="marginal",times,ROC=TRUE,iid=FALSE){
  # {{{ check some inputs
  if (length(delta)!=length(T) | length(marker)!=length(T) | length(delta)!=length(T)){
    stop("lengths of vector T, delta and marker have to be equal\n") }
  if (missing(times)){
    stop("Choose at least one time for computing the time-dependent AUC\n") } 
  if (!weighting %in% c("marginal","cox","aalen")){
    stop("the weighting argument must be marginal (default), cox or aalen.\n") }  
  if (weighting %in% c("cox","aalen") & !missing(other_markers) & !("matrix" %in% class(other_markers))){
    stop("argument other_markers must be a matrix\n") }
  if (weighting %in% c("cox","aalen") & !missing(other_markers)){
    if(!nrow(other_markers)==length(marker))  stop("lengths of vector T, delta, marker and number of rows of other_markers have to be equal\n")
  }
  # }}}
  # {{{ check if there are missing values, and delete rows with missing values
  if (weighting %in% c("cox","aalen") & !missing(other_markers) ){
    is_not_na<-as.logical(apply(!is.na(cbind(T,delta,marker,other_markers)),1,prod))
    T<-T[is_not_na]
    delta<-delta[is_not_na]
    marker<-marker[is_not_na]
    other_markers<-as.matrix(other_markers[is_not_na,])
  }else{
    is_not_na<-as.logical(apply(!is.na(cbind(T,delta,marker)),1,prod)) 
    T<-T[is_not_na]
    delta<-delta[is_not_na]
    marker<-marker[is_not_na]
  }
  # }}} 
  start_computation_time<-Sys.time()
  # {{{ create some usefull objects
  n<-length(T)
  n_marker<-length(unique(marker))
  n_times<-length(times)
  if (n_times==1){times<-c(0,times)
  n_times<-2}           # trick to use ipcw.cox() even if there is only one time
  times<-times[order(times)]
  times_names<-paste("t=",times,sep="")
  # }}}
  # {{{ output initialisation
  AUC_1<-rep(NA,n_times)
  AUC_2<-rep(NA,n_times)
  CumInci<-rep(NA,n_times)
  surv<-rep(NA,n_times)
  names(AUC_1)<-times_names
  names(AUC_2)<-times_names
  names(CumInci)<-times_names
  names(surv)<-times_names
  Stats<-matrix(NA,nrow=n_times,ncol=4)
  colnames(Stats)<-c("Cases","survivor at t","Other events at t","Censored at t")
  rownames(Stats)<-times_names
  # }}}
  # {{{  computation of weights (1/2)
  # we need to order to use the pec::ipcw() fonction
  order_T<-order(T)
  T <- T[order_T]
  delta <- delta[order_T]
  marker<- marker[order_T]
  # use ipcw function from pec package
  if(weighting=="marginal"){
    weights <- pec::ipcw(Surv(failure_time,status)~1,data=data.frame(failure_time=T,status=as.numeric(delta!=0)),method="marginal",times=times,subjectTimes=T,subjectTimesLag=1)
  }
  if(weighting=="cox"){
    if (missing(other_markers)){marker_censoring<-marker } 
    other_markers<-other_markers[order_T,]
    marker_censoring<-cbind(marker,other_markers)
    colnames(marker_censoring)<-paste("X", 1:ncol(marker_censoring), sep="")
    fmla <- as.formula(paste("Surv(T,status)~", paste(paste("X", 1:ncol(marker_censoring), sep=""), collapse= "+")))
    data_weight<-as.data.frame(cbind(data.frame(T=T,status=as.numeric(delta!=0)),marker_censoring))
    weights <- pec::ipcw(fmla,data=data_weight,method="cox",times=as.matrix(times),subjectTimes=data_weight[,"T"],subjectTimesLag=1)
  }
  if(weighting=="aalen"){
    if (missing(other_markers)){marker_censoring<-marker }
    other_markers<-other_markers[order_T,]
    marker_censoring<-cbind(marker,other_markers)
    colnames(marker_censoring)<-paste("X", 1:ncol(marker_censoring), sep="")
    fmla <- as.formula(paste("Surv(T,status)~", paste(paste("X", 1:ncol(marker_censoring), sep=""), collapse= "+")))
    data_weight<-as.data.frame(cbind(data.frame(T=T,status=as.numeric(delta!=0)),marker_censoring))
    weights <- pec::ipcw(fmla,data=data_weight,method="aalen",times=as.matrix(times),subjectTimes=data_weight[,"T"],subjectTimesLag=1)
  }
  # we order by marker values (in order to compute Se and Sp)
  order_marker<-order(-marker)
  Mat_data<-cbind(T,delta,marker)[order_marker,]
  colnames(Mat_data)<-c("T","delta","marker")
  # Create some weights
  Weights_cases_all<-1/(weights$IPCW.subjectTimes*n)
  Weights_cases_all<-Weights_cases_all[order_marker]
  # }}}
  # {{{ Make TP and FP outputs if needed
  if(ROC==TRUE){ 
    FP_1<-matrix(NA,nrow=(n_marker+1),ncol=n_times)
    TP<-matrix(NA,nrow=(n_marker+1),ncol=n_times)
    FP_2<-matrix(NA,nrow=(n_marker+1),ncol=n_times)
    colnames(FP_1)<-times_names
    colnames(TP)<-times_names
    colnames(FP_2)<-times_names
  } else{FP_1<-NA
  FP_2<-NA
  TP<-NA}
  # }}}
  # {{{ loop on all timepoints t
  for(t in 1:n_times){
    Cases<-(Mat_data[,"T"]< times[t] &  Mat_data[,"delta"]==cause)
    Controls_1<-(Mat_data[,"T"]> times[t] )
    Controls_2<-(Mat_data[,"T"]< times[t] &  Mat_data[,"delta"]!=cause & Mat_data[,"delta"]!=0)  
    if (weights$method!="marginal"){ 
      Weights_controls_1<-1/(weights$IPCW.times[,t]*n)  }
    else{
      Weights_controls_1<-rep(1/(weights$IPCW.times[t]*n),times=n)
    }
    Weights_controls_1<-Weights_controls_1[order_marker] 
    Weights_cases<-Weights_cases_all 
    Weights_controls_2<-Weights_cases_all
    Weights_cases[!Cases]<-0
    Weights_controls_1[!Controls_1]<-0
    Weights_controls_2[!Controls_2]<-0
    den_TP_t<-sum(Weights_cases)
    den_FP_1_t<-sum(Weights_controls_1)
    den_FP_2_t<-sum(Weights_controls_2)+sum(Weights_controls_1)
    if(den_TP_t!=0){  
      TP_tbis<-c(0,cumsum(Weights_cases))/den_TP_t
      TP_t<-TP_tbis[!duplicated(marker[order_marker])]
    }
    else TP_t<-NA
    if(den_FP_1_t!=0){
      FP_1_tbis<-c(0,cumsum(Weights_controls_1))/den_FP_1_t
      FP_1_t<-FP_1_tbis[!duplicated(marker[order_marker])]}
    else FP_1_t<-NA
    if(den_FP_2_t!=0){
      FP_2_tbis<-c(0,cumsum(Weights_controls_1)+cumsum(Weights_controls_2))/den_FP_2_t
      FP_2_t<-FP_2_tbis[!duplicated(marker[order_marker])]}
    else FP_2_t<-NA
    # internal fonction to compute an area under a curve by trapezoidal rule
    AireTrap<-function(Abs,Ord){
      nobs<-length(Abs)
      dAbs<-Abs[-1]-Abs[-nobs]
      mil<-(Ord[-nobs]+Ord[-1])/2
      area<-sum(dAbs*mil)
      return(area)
    }
    if ( den_TP_t*den_FP_1_t != 0){AUC_1[t]<-AireTrap(FP_1_t,TP_t)}
    else AUC_1[t]<-NA
    if ( den_TP_t*den_FP_2_t != 0){AUC_2[t]<-AireTrap(FP_2_t,TP_t)}
    else AUC_2[t]<-NA
    if(ROC==TRUE){ 
      TP[,t]<-TP_t
      FP_1[,t]<-FP_1_t
      FP_2[,t]<-FP_2_t
    }  
    CumInci[t]<-c(den_TP_t)
    surv[t]<-c(den_FP_1_t)
    Stats[t,]<-c(sum(Cases),sum(Controls_1),sum(Controls_2),n-sum(Cases)-sum(Controls_1)-sum(Controls_2))
  }
  # }}}
  inference<-NA
  if (iid==TRUE){   
    if(weighting!="marginal"){
      stop("Error : Weighting must be marginal for computing the iid representation \n Choose iid=FALSE or weighting=marginal in the input arguments")
    } 
    else{      
      # create iid representation required for inference procedures
      out_iid<-vector("list", n_times)
      names(out_iid)<-paste("t=",times,sep="")
      vect_iid_comp_time<-rep(NA,times=n_times)
      names(vect_iid_comp_time)<-paste("t=",times,sep="")
      mat_iid_rep<-matrix(NA,nrow=n,ncol=n_times)
      colnames(mat_iid_rep)<-paste("t=",times,sep="")  
      mat_iid_rep_star<-matrix(NA,nrow=n,ncol=n_times)
      colnames(mat_iid_rep_star)<-paste("t=",times,sep="")  
      vetc_se<-rep(NA,times=n_times)
      names(vetc_se)<-paste("t=",times,sep="")
      vetc_sestar<-rep(NA,times=n_times)    
      names(vetc_sestar)<-paste("t=",times,sep="")   
      # compute iid for Kaplan Meier
      MatInt0TcidhatMCksurEff <- Compute.iid.KM(times=T,status=delta)
      for (j in 1:n_times){
        #compute iid representation when AUC can be computed
        if(!is.na(AUC_1[j]) | !is.na(AUC_2[j])){
          out_iid[[j]]<-compute_iid_decomposition(t=times[j],n=n,cause=cause,F01t=CumInci[j],St=surv[j],weights,T,delta,marker,MatInt0TcidhatMCksurEff=MatInt0TcidhatMCksurEff)} 
        else{
          out_iid[[j]]<-NA}
        #browser()
        #save output for inference for AUC_1 when AUC_1 can be computed        
        if(!is.na(AUC_1[j])){
          mat_iid_rep_star[,j]<-out_iid[[j]]$iid_representation_AUCstar
          vetc_sestar[j]<-out_iid[[j]]$seAUCstar
          vect_iid_comp_time[j]<-out_iid[[j]]$computation_times               
        }
        #save output for inference for AUC_2 when AUC_2 can be computed         
        if(!is.na(AUC_2[j])){
          mat_iid_rep[,j]<-out_iid[[j]]$iid_representation_AUC
          vetc_se[j]<-out_iid[[j]]$seAUC
          vect_iid_comp_time[j]<-out_iid[[j]]$computation_times               
        }   
      }
      inference<-list(mat_iid_rep_2=mat_iid_rep,
                      mat_iid_rep_1=mat_iid_rep_star,
                      vect_sd_1=vetc_sestar,
                      vect_sd_2=vetc_se,
                      vect_iid_comp_time=vect_iid_comp_time
      )
    }
  }
  stop_computation_time<-Sys.time() 
  # output if there is competing risks or not
  if (max(Stats[,3])==0){
    out <- list(TP=TP,FP=FP_1,AUC=AUC_1,times=times,
                CumulativeIncidence=CumInci,survProb=surv,n=n,Stats=Stats[,c(1,2,4)],weights=weights,
                inference=inference,computation_time=difftime(stop_computation_time,start_computation_time,units="secs"),iid=iid)
    class(out) <- "ipcwsurvivalROC"
    out
  }else{
    out <- list(TP=TP,FP_1=FP_1,AUC_1=AUC_1,FP_2=FP_2,AUC_2=AUC_2,times=times,
                CumulativeIncidence=CumInci,survProb=surv,n=n,Stats=Stats,weights=weights,
                inference=inference,computation_time=difftime(stop_computation_time,start_computation_time,units="secs"),iid=iid)
    class(out) <- "ipcwcompetingrisksROC"
    out    
  }
}

##############################################
# Function to calculate dynamic BS values 
##############################################

## {{ Inpus:
## timepoints, a vector of time points for which we aim to compute the iid representations
## times, a vector of observed (censored) time-to-event
## status, a vector of event indicators with 1=uncensored, 0=censored
## pred, a matrix dynamic prediction results, e.g. pi(t,s), 
##       with each row for each subject in the testing data set, and
##       each column for each prediction time.
## cause, a scalar (1 or 2), cause for which we aim to compute the expected Brier score estimator
## compute.iid, a logical value that indicates if we want to need to compute iid decomposition of KM estimator.
##      iid = TRUE is required for computation of all inference procedures 
# }}

# {{ Outputs :
# a matrix with iid representation for all prediction time points
# }}

BS <- function(timepoints, times, status, pred, cause=1, compute.iid=TRUE){ 
  start_computation_time <- Sys.time()
  # define useful objects
  n <- length(times)
  n_times <- length(timepoints)
  timepoints <- timepoints[order(timepoints)]
  times_names <- paste("t=",timepoints,sep="")
  # output initialisation 
  BS <- rep(NA,n_times)
  CumInci <- rep(NA,n_times)
  surv <- rep(NA,n_times)
  Stats <- matrix(NA,nrow=n_times,ncol=4)
  hit1_all <- matrix(NA,nrow=n,ncol=n_times)
  hit2_all <- matrix(NA,nrow=n,ncol=n_times)
  epsilon_i <- matrix(NA,nrow=n,ncol=n_times)
  #adds name to outputs
  names(BS) <- times_names
  names(CumInci) <- times_names
  names(surv) <- times_names
  colnames(Stats) <- c("Cases","survivor at t","Other events at t","Censored at t")
  rownames(Stats) <- times_names
  colnames(epsilon_i) <- times_names
  colnames(hit1_all) <-  times_names
  colnames(hit2_all)  <- times_names 
  # we need to order to use the ipcw() function of the `pec` package
  order_T <- order(times)
  times <-  times[order_T]
  delta  <-  status[order_T]
  pred <-  pred[order_T,,drop=FALSE]
  #compute KM weights
  weights <- ipcw(Surv(failure_time,status)~1,
                  data=data.frame(failure_time=times,status=as.numeric(delta!=0)),
                  method="marginal",times=timepoints,subjectTimes=times,subjectTimesLag=1)
  Mat_data <- cbind(times,delta,as.numeric(delta==0))
  colnames(Mat_data) <- c("T","delta","indic_Cens")
  # computate weights of cases
  Weights_cases_all <- 1/(weights$IPCW.subjectTimes*n)
  # compute KM censoring estimator iid representation
  if (compute.iid){ MatInt0TcidhatMCksurEff <- Compute.iid.KM(times,delta!=0)}
  # loop on all time points
  for(t in 1:n_times){
    Cases <- (Mat_data[,"T"]<= timepoints[t] &  Mat_data[,"delta"]==cause)
    Controls_1 <- (Mat_data[,"T"]> timepoints[t] )
    Controls_2 <- (Mat_data[,"T"]<= timepoints[t] &  Mat_data[,"delta"]!=cause & Mat_data[,"delta"]!=0)  
    # compute weights
    Weights_controls_1 <- rep(1/(weights$IPCW.times[t]*n),times=n)
    Weights_cases <- Weights_cases_all
    Weights_controls_2 <- Weights_cases_all
    Weights_cases[!Cases] <- 0
    Weights_controls_1[!Controls_1] <- 0
    Weights_controls_2[!Controls_2] <- 0   
    # compute outputs
    CumInci[t] <- c(sum(Weights_cases))
    surv[t] <- c(sum(Weights_controls_1))
    Stats[t,] <- c(sum(Cases),sum(Controls_1),sum(Controls_2),n-sum(Cases)-sum(Controls_1)-sum(Controls_2)) 
    hit1_all[,t] <- (Weights_controls_1*((pred[,t])^2))*n
    hit2_all[,t] <- (Weights_cases*((1-pred[,t])^2) + Weights_controls_2*((pred[,t])^2))*n
    BS[t] <- (sum(hit1_all[,t]) +sum(hit2_all[,t]))/n
    if (compute.iid){
      # compute 
      Int0tdMCsurEffARisk <- MatInt0TcidhatMCksurEff[max(which(Mat_data[,"T"]<=timepoints[t])),]   
      #browser()
      epsilon_i[,t] <- hit1_all[,t]+hit2_all[,t]-BS[t] + mean(hit1_all[,t])*Int0tdMCsurEffARisk +  colMeans(MatInt0TcidhatMCksurEff*hit2_all[,t])
    }
  } 
  # compute mean and sd of iid representation
  sd_all <- rep(NA,n_times)
  mean_all <- rep(NA,n_times)
  if (compute.iid){sd_all <- apply(epsilon_i,2,sd)/sqrt(n)
  mean_all <- apply(epsilon_i,2,mean)}
  # compute a table to print 
  print.tab <- cbind(Stats,BS,sd_all,mean_all)
  colnames(print.tab) <- c(colnames(Stats),"BS","sd","mean_iid")
  # compute the computation time
  stop_computation_time <- Sys.time()
  computation_time=difftime(stop_computation_time,start_computation_time,units="secs")
  
  out <- list(BS=BS,iid=epsilon_i,sd=sd_all,res=(hit1_all+hit2_all),
              CumulativeIncidence=CumInci,survProb=surv,n=n,Stats=Stats,print.tab=print.tab,timepoints=timepoints,
              computation_time=difftime(stop_computation_time,start_computation_time,units="secs"))
  class(out) <- "ipcwEBS"
  out 
}

##############################################################
# Set the prediction time t and prediction time window \Delta t
##############################################################
T.start = c(0, 0.2)
tt = c(0.2, 0.4, 0.6, 0.8, 1)

#######################################
## Collect dynamic prediction outputs 
#######################################
# For tt = 0.2
TT0_1_02 <- TT0_2_02 <- TT2_1_02 <- TT2_2_02 <- rep(NA, n.sample)
for (i in 1:n.sample){
    TT0_1_02[i] <- mean(PredOut[[i]]$cif.out[[1]]$cond.cif.t1[,1], na.rm = T)
    TT0_2_02[i] <- mean(PredOut[[i]]$cif.out[[1]]$cond.cif.t2[,1], na.rm = T)
    TT2_1_02[i] <- mean(PredOut[[i]]$cif.out[[2]]$cond.cif.t1[,1], na.rm = T)
    TT2_2_02[i] <- mean(PredOut[[i]]$cif.out[[2]]$cond.cif.t2[,1], na.rm = T)
 }
TT02 <- data.frame(cbind(TT0_1_02, TT0_2_02, TT2_1_02, TT2_2_02))
TT02_2 <- TT02
TT02_2[which(TT02 < 0, arr.ind=TRUE)] <- 0
TT02_2[which(TT02 > 1, arr.ind=TRUE)] <- 1
  
  
# For tt = 0.4
TT0_1_04 <- TT0_2_04 <- TT2_1_04 <- TT2_2_04 <- rep(NA, n.sample)
for (i in 1:n.sample){
    TT0_1_04[i] <- mean(PredOut[[i]]$cif.out[[1]]$cond.cif.t1[,2], na.rm = T)
    TT0_2_04[i] <- mean(PredOut[[i]]$cif.out[[1]]$cond.cif.t2[,2], na.rm = T)
    TT2_1_04[i] <- mean(PredOut[[i]]$cif.out[[2]]$cond.cif.t1[,2], na.rm = T)
    TT2_2_04[i] <- mean(PredOut[[i]]$cif.out[[2]]$cond.cif.t2[,2], na.rm = T)
 }
TT04 <- data.frame(cbind(TT0_1_04, TT0_2_04, TT2_1_04, TT2_2_04))
TT04_2 <- TT04
TT04_2[which(TT04 < 0, arr.ind=TRUE)] <- 0
TT04_2[which(TT04 > 1, arr.ind=TRUE)] <- 1
  
  
# For tt = 0.6
TT0_1_06 <- TT0_2_06 <- TT2_1_06 <- TT2_2_06 <- rep(NA, n.sample)
for (i in 1:n.sample){
    TT0_1_06[i] <- mean(PredOut[[i]]$cif.out[[1]]$cond.cif.t1[,3], na.rm = T)
    TT0_2_06[i] <- mean(PredOut[[i]]$cif.out[[1]]$cond.cif.t2[,3], na.rm = T)
    TT2_1_06[i] <- mean(PredOut[[i]]$cif.out[[2]]$cond.cif.t1[,3], na.rm = T)
    TT2_2_06[i] <- mean(PredOut[[i]]$cif.out[[2]]$cond.cif.t2[,3], na.rm = T)
  }
TT06 <- data.frame(cbind(TT0_1_06, TT0_2_06, TT2_1_06, TT2_2_06))
TT06_2 <- TT06
TT06_2[which(TT06 < 0, arr.ind=TRUE)] <- 0
TT06_2[which(TT06 > 1, arr.ind=TRUE)] <- 1
  
# For tt = 0.8
TT0_1_08 <- TT0_2_08 <- TT2_1_08 <- TT2_2_08 <- rep(NA, n.sample)
for (i in 1:n.sample){
    TT0_1_08[i] <- mean(PredOut[[i]]$cif.out[[1]]$cond.cif.t1[,4], na.rm = T)
    TT0_2_08[i] <- mean(PredOut[[i]]$cif.out[[1]]$cond.cif.t2[,4], na.rm = T)
    TT2_1_08[i] <- mean(PredOut[[i]]$cif.out[[2]]$cond.cif.t1[,4], na.rm = T)
    TT2_2_08[i] <- mean(PredOut[[i]]$cif.out[[2]]$cond.cif.t2[,4], na.rm = T)
  }
TT08 <- data.frame(cbind(TT0_1_08, TT0_2_08, TT2_1_08, TT2_2_08))
TT08_2 <- TT08
TT08_2[which(TT08 < 0, arr.ind=TRUE)] <- 0
TT08_2[which(TT08 > 1, arr.ind=TRUE)] <- 1
  
# For tt = 1.0
TT0_1_010 <- TT0_2_010 <- TT2_1_010 <- TT2_2_010 <- rep(NA, n.sample)
for (i in 1:n.sample){
    TT0_1_010[i] <- mean(PredOut[[i]]$cif.out[[1]]$cond.cif.t1[,5], na.rm = T)
    TT0_2_010[i] <- mean(PredOut[[i]]$cif.out[[1]]$cond.cif.t2[,5], na.rm = T)
    TT2_1_010[i] <- mean(PredOut[[i]]$cif.out[[2]]$cond.cif.t1[,5], na.rm = T)
    TT2_2_010[i] <- mean(PredOut[[i]]$cif.out[[2]]$cond.cif.t2[,5], na.rm = T)
  }
TT010 <- data.frame(cbind(TT0_1_010, TT0_2_010, TT2_1_010, TT2_2_010))
TT010_2 <- TT010
TT010_2[which(TT010 < 0, arr.ind=TRUE)] <- 0
TT010_2[which(TT010 > 1, arr.ind=TRUE)] <- 1
  
# Summarize prediction output and the true survival times and event status
# Pass them to the functions defined above
pred.times <- T.start   # vector of prediction times (t)
windowt <- tt  # vector of prediction time windows (\Delta t)


## Matrix to store dynamic AUC and BS values from the 1st competing risk event
AUCst.1 <- matrix(NA,ncol = length(pred.times), nrow = length(windowt))
rownames(AUCst.1) <- windowt
colnames(AUCst.1) <- pred.times
BrierS.s.1 <- matrix(NA,ncol = length(pred.times), nrow = length(windowt))
rownames(BrierS.s.1) <- windowt
colnames(BrierS.s.1) <- pred.times
## Matrix to store dynamic AUC and BS values from the 2nd competing risk event
AUCst.2 <- matrix(NA,ncol = length(pred.times), nrow = length(windowt))
rownames(AUCst.2) <- windowt
colnames(AUCst.2) <- pred.times
BrierS.s.2 <- matrix(NA,ncol = length(pred.times), nrow = length(windowt))
rownames(BrierS.s.2) <- windowt
colnames(BrierS.s.2) <- pred.times


TT <- list()
TT[[1]] <- cbind(pred.subject$dat.terminal[,2:3], TT02_2)
TT[[2]] <- cbind(pred.subject$dat.terminal[,2:3], TT04_2)
TT[[3]] <- cbind(pred.subject$dat.terminal[,2:3], TT06_2)
TT[[4]] <- cbind(pred.subject$dat.terminal[,2:3], TT08_2)
TT[[5]] <- cbind(pred.subject$dat.terminal[,2:3], TT010_2)

######################################################################
# Calculate dynamic AUC and BS values using the functions defined above
######################################################################
for (k in 1:length(windowt)){
    
    rownames(TT[[k]]) <- 1:50
    d <- TT[[k]]
    colnames(d)[1] <- "time"
    colnames(d)[2] <- "status"

    for (s in pred.times){
      print(paste("computation for prediction time s=", s))
      # Create prediction data set
      d.s <- d[,c("time","status")]
      d.s$Pred.s.1 <- d[,paste("TT",(10*s),"_1_0",(10*windowt[k]),sep="")]
      d.s$Pred.s.2 <- d[,paste("TT",(10*s),"_2_0",(10*windowt[k]),sep="")]
      d.s <- d.s[d.s$time>s,]
      d.s$time.s <- d.s$time-s
      # AUC and BS for predicting the 1st competing risk event
      # Estimate ROC curve and AUC
      ROC.s.1 <- timeROC(T=d.s$time.s,
                           delta=d.s$status,
                           marker=d.s$Pred.s.1,
                           cause=1,weighting="marginal",
                           times=c(windowt[k]),
                           iid=FALSE)
      # Estimate expected Brier score
      BS.s.1 <- BS(timepoints=c(windowt[k]),
                     times=d.s$time.s,
                     status=d.s$status,
                     pred=as.matrix(d.s$Pred.s.1),
                     cause=1)
      # Save useful results
      BrierS.s.1[k,which(s==pred.times)] <- BS.s.1$BS # BS estimate
      if (is.null(ROC.s.1$AUC_2[2])){
        AUCst.1[k,which(s==pred.times)] <- ROC.s.1$AUC[2] # AUC estimate
      } else {AUCst.1[k,which(s==pred.times)] <- ROC.s.1$AUC_2[2]}
      
      # AUC and BS for predicting the 1st competing risk event
      # Estimate ROC curve and AUC
      ROC.s.2 <- timeROC(T=d.s$time.s,
                            delta=d.s$status,
                            marker=d.s$Pred.s.2,
                            cause=2,weighting="marginal",
                            times=c(windowt[k]),
                            iid=FALSE)
      # Estimate expected Brier score
      BS.s.2 <- BS(timepoints=c(windowt[k]),
                      times=d.s$time.s,
                      status=d.s$status,
                      pred=as.matrix(d.s$Pred.s.2),
                      cause=2)
      # Save useful results
      BrierS.s.2[k,which(s==pred.times)] <- BS.s.2$BS # BS estimate
      if (is.null(ROC.s.2$AUC_2[2] )){
        AUCst.2[k,which(s==pred.times)] <- ROC.s.2$AUC[2] # AUC estimate
      } else {AUCst.2[k,which(s==pred.times)] <- ROC.s.2$AUC_2[2]}
    }
  }

########################################################################################
# Get dynamic AUC and BS values at different prediction times and prediction time windows
########################################################################################
AUCst.1
AUCst.2
BrierS.s.1
BrierS.s.2

