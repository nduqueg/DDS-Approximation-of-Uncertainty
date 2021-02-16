library(ggplot2)
library(zoo)
library(parallel)
library(foreach)
library(doSNOW)
library(tcltk)

dds_au<-function(xBounds.df, numIter, numBeh, per.m.dds, threshold, r = 0.2, OBJFUN, numCores=NA, packages= NA, ...){
  # INPUTS:
  # xBounds.df must be a dataframe with 1st column as minimum and 2nd column as maximum of the parameter range.
  # numIter is an integer which defines the total of model simulations.
  # numBeh is an integer which defines the number of independent DDS and the maximum number of behavioural threshold.
  # per.m.dds is a double which indicates the percentage that affects the computational budget of each independent DDS.
  # threshold is a double which classifies the parameter sets in non-behavioural or behavioural.
  # r is a double between 0 and 1, the default value is 0.2 (DDS parameter of perturbation).
  # OBJFUN is a function which returns a scalar value, for which we are trying to minimize.
  #
  # OUTPUTS: A list that cointains two more list. The first one contains the behavioural sets of parameters and their respective
  # objective function value, and the second one cointains all the simulations classified by the independent DDS number
  # 
  
  # Format xBounds.df and BoundsIter.df colnames
  colnames(xBounds.df)<-c("min", "max")
  
  # Define the number of models evaluations per DDS and ramdomize them around a general number
  m.dds<-round(numIter/numBeh)
  numSearch<-runif(numBeh, m.dds - per.m.dds*m.dds , m.dds + per.m.dds*m.dds)
  
  # Smart Round the values for the m.dds randomly generated
  numSearch <- round(numSearch / sum(numSearch) * numIter)
  deviation <- numIter - sum(numSearch)
  for (. in seq_len(abs(deviation))) {  # adds one to some entries to complete the values
    numSearch[i] <- numSearch[i <- sample(numBeh, 1)] + sign(deviation)
  } 
  
  
  x_init<-apply(xBounds.df, 1, function(x) runif(numBeh, x[1], x[2]))  # initialize the initial first guess
  
  if (is.na(numCores)){  # define parallel processing
    numCores <- 1
  } else if(numCores > detectCores()){
    numCores <- detectCores(); numCores <- numCores -2
  }
  cl <- makeSOCKcluster(numCores)
  registerDoSNOW(cl)
  
  # creation local temporal simulation directory
  if (is.na(dir)){
    dir.create( paste( getwd(), "/local_DDS", sep=""))
    dir <- paste( getwd(), "/local_DDS", sep="")
  }else{
    if( !dir.exists(dir)){ dir.create(dir)}
  }
  #  Parallel processing of the numBeh DDS
  pb <- tkProgressBar(max= numBeh, title="Percentage of independent DDS excecuted")
  progress <- function(n) setTkProgressBar(pb, n, label=paste(n/numBeh*100,"%"), width= 500)
  opts <- list(progress=progress)
  output.list <- foreach (i= 1:numBeh, .export = ls(globalenv()), 
                          .packages= packages,
                          .options.snow=opts) %dopar% {

    all.sim.DDS <- dds(xBounds.df, numSearch[i], iniPar= as.vector(x_init[i,]), r = r, OBJFUN, ...)  # Perform the independent DDS
    
    # Identification of the behavioral sets of parameters
    if ( tail(all.sim.DDS[[2]],1) <= threshold){
      Param <- tail(all.sim.DDS[[1]], 1)
      Obj.Fun <- tail(all.sim.DDS[[2]], 1)
      
    }else{
      Param <- matrix(NA, ncol=ncol(all.sim.DDS[[1]]), nrow=1)
      Obj.Fun <- NA
      
    }
    # local saving 
    save(all.sim.DDS, list="all.sim.DDS", file= paste(dir, "/DDS_",i, ".RData",sep=""))
    
    DDS.output <- list(Param= Param, Obj.Fun= Obj.Fun, All.sim= all.sim.DDS)
  }
  stopCluster(cl)
  close(pb)
  
  all.sim<-list() # initialize the variable where all results will be
  Behavioural<-list()  # initialize the variable wehere only the behavioural will be
  Behavioural[["Param"]]<-as.data.frame(matrix(NA, ncol=nrow(xBounds.df), nrow=numBeh))
  Behavioural[["Obj.Fun."]]<-vector()
  for (i in 1:numBeh){  # extraction of the behavioral so as to have them isolated
    Behavioural[["Param"]][i,] <- output.list[[i]][["Param"]]
    Behavioural[["Obj.Fun."]][i] <- output.list[[i]][["Obj.Fun"]]
    all.sim[[i]] <- output.list[[i]][["All.sim"]]
  }
  
  
  output<-list(Behavioural= Behavioural, all.sim= all.sim)
  return(output)
}

reclas_threshold<-function(all.sim, numBeh, threshold){
  # INPUTS:
  # all.sim is a list containing the second result from the output of the function dds_au.R
  # numBeh is the number of independent DDS used in the funtion dds_au.R
  # threshold is a double defining the new threshold to reclassified the later simulations
  #
  # OUTPUTS:
  # A list that contains the new behavioural sets of parameters
  #
  
  Behavioural<-list()  # initialize the variable wehere only the behavioural will be
  
  for ( i in 1:numBeh){  # reclasiffication
    if (tail(all.sim[[i]][[2]],1) <= threshold){
      Behavioural[["Param"]][i,]<-tail(all.sim[[i]][[1]],1)
      Behavioural[["Obj.Fun."]][i]<-tail(all.sim[[i]][[2]],1)
    }else{
      Behavioural[["Param"]][i,]<-matrix(NA, ncol=ncol(all.sim[[i]][[1]]), nrow=1)
      Behavioural[["Obj.Fun."]][i]<-NA
    }
  }
  
  return(Behavioural)
}

dds<<-function(xBounds.df, numIter,iniPar=NA, r = 0.2, OBJFUN, ...){
  # INPUTS:
  # xBounds.df must be a dataframe with 1st column as minimum, 2nd column as maximum
  # numIter is an integer
  # iniPar is a vector, containing an optional initial parameters
  # r is a double between 0 and 1, the default value is 0.2
  # OBJFUN is a function which returns a scalar value, for which we are trying to minimize.  
  #
  # OUTPUTS:
  # outputs.df is a four entry list, containing x_best, y_best, x_test and y_test, as they evolve over numIter iterations.
  
  # Format xBounds.df colnames
  colnames(xBounds.df)<-c("min", "max")
  # Generate initial first guess
  #xBounds.df<-data.frame(col1 = rep(10,10), col2=rep(100, 10))
  if (is.na(iniPar[1])){  # identification of initial parameters
    x_init<-apply(xBounds.df, 1, function(x) runif(1, x[1], x[2]))
  }else{
    x_init<-as.numeric(iniPar)
  }
  
  x_best<-data.frame(x=x_init)
  x_test<-data.frame(x=x_init)
  
  # Evaluate first cost function
  y_init<-OBJFUN(x_init, ...)
  y_test<-y_init
  y_best<-y_init
  
  # Select which entry to peturb at each iteration
  xDims<-nrow(xBounds.df)
  Prob<-matrix(1-log(1:numIter)/log(numIter),ncol=1) # Returns numIter length list of entries to be peturbed
  peturbIdx<-apply(t(apply(Prob,1, function(x) as.logical(rbinom(xDims, 1, x)))), 1, which)
  # identify where it is not changing any parameter and assign one ramdomly
  Correct.Peturb<-which( unlist( lapply(peturbIdx,sum) ) ==0)
  peturbIdx[Correct.Peturb]<-sample(1:xDims, length(Correct.Peturb), replace=T)
  
  # Peturb each entry by N(0,1)*r(x_max - x_min) reflecting if @ boundaries
  sigma<-xBounds.df$max - xBounds.df$min
  
  pb1<-txtProgressBar(style=3)
  for (i in 2:numIter){
    setTxtProgressBar(pb1,i/numIter,title=paste(i/numIter,"% Calibrando"))
    # Set up test x
    x_test[,i]<-as.matrix(x_best)
    
    # Get entries we will peturb
    idx<-peturbIdx[[i]]
    
    # Initialize vector of peturbations initially zeros with same length of x so we will add this vector to peturb x
    peturbVec<-rep(0, nrow(x_test[,i]))
    # Generate the required number of random normal variables
    N<-rnorm(nrow(x_test[,i]), mean=0, sd=1)
    
    # Set up vector of peturbations
    peturbVec[idx]<-r*N[idx]*sigma[idx]
    
    # Temporary resulting x value if we peturbed it
    x_test[,i]<-x_test[,i] + peturbVec  
    
    # Find the values in testPeturb that have boundary violations.
    B.Vio.min.Idx<-which(x_test[,i] < xBounds.df$min)
    B.Vio.max.Idx<-which(x_test[,i] > xBounds.df$max)
    
    # Correct them by mirroring set them to the minimum or maximum values
    x_test[B.Vio.min.Idx,i]<-xBounds.df$min[B.Vio.min.Idx] + ( xBounds.df$min[B.Vio.min.Idx] - x_test[B.Vio.min.Idx,i] )
    set.min<- B.Vio.min.Idx[ x_test[B.Vio.min.Idx,i] > xBounds.df$max[B.Vio.min.Idx] ]  # which are still out of bound
    x_test[ set.min, i]<- xBounds.df$min[set.min]
    
    x_test[B.Vio.max.Idx,i]<-xBounds.df$max[B.Vio.max.Idx] - ( x_test[B.Vio.max.Idx,i] - xBounds.df$max[B.Vio.max.Idx] )
    set.max<- B.Vio.max.Idx[ x_test[B.Vio.max.Idx,i] < xBounds.df$min[B.Vio.max.Idx] ]
    x_test[ set.max, i]<- xBounds.df$max[set.max]
    
    # Evaluate objective function
    y_test[i]<-OBJFUN(x_test[,i], ...)
    
    y_best<-min(c(y_test[i], y_best))
    bestIdx<-which.min(c(y_test[i], y_best))
    x_choices<-cbind(x_test[,i], as.matrix(x_best))
    x_best<-x_choices[,bestIdx]
  }
  close(pb1)
  
  output.list<-list(X_BEST=t(x_best), Y_BEST=y_best, X_TEST=t(x_test), Y_TEST=y_test)
  return(output.list)
}

post_dds_au<-function(Behavioural, all.sim, obs, dates, MODEL.FUN, ...){
  # INPUTS:
  # 
  # Behavioural is a list containing the first result from the output of the function dds_au.R
  # all.sim is a list containing the second result from the output of the function dds_au.R
  # obs is a vector containing the streamflow observations
  # dates is a dataframe containing the dates of begining and ending of the observations and simulations, and the time interval
  # MODEL.FUN is a function which returns a vector with the simulated streamflow values.  
  # 
  
  ntimes<-length(obs)
  Qsim<-matrix(NA, ncol= nrow(Behavioural[["Param"]]), nrow= ntimes)
  
  # create the simulated timeseries with the behavioural sets of parameters
  for (i in 1:nrow(Behavioural[["Param"]])){
    
    if(sum(is.na(Behavioural[["Param"]][i,]))==0){
      Qsim[,i]<-MODEL.FUN(Behavioural[["Param"]][i,], ...)
    }
  }
  
  
  # calculation of prediction Bounds
  Bounds<-as.data.frame(matrix(NA, nrow=ntimes, ncol=3))
  for (j in 1:ntimes){
    Bounds[j,1]<-min(Qsim[j,],na.rm=T)
    Bounds[j,2]<-mean(Qsim[j,],na.rm=T)
    Bounds[j,3]<-max(Qsim[j,],na.rm=T)
  }
  Bounds.Stats<-list()  # Bounds statistics
  Bounds.Stats[["Mean.Width"]]<-mean(Bounds[,3]-Bounds[,1],na.rm=T)
  Bounds.Stats[["Max.Width"]]<-max(Bounds[,3]-Bounds[,1],na.rm=T)
  Bounds.Stats[["Min.Width"]]<-min(Bounds[,3]-Bounds[,1],na.rm=T)
  Bounds.Stats[["Countaning.R"]]<-sum( obs<=Bounds[,3] & obs>=Bounds[,1], na.rm=T )/ntimes*100
  Bounds.Stats[["Avg.Dev.Amp"]]<-sum( abs( Bounds[,2]-obs ), na.rm=T )/ntimes
  
  # data and bound graphs
  dates.v<-seq(as.Date(dates[1,1]), as.Date(dates[2,1]), by=dates[3,1])
  data<-as.data.frame(cbind(dates.v, Bounds,obs));colnames(data)<-c("Date","Min","Mean","Max","Obs")
  
  
  p<-ggplot(data=data)
  print(p+
          geom_ribbon(aes(x=Date,ymin=Min,ymax=Max,fill="Unc. Bound"),alpha=0.4)+
          scale_fill_manual(values="blue")+
          geom_line(aes(x=Date,y=Obs,col="Qobs"))+scale_color_manual(values="blue")+
          ylim(0,max(obs,na.rm=T))+theme(legend.title = element_blank())+
          ggtitle("Simulations uncertainty bounds and measurements",
                                                        sub="DDS-AU and observed data")+
    theme_bw()
  )
  
  output<-list(data,Bounds.Stats)
  return(output)
}
