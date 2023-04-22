# compile and load necessary c-code for this script
# This requires an Rtools installation

if(file.exists("include.so")) file.remove("include.so")
if(file.exists("include.o")) file.remove("include.o")
system("R CMD SHLIB include.c")

if(.Platform$OS.type == "unix"){
  dyn.load("include.so")
} else if(.Platform$OS.type == "windows"){
  dyn.load("include.dll")
} else{
  warning("please check the name of your compiled shared library in this folder and enter it int include.R")
  # HERE
}

# This script uses some helper functions taken from LaplacesDemon. 
# These are in this separate file
source("LDfunctions.R")

#' Transform a vector to a predefined range.
#' 
#' @param x The vector to be rescaled/transformed.
#' @param f A 2 element vector specifying the range for rescaling on the scale of x
#' @return A vector of numbers scaled between 0 and 1.
#' @details The function is defined as mynormalise2<-function(x,f) (x - f[1] ) / (f[2]-f[1])
#' @examples
#' x<- 1:10
#' mynormalise2(x, c(0,1) )
mynormalise2<-function(x,f) (x - f[1] ) / (f[2]-f[1])

mytruncate <- function(x, a, b) pmin(pmax(x, a), b)    

#' Produce a data input list for running the TTR-SDM-FQR or TTR-SDM-RED, both of which use a Farquhar style photosynthesis model to describe carbon assimilation
#' 
#' @param SD A data frame with the input data structured as the example data set. This function is for runnning the TTR-SDM using a photosynthesis model. 
#' @param C3 Is the plant C3 or C4, default C3=TRUE.
#' @param CA Atmospheric CO2 partial pressure Pa, default CA=40.
#' @return A list containing parameters and data for runnning TTR-SDM-FQR or TTR-SDM-RED.
#' @details 
#' The function prodcdata uses sd to create the input data structure needed for the TTR model. This function creates a list that arranges the environmental forcing data in the following format: each variable (e.g. the mean monthly temperature) is packed into a vector. The first 12 data points represent the 12 monthly values of that variable for the first site. The next 12 data points are the 12 monthly values for the second site, and so on. For soil nitrogen only a single value is available (the data are total soil N) in which case this value is simply repeated 12 times to give a constant monthly value. The values of the environmental variables are additionally rescaled using the \code{mynormalise2} function in \code{prodcdata} to a scale approximately between 0 and 1. The normalisation parameters can be user-defined, but should be aligned within the bounds of the optimisation routine. This transformation should be noted so that a back transformation can be applied. That is, the model parameters will be on this normalised scale if you want them for example back on a degrees Celsius scale, you need to apply a back transformation. 
prodcdata.chelsa<-function(SD,C3=TRUE,CA=40) {
  # these normalisation constants are important later for backtransforming
  fTMEAN<-c(-60,60)
  fTMIN<-c(-60,60)
  fTMAX<-c(-60,60)
  fM<-c(0,100)
  fN<-c(0,100) 
  fPHOTO<-c(0,30)

  n<-dim(SD)[1]
  
  names<-paste("tavg",1:12,sep="_")
  TMEANseq<-mynormalise2( as.vector(t(as.matrix(SD[,names]))),fTMEAN )
  names<-paste("tmin",1:12,sep="_")
  TMINseq<-mynormalise2( as.vector(t(as.matrix(SD[,names]))), fTMIN )
  names<-paste("swc",1:12,sep="_")
  Mseq<-mynormalise2( as.vector(t(as.matrix(SD[,names]))),fM ) 
  Nseq <- Mseq*0 + 1  #assume Nsoil is constant
  names<-paste("tmax",1:12,sep="_")
  TMAXseq <- mynormalise2( as.vector(t(as.matrix(SD[,names]))), fTMAX )
  
  
  names<-paste("tmax",1:12,sep="_")
  TLEAF<-as.vector(t(as.matrix(SD[,names]))) + 273
  names<-paste("ppfd",1:12,sep="_")
  PAR<-as.vector(t(as.matrix(SD[,names]))) 
  
  if(C3==TRUE) {
    PHOTO<-C3loop(tlseq=TLEAF,PARseq=PAR,caseq=rep(CA,n*12),p3)
  }
  if(C3==FALSE) {
    PHOTO<-C4loop(tlseq=TLEAF,PARseq=PAR,caseq=rep(CA,n*12),p4)
  }
  
  Aseq<-mynormalise2( pmax(0,PHOTO), fPHOTO )*0.2 #net photo < 0 simply set to zero!
  
  c.data <- list(SITES=length(SD[,1]),
                 obs=SD$Presence,
                 steps=300,
                 initials=0.1,
                 TMEAN=TMEANseq,
                 TMIN=TMINseq,
                 TMAX=TMAXseq,
                 M=Mseq,
                 N=Nseq,
                 PAR=PAR,
                 TLEAF=TLEAF,
                 PHOTO=PHOTO,
                 A=Aseq,
                 lon=SD$lon,
                 lat=SD$lat)
  
  #parameters in the TTR model that are treated as constants                         
  c.data$Kl<-0.1 
  c.data$gs<-20 
  c.data$gr<-20 
  c.data$KM<-2.5 
  c.data$A0<-0.1
  c.data$N0<-0.02
  c.data$KA<-3 
  c.data$Jc<-0.1
  c.data$Jn<-0.01
  c.data$q<-1
  c.data$RHOc<-1
  c.data$RHOn<-1
  c.data$Fc<-0.5
  c.data$Fn<-0.025
  return(c.data)
}

#' Produce an input data frame for use in prodcdata.chelsa
#' 
#' @param pts coordinates as expected by terra::extract
#' @param CHELSA.files vector of full path names of CHELSA climatology input files
#' @param CA Atmospheric CO2 partial pressure Pa, default CA=40.
#' @return A list containing parameters and data for runnning TTR-SDM-FQR or TTR-SDM-RED.
#' @details 
#' This function is not used in the tutorial, but is provided as an example of how to generate a data frame for prodcdata.chelsa. The user would need to ensure that the appropriate input data files are available. 
make.sd.data <- function( pts,
                          CHELSA.files,
                          FC.file,
                          WP.file)
{
#-1- CHELSA climate layers; this will vary with which CHELSA product one uses!
    CHELSA <- rast(CHELSA.files)
    TMAX <- CHELSA[[grep("tasmax",names(CHELSA))]]
    TMIN <- CHELSA[[grep("tasmin",names(CHELSA))]]
    TAVG <- CHELSA[[grep("tas_",names(CHELSA))]]
    RAIN <- CHELSA[[grep("pr",names(CHELSA))]]
    PENMAN <- CHELSA[[grep("penman",names(CHELSA))]]    
    RSDS   <- CHELSA[[grep("rsds",names(CHELSA))]]    
    VPD    <- CHELSA[[grep("vpd",names(CHELSA))]]    
    rm(CHELSA)
    
#-2- Field capacity & Wilting point
    FC   <- rast(FC.file)
    WP <- rast(WP.file)
                      
#-3- extract to data frame
    Tmax    <- extract(TMAX, pts,ID=FALSE)
    Tmin    <- extract(TMIN, pts,ID=FALSE)
    Tavg    <- extract(TAVG, pts,ID=FALSE)
    Rain    <- extract(RAIN, pts,ID=FALSE)
    rsds    <- extract(RSDS,pts,ID=FALSE)        
    pet     <- extract(PENMAN,pts,ID=FALSE)        
    vpd     <- extract(VPD,pts,ID=FALSE) #not used!       
    wp      <- unlist(extract(WP, pts,ID=FALSE))  
    fc      <- unlist(extract(FC, pts,ID=FALSE))  
    
#-4- unit conversion
    # rsds is MJ/m2/d: 
    # 86400 for per day to per second
    # 1e06 from Mega W to W
    rsdsW <- rsds / 86400 * 1e06 
    J_to_mol=4.6;
    frac_PAR=0.5;
    ppfd <- rsdsW * frac_PAR * J_to_mol # conversion from rsds W/m2 to PAR umol/m2/s

#-5- estimate monthly soil water    
    sw<-make.swc(Rain, pet, fc, wp)

#-6- put everything in a data frame    
    sd<-cbind(as.data.frame(pts[,c(1,2)]),Tmax,Tavg,Tmin,sw,ppfd,vpd,pet,Rain,wp,fc)
    names(sd)[1:2] <- c("lon","lat")
    names(sd)[3:14] <- paste("tmax", 1:12, sep="_")
    names(sd)[15:26] <- paste("tavg", 1:12, sep="_")
    names(sd)[27:38] <- paste("tmin", 1:12, sep="_")
    names(sd)[39:50] <- paste("swc", 1:12, sep="_")
    names(sd)[51:62] <- paste("ppfd", 1:12, sep="_")
    names(sd)[63:74] <- paste("vpd", 1:12, sep="_")
    names(sd)[75:86] <- paste("pet", 1:12, sep="_")
    names(sd)[87:98] <- paste("prec", 1:12, sep="_")
    names(sd)[99] <- "wp"
    names(sd)[100] <- "fc"
    return(sd)
}   

#' Run a simple soil water balance model.
#' 
#' @param RAIN monthly rainfall data from CHELSA
#' @param PET montly potential evapotranspiration data from CHELSA
#' @param FC soil field capacity
#' @param WP soil wilting point
#' @return Monthly soil water content
#' @details 
#' This function is not used in the tutorial, but is called by make.sd.data.
make.swc<-function(RAIN, PET, FC, WP)
{
  soilwater.extract <- matrix(nrow=length(FC), ncol=12) 
  YEARS = 3
  SWC.m<-matrix(0, ncol = 12*YEARS+1, nrow=length(FC))
  SWC.m[,1]<-FC*1000

  for (r in 1:nrow(SWC.m)) { # for each cell
    i<-2
    for(y in 1:YEARS) { # for each year 
      for(m in 1:12) { # for each month
        EP <- RAIN[r,m]  # Effective precipitation mm/month
        Kswc <- trap1(SWC.m[r,i-1],WP[r]*1000,FC[r]*1000) # constrained to 0..1
        AET <- Kswc * PET [r,m]
        SWC.m[r,i] <- SWC.m[r,i-1] + EP - AET 
        SWC.m[r,i] <- min(FC[r]*1000,SWC.m[r,i]) # SWC cannot exceed FC, excess is lost to runoff/drainage
        SWC.m[r,i] <- max(WP[r]*1000,SWC.m[r,i]) # SWC cannot go below WP
        i<-i+1
      } 
    }
  }

  # Take last 12 months for output
  SWC.m <- SWC.m[, seq(ncol(SWC.m)-11,ncol(SWC.m))]
  # Scale between 0 and 100, representing WP and FC
  for(r in 1:nrow(SWC.m)){
    Y = trap1(SWC.m[r,],WP[r]*1000,FC[r]*1000)*100
    # Write results in matrix
    if( is.na(sum(Y)) ) {
      soilwater.extract[r, ] <- -999
    } else {
      soilwater.extract[r, ] <- round(Y)
    }
  }
  soilwater.extract[soilwater.extract==-999] <- NA
  return(soilwater.extract)
}


#' Run the TTR-SDM-RED model for multiple sites forced by environmental data.
#' 
#' @param args A list created by \code{prodcdata.A}.
#' @return Predicted total plant biomass for each site.
#' @details The function is not necessarily directly called by the user. 
#' @examples
#' \dontrun{ ThornLoopAR.C(args) }
'ThornLoopAR.C' <- 
  function(args){
  abu<-rep(0,length=args$SITES)
  res <- .C("ThornLoopAR", # what the function is called in the c code
    as.integer(args$steps),  # number steps to run
    as.double(args$initials),  # initial pool size
    as.integer(args$SITES), # number of sites	
    as.double(args$TMAX),  
    as.double(args$TMEAN),  
    as.double(args$TMIN),  
    as.double(args$M),  
    as.double(args$N),  
    as.double(args$A),  
    as.double(args$Kl),  
    as.double(args$gs),  
    as.double(args$gr),  
    as.double(args$KM),  
    as.double(args$A0),  
    as.double(args$N0),  
    as.double(args$KA),  
    as.double(args$Jc),  
    as.double(args$Jn),  
    as.double(args$q),  
    as.double(args$RHOc),  
    as.double(args$RHOn),  
    as.double(args$Fc),  
    as.double(args$Fn),  
    as.double(args$ma1),   
    as.double(args$ma2),  
    as.double(args$tn1),  
    as.double(args$tn2),  
    as.double(args$mn1),  
    as.double(args$mn2),  
    as.double(args$mn3),  
    as.double(args$mn4),  
    as.double(args$tg1),  
    as.double(args$tg2),  
    as.double(args$tg3),   
    as.double(args$tg4),   
    as.double(args$mg1),   
    as.double(args$mg2),   
    as.double(args$tr1),  
    as.double(args$tr2),  
    as.double(args$na1),   
    as.double(args$na2),   
    ABUSITE = as.double(abu))
  
  return(list(abu=res$ABUSITE))	
}

#' A wrapper function for calling the TTR model.
#' 
#' @param Several parameters needed by the TTR model.
#' @return Predicted total plant biomass for each site.
#' @details The function is not necessarily directly called by the user. 
#'  This function calls ThornSpaceAR.C which is an R wrapper to the c code 
#'  that runs the TTR model. This is long and ugly, but it makes subsequent 
#'  code more compact
callTTR <- function( steps,
                     initials,
                     sites,
                     TMAX,
                     TMEAN,
                     TMIN,
                     M,
                     N,
                     A,
                     Kl,
                     gs,
                     gr,
                     KM,
                     A0,
                     N0,
                     KA,
                     Jc,
                     Jn,
                     q,
                     RHOc,
                     RHOn,
                     Fc,
                     Fn,
                     ma1,   #from here fitted
                     ma2,  
                     tn1,  
                     tn2,  
                     mn1,  
                     mn2,  
                     mn3,  
                     mn4,  
                     tg1,  
                     tg2,  
                     tg3,   
                     tg4,   
                     mg1,   
                     mg2,   
                     tr1,  
                     tr2,
                     na1,
                     na2
)  {
  
  args=list()
  args$steps = steps
  args$initials = initials*0.5
  args$SITES = sites
  args$TMAX = TMAX
  args$TMEAN = TMEAN
  args$TMIN = TMIN
  args$M = M
  args$N = N
  args$A = A
  args$Kl = Kl
  args$gs = gs
  args$gr = gr
  args$KM = KM
  args$A0 = A0
  args$N0 = N0
  args$KA = KA
  args$Jc = Jc
  args$Jn = Jn
  args$q = q
  args$RHOc = RHOc
  args$RHOn = RHOn
  args$Fc = Fc
  args$Fn = Fn
  args$ma1 = ma1    
  args$ma2 = ma2   
  args$tn1 = tn1
  args$tn2 = tn2  
  args$mn1 = mn1  
  args$mn2 = mn2  
  args$mn3 = mn3   
  args$mn4 = mn4  
  args$tg1 = tg1  
  args$tg2 = tg2  
  args$tg3 = tg3   
  args$tg4 = tg4   
  args$mg1 = mg1   
  args$mg2 = mg2   
  args$tr1 = tr1  
  args$tr2 = tr2  
  args$na1 = na1  
  args$na2 = na2         
  abu<-ThornLoopAR.C(args)$abu
  abu[which( is.nan(abu))] <- 0.0
  abu[which( is.na(abu))]  <- 0.0
  abu[which( !is.finite(abu))] <- 0.0
  abu[abu==0]  <- 0.0
  return(abu)
}        





#' Format the parameters as a list
#' 
#' @param x Vector of parameters.
#' @return List of named parameters
#' @details The function is not necessarily directly called by the user. 
as.parm.names.list <- function(x){
  #get the total number of elements of the list
  out <- rep("", sum(unlist(lapply(x,length))))
  ipar <- 1
  for (i in 1:length(x)) {
    #for vectors, add [<index>] after the name
    par <- x[[i]]
    if(length(par) > 1){
      for(j in 1:length(par)){
        out[[ipar]] <- paste(names(x)[[i]], "[", ipar, "]", sep='')
        ipar <- ipar + 1
      }
    }
    #for single values, copy the name
    else {
      out[[ipar]] <- names(x)[[i]]
      ipar <- ipar + 1
    }
  }
  return(out)
}

#' Define a wrapper function for taking data formatted by prodcdata.chelsa and formatting for 
#' DEoptim
#' 
#' @param data.A A list preoduced by prodcdata.chelsa
#' @return A list formatted for DEoptim
make.MyData <- function(data.A) {
  # upper and lower bounds for parameters
  lo <- c(0, 0, 
          0, 0, 
          0, 0, 0, 0, 
          0, 0, 0, 0, 
          0, 0, 
          0, 0,
          0, 0,
          15/1000)  
  up <- c(1,  0.75,  
          1,  0.75,  
          1,  0.75,  1,  1,  
          1,  0.75,  1,  1, 
          1,  0.75,   
          1,  0.75,
          1,  5.0,
          5) 
  
  mu<-up*0
  
  # parameter names for convience 
  names<-c("ma1","ma2","tn1","tn2","mn1","mn2","mn3","mn4",
           "tg1","tg2","tg3","tg4","mg1","mg2",
           "tr1","tr2", "na1","na2",
           "m")
  
  bound.mat <- cbind(lo,up,mu)
  rownames(bound.mat)<-names
  round(bound.mat,4)
  
  # create parameter names
  parm.names <- as.parm.names.list(list(parset=rep(0,19)))
  
  #-put the forcing data in a matrix
  X <- as.matrix(with(data.A,cbind(A,TMAX,TMEAN,TMIN,M,N)))
  #-put the response data (presence/absence) in a vector
  y <- as.numeric(data.A$obs)
  steps <-data.A$steps # number of time steps
  sites <-data.A$SITES # number of sites
  
  #-take the TTR constants out of data.A
  MyData <- data.A
  MyData <- MyData[names(MyData) %in% c("Kl","gs","gr","KM","A0","N0","KA","Jc","Jn","q","RHOc","RHOn","Fc","Fn")]      
  
  #-join all the required elements together
  MyData <- c(MyData, list(steps=steps, sites=sites, X=X, parm.names=parm.names,
                           pos.parset=1:19,y=y,
                           mu=bound.mat[,"mu"],lo=bound.mat[,"lo"],up=bound.mat[,"up"]) )
  return(MyData)
}


#' Return the biomass predicted by the model
#' 
#' @param parm A vector of parameters describing the species physiological prefererences. 
#' @param Data list of input forcing data created by \code{make.MyData}. 
#' @return A vector of predicted biomass for each site
ModelABU <- function(parm, Data)
{
  #--parameters
  parset <- parm[Data$pos.parset]
  for(i in 1:length(Data$pos.parset))
  {
    parset[i] <- mytruncate( parm[i], Data$lo[i], Data$up[i] )
  }
  x0 <- 0.2
  
  #--parameters for TTR
  ma1<-parset[1]    
  ma2<-ma1+parset[2]  
  tn1<-parset[3]     
  tn2<-tn1+parset[4]    
  mn1<-parset[5]   
  mn2<-mn1+parset[6]     
  mn3<-mn2+parset[7]     
  mn4<-mn3+parset[8]     
  tg1<-parset[9]     
  tg2<-tg1+parset[10]    
  tg3<-tg2+parset[11]    
  tg4<-tg3+parset[12]    
  mg1<-parset[13]    
  mg2<-mg1+parset[14]     
  tr1<- parset[15]     
  tr2<- tr1 + parset[16]     
  na1<- parset[17]     
  na2<- na1 + parset[18]     
  #--call TTR, this is a wrapper function to the c code
  #--it returns a modelled time series of biomass values
  x <- callTTR(Data$steps,x0,Data$sites,Data$X[,"TMAX"],Data$X[,"TMEAN"],Data$X[,"TMIN"],
               Data$X[,"M"],Data$X[,"N"],Data$X[,"A"],
               Data$Kl,Data$gs,Data$gr, 
               Data$KM,Data$A0,Data$N0,Data$KA,Data$Jc,Data$Jn,
               Data$q,Data$RHOc,Data$RHOn,Data$Fc,Data$Fn,
               ma1,ma2,tn1,tn2,mn1,mn2,mn3,mn4,
               tg1,tg2,tg3,tg4,mg1,mg2,tr1,tr2,na1,na2)
  x<- x/parset[19]                    
  
  #--return elements 
  p01 <- invcloglog(log(x)) + 1e-22 
  return.list = list( p01=p01,  x=x)
  return(return.list)
} 

#' Estimate the log-likihood for the model 
#' 
#' @param parm A vector of parameters describing the species physiological prefererences. 
#' @param Data list of input data created by \code{make.MyData}. 
#' @return The negative log-likelihood of the parameters given the data.
#' @details This function is designed to be used by \code{DEoptim}.
ModelDE <- function(parm, Data)
{
  #--parameters
  parset <- parm[Data$pos.parset]
  for(i in 1:length(Data$pos.parset))
  {
    parset[i] <- mytruncate( parm[i], Data$lo[i], Data$up[i] )
  }
  x0 <- 0.2
  
  #--parameters for TTR
  ma1<-parset[1]    
  ma2<-ma1+parset[2]  
  tn1<-parset[3]     
  tn2<-tn1+parset[4]    
  mn1<-parset[5]   
  mn2<-mn1+parset[6]     
  mn3<-mn2+parset[7]     
  mn4<-mn3+parset[8]     
  tg1<-parset[9]     
  tg2<-tg1+parset[10]    
  tg3<-tg2+parset[11]    
  tg4<-tg3+parset[12]    
  mg1<-parset[13]    
  mg2<-mg1+parset[14]     
  tr1<- parset[15]     
  tr2<- tr1 + parset[16]   
  na1<- parset[17]     
  na2<- na1 + parset[18]     
  #--call TTR, this is a wrapper function to the c code
  #--it returns a modelled time series of biomass values
  x <- callTTR(Data$steps,x0,Data$sites,Data$X[,"TMAX"],Data$X[,"TMEAN"],Data$X[,"TMIN"],
               Data$X[,"M"],Data$X[,"N"],Data$X[,"A"],
               Data$Kl,Data$gs,Data$gr, 
               Data$KM,Data$A0,Data$N0,Data$KA,Data$Jc,Data$Jn,
               Data$q,Data$RHOc,Data$RHOn,Data$Fc,Data$Fn,
               ma1,ma2,tn1,tn2,mn1,mn2,mn3,mn4,
               tg1,tg2,tg3,tg4,mg1,mg2,tr1,tr2,na1,na2)
  x<- x/parset[19]                    
  
  #--Log-Likelihood  
  p01 <- invcloglog(log(x)) +1e-22
  LL <- sum( dbern(Data$y,p01,log=TRUE) )
  return( -1*LL)
}

#' Return some model evaluation statistics us ROCR
#' 
#' @param bestmem A vector of parameters describing the species physiological prefererences. 
#' @param Data A vector of observed occurrences.
#' @param PLOT Should a AUC plot be produced (TRUE or FALSE).
#' @return A vector of model performance metrics.
#' @details Adatpted from https://www.r-bloggers.com/2014/12/a-small-introduction-to-the-rocr-package/
get_eval <- function(bestmem,Data,PLOT=TRUE) {
# code from 
# https://www.r-bloggers.com/2014/12/a-small-introduction-to-the-rocr-package/
    opt.cut = function(perf, pred){
        cut.ind = mapply(FUN=function(x, y, p){
            d = (x - 0)^2 + (y-1)^2
            ind = which(d == min(d))
            c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
                cutoff = p[[ind]], ind=ind)
        }, perf@x.values, perf@y.values, pred@cutoffs)
    }

    Model_out <- ModelABU(bestmem,Data) 

    pred <- ROCR::prediction(as.numeric(Model_out$p01),Data$y)
    roc.perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
    stats <- opt.cut(roc.perf, pred)
    auc.perf <- ROCR::performance(pred, measure = "auc")

    if(PLOT==TRUE) {
        pdf(file="./roc.pdf",width=5,height=5) 
        plot(roc.perf, main = paste0("AUC=",round(unlist(auc.perf@y.values),3) ) )
        abline(a=0, b= 1,lty=2)
        dev.off()
    }

    ret <- c ( auc = unlist(auc.perf@y.values),
                cutoff = as.numeric(unlist(stats["cutoff",1])),
                sensitivity = as.numeric(unlist(stats["sensitivity",1])),
                specificity = as.numeric(unlist(stats["specificity",1])),
                tp = unlist(pred@tp)[stats["ind",1]],
                tn = unlist(pred@tn)[stats["ind",1]],
                fp = unlist(pred@fp)[stats["ind",1]],
                fn = unlist(pred@fn)[stats["ind",1]])
    return( ret )
}



##---------the C4 photosyntehsis parameters-stored in a list---------------#
p4<-list()
p4$Vcmax25   <- 32  
p4$Vpmax25   <- 191
p4$ps2a25    <- 0.069 
p4$ps2b25    <- 0.701 
p4$ps2c25    <- 773 
p4$sprime    <- 0.319
p4$x         <- 0.404 
p4$alpha     <- 0.001
p4$Vpr       <- 80 
p4$Kc25      <- 125.143  
p4$Ko25      <- 29199
p4$Kp25      <- 43.286   
p4$Sco25     <- 1310  
p4$gbs25     <- 0.038
p4$Rd25      <- 0.555
p4$Rm25      <- 0.5*p4$Rd25     
p4$a025      <- 0.047
p4$EVcmax <- 72.178 
p4$EVpmax <- 46.686 
p4$Eps2c  <- 19.015
p4$Eps2  <- 31.206
p4$Sps2  <- 0.619
p4$Hps2  <- 191.858
p4$EKc  <- 65.026     
p4$EKo  <- 10.268     
p4$EKp  <- 41.666     
p4$ESco <- -29.982
p4$Egbs <- 15.152
p4$ERd  <- 67.113
p4$Ea0  <- 1.63
p4$z  <- 1
p4$m<-4      
p4$b<-0.04   
p4$g0 <- 0.01 
p4$g1 <- 1.62 
p4$bb <- 0 #if using bally berry make this 1 else medlyn



##---------the C3 photosynthesis parameters-stored in a list---------------#
p3<-list()
p3$Vcmax25   <- 73 
p3$ps2a25    <- 0.151 
p3$ps2b25    <- 0.834 
p3$ps2c25    <- 781 
p3$sprime    <- 0.315
p3$Kc25     <-29.829    
p3$Ko25     <-178000 
p3$Gstar25  <-3.970 
p3$gm25     <-3.650   
p3$Rd25     <-1.264
p3$EVcmax   <- 49.251
p3$Eps2c    <- 44.329
p3$Eps2    <- 26.200
p3$Sps2    <- 0.648
p3$Hps2    <- 198.771
p3$EKc      <- 65.051
p3$EKo      <- 34.192
p3$EGstar   <- 24.045
p3$Egm      <- 50.586
p3$ERd      <- 49.393
p3$m<-9      
p3$b<-0.01   
p3$g0 <- 0.01 
p3$g1 <- 4    
p3$bb <- 0 #if using bally berry make this 1 else medlyn




VPDcalc <- function(t, rh)
{
  # VPD needed by Medlyn equation
  # rh as a proportion
  # es and ea in Pa
  # return VPS in kPa
  #Jones. Plants and Microclimate. Second edition. page 110.
  es = 613.75*exp(17.502*t/(t+240.97)); #saturated vapour pressure Pa
  ea = rh*es; #actual vapour pressure Pa
  return ((es-ea)/1000);
}



##-------temperature dependency functions-----------------------##   

pf<-function(KTemp,E,S,H,y25,KTref=298,R=0.008314) { #peaked Arrhenius function
  firstexp<-exp(((KTemp-KTref)*E)/(KTref*R*KTemp))
  topexp<-(1+exp((KTref*S-H)/(KTref*R)))
  bottomexp<-(1+exp((KTemp*S-H)/(KTemp*R))) 
  r<-y25*firstexp*(topexp/bottomexp)
  return(r)
}

npf<-function(KTemp,E,y25,KTref=298,R=0.008314) { #non-peaked Arrhenius function
  y25*exp((KTemp-KTref)*E/(KTref*R*KTemp))
}


##--------C3 model equations--------------------------------------#
acC3 <- function(  cm, om, tl, pp ) {
  # Rubisco-limited photosynthetic rate (umol/(m^2s^1))
  Vcmax<- npf(tl,pp$EVcmax,pp$Vcmax25)
  Ko<-npf(tl,pp$EKo,pp$Ko25) 
  Kc<-npf(tl,pp$EKc,pp$Kc25) 
  Gstar<-npf(tl,pp$EGstar,pp$Gstar25)
  Rd<-npf(tl,pp$ERd,pp$Rd25)
  gm<-npf(tl,pp$Egm,pp$gm25)
  a1 <- (-1/gm)
  b1 <- (Vcmax - Rd)/gm + cm + Kc*(1 + om/Ko)
  c1 <- Rd*(cm + Kc*(1 + om/Ko)) - Vcmax*(cm - Gstar)
  tosqrt <- pmax(0, b1^2 - 4*a1*c1 )
  Ac <- (-b1+sqrt(tosqrt)) / (2*a1)
  return(Ac)
}

ajC3 <- function( cm, tl, PAR, pp) {
  # Electron transport - limited photosynthetic rate (umol/(m^2s^1))
  Gstar<-npf(tl,pp$EGstar,pp$Gstar25)
  Rd<-npf(tl,pp$ERd,pp$Rd25)
  gm<-npf(tl,pp$Egm,pp$gm25)
  #---ps2 method------------------------------------------------------#
  ps2a<- pp$ps2a25 #
  ps2b<- pp$ps2b25 #
  ps2c<- npf(tl,pp$Eps2c,pp$ps2c25)  #
  ps2phi25 <-  ps2a + (ps2b - ps2a) * exp( - PAR / ps2c )
  ps2phi<- pf(tl,pp$Eps2,pp$Sps2,pp$Hps2,ps2phi25)  #
  J <- pp$sprime*PAR*ps2phi  # no x, for C3 ie no partitioning between C3 and C4 cycles
  #---ps2 method------------------------------------------------------#
  a2<-(-1/gm)
  b2<-((J/4) - Rd)/gm + cm + 2*Gstar
  c2<-Rd*(cm+2*Gstar) - (J/4)*(cm - Gstar)
  tosqrt <- pmax(0, b2^2 - 4*a2*c2 )
  Aj <- (-b2+sqrt(tosqrt)) / (2*a2)
  return(Aj)
}


aC3<-function( cm, om, tl, PAR, pp) { # combined C3 function
  ac <- acC3(  cm, om, tl, pp )
  aj <- ajC3(  cm, tl, PAR, pp )
  return( pmin(ac,aj) ) 
}


optC3<-function( cm, om,tl,PAR, p3, gsc, cs) {
  A<-aC3( cm, om, tl, PAR, p3 )
  Ad<- 10*(cs-cm) / (1/gsc)  # *10 coverts from Pa to ppm (umol/mol)
  A<-aC3( cm, om, tl, PAR, p3 )
  d<-abs(A-Ad)
}

#' Main loop to run the C3 photosynthesis model for a set of sites with monthly forcing data.
#'
#' @param tlseq A vector of leaf temperatures in Kelvin
#' @param PARseq A vector of photosynthetic active radiation values in micro mol per m2 per second
#' @param caseq A vector of atmospheric CO2 values in Pa
#' @param p3 A list specifying the parameters of the photosynthesis model
#' @return A vector of photosynthetic values
C3loop<-function(tlseq,PARseq,caseq,p3) {
  oa<-21000
  a<-rep(5,length(tlseq))
  cld = 0.01 #characteristic leaf dimension in m
  u = 1.5 #wind velocity at reference height in m/s
  gbw = 0.223*sqrt(u/cld)*10 # gbw is boundary layer conductance to water Bonan eq 15.3 pg 234 in mol per m2 per second, *10 to get to per Pa
  h = 0.5 #relative humdity 
  nsites<-length(tlseq)/12 #loop to run each site for 24 months
  for(site in 1:nsites) {
    for(month in 1:12) {
      i <- (site-1)*12 + month 
      tl <-tlseq[i] 
      PAR<-PARseq[i]     
      a.temp <- aC3( caseq[i]*0.9, oa, tl, PAR, p3 ) 
      cs = caseq[i] 
      # p3$bb is a switch to either use bb or medlyn equation:
      VPD = VPDcalc(tl-273,h)
      gsc = p3$bb*((p3$m*a.temp*h/(cs*10)+p3$b)/1.6) + (1-p3$bb)*(p3$g0 + 1.6*(1.0 + p3$g1/sqrt(VPD)) * a.temp/(cs*10)) 
      gsc = pmin(1, gsc)
      cm = optimise(optC3,interval=c(0.1,caseq[i]),om=oa,tl=tl,PAR=PAR, p3=p3, gsc=gsc, cs=cs,tol=0.2)$minimum 
      a[i]<-aC3( cm, oa, tl, PAR, p3 )
    } } 
  return(a)    
}


##--------C4 model equations--------------------------------------#    
vpf <- function( cm, tl, pp ){ 
  # CO2 concentrating flux (umol/m2/s)
  # von Caemmerer eq 4.19
  Kp<-npf(tl,pp$EKp, pp$Kp25)
  vpmax <- npf(tl,pp$EVpmax,pp$Vpmax25)  
  return( pmin((cm*vpmax)/(cm + Kp), pp$Vpr) )
}

acC4 <- function(  cm, om, tl, pp ) {
  # Rubisco-limited photosynthetic rate (umol/(m^2s^1))
  # Cm and Om are (mesophyll partial pressures)
  # gbs is bundle sheath conductance to CO2
  Vcmax <- npf(tl,pp$EVcmax,pp$Vcmax25)  # Maximum carboxylation rate (umol/(m^2s))
  Ko    <- npf(tl,pp$EKo,pp$Ko25) #Michaelis-menten coefficient for O2
  Kc    <- npf(tl,pp$EKc,pp$Kc25) #Michaelis-menten coefficient for CO2
  sco   <- npf(tl,pp$ESco,pp$Sco25)
  gammastar<-0.5/sco
  Rd    <- npf(tl,pp$ERd,pp$Rd25)
  Rm    <- Rd*0.5
  gbs   <- npf(tl,pp$Egbs,pp$gbs25)
  a0    <- npf(tl,pp$Ea0,pp$a025)
  Vp    <- vpf(cm, tl, pp) 
  a1 <- 1 - p4$alpha/a0 * Kc/Ko
  b1 <- -1 * ( (Vp - Rm + gbs * cm) + (Vcmax - Rd) +  gbs*(Kc*(1+om/Ko)) + (p4$alpha/a0 * (gammastar*Vcmax + Rd * Kc/Ko ) ) )
  c1 <- (Vcmax - Rd) * (Vp - Rm + gbs*cm) -  (Vcmax*gbs*gammastar*om + Rd*gbs*Kc*(1+om/Ko)) 
  tosqrt <- pmax(0, b1^2 - 4*a1*c1 ) 
  Ac<- (-b1 - sqrt( tosqrt ) ) / (2 * a1)
  return(Ac)
}


ajC4 <- function( cm, om, tl, PAR, pp) {
  #Electron transport - limited photosynthetic rate (umol/(m^2s^1))
  sco <- npf(tl,pp$ESco,pp$Sco25)
  gammastar <- 0.5/sco
  Rd  <- npf(tl,pp$ERd,pp$Rd25)    
  Rm  <- Rd*0.5
  gbs <- npf(tl,pp$Egbs,pp$gbs25)
  a0  <- npf(tl,pp$a025,pp$Ea0)
  z   <- pp$z
  #---ps2 method------------------------------------------------------#
  ps2a<- pp$ps2a25 #
  ps2b<- pp$ps2b25 #
  ps2c<- npf(tl,pp$Eps2c,pp$ps2c25)  #
  ps2phi25 <-  ps2a + (ps2b - ps2a) * exp( - PAR / ps2c ) 
  ps2phi<- pf(tl,pp$Eps2,pp$Sps2,pp$Hps2,ps2phi25)  #
  Jatp <- pp$sprime*PAR*ps2phi / (1 - pp$x ) 
  #---ps2 method------------------------------------------------------#
  a2 <- 1 - 7 * gammastar * pp$alpha / (3 * a0)
  b2 <- -1* ( (pp$x*Jatp * z / 2 - Rm + gbs*cm)
              + ( (1 - pp$x) * Jatp * z / 3 - Rd) 
              + gbs*( 7* gammastar * om  / 3)
              + pp$alpha * gammastar / a0 * ( (1 - pp$x) * Jatp * z / 3 + 7*Rd / 3 )    )
  c2 <- ( ( pp$x*Jatp * z /2 - Rm + gbs*cm ) * ( (1 - pp$x) * Jatp * z / 3 - Rd) ) 
  - gbs*gammastar*om*( (1 - pp$x) * Jatp * z / 3 +  7*Rd / 3 ) 
  tosqrt <- pmax(0, b2^2  - 4*a2*c2) 
  Aj <- (-b2 - sqrt( tosqrt ) ) / (2 * a2)
  return(Aj)
}


aC4<-function( cm, obs, tl, PAR, pp) { # combined C4 function
  ac <- acC4(  cm, obs, tl, pp )  
  aj <- ajC4(  cm, obs, tl, PAR, pp ) 
  return( pmin(ac,aj) ) 
}


optC4<-function( cm, om,tl,PAR, p4, gsc, cs) {
  Ad<- 10*(cs-cm) / (1/gsc) # *10 coverts from Pa to ppm (umol/mol)
  A<-aC4( cm, om, tl, PAR, p4 )
  d<-abs(A-Ad)
}


#' Main loop to run the C4 photosynthesis model for a set of sites with monthly forcing data.
#'
#' @param tlseq A vector of leaf temperatures in Kelvin
#' @param PARseq A vector of photosynthetic active radiation values in micro mol per m2 per second
#' @param caseq A vector of atmospheric CO2 values in Pa
#' @param p4 A list specifying the parameters of the photosynthesis model. If bb=1 then the Bally Berry model is used if bb=0 the Medlyn model is used.
#' @return A vector of photosynthetic values
C4loop<-function(tlseq,PARseq,caseq,p4) {
  oa<-21000
  a<-rep(5,length(tlseq))
  cld = 0.01 #characteristic leaf dimension in m
  u = 1.5 #wind velocity at reference height in m/s
  gbw = 0.223*sqrt(u/cld)*10 # gbw is boundary layer conductance to water Bonan eq 15.3 pg 234 in mol per m2 per second, *10 to get to per Pa
  h = 0.5 #relative humdity
  nsites<-length(tlseq)/12 #loop to run each site for 24 months
  for(site in 1:nsites) {
    for(month in 1:12) {
      i <- (site-1)*12 + month 
      tl <-tlseq[i] 
      PAR<-PARseq[i]     
      a.temp <- aC4( caseq[i]*0.9, oa, tl, PAR, p4 ) 
      cs = caseq[i]
      # p4$bb is a switch to either use bb or medlyn equation:
      VPD = VPDcalc(tl-273,h)
      gsc = p4$bb*((p4$m*a.temp*h/(cs*10)+p4$b)/1.6) + (1-p4$bb)*(p4$g0 + 1.6*(1.0 + p4$g1/sqrt(VPD)) * a.temp/(cs*10)) 
      gsc = pmin(1, gsc)
      cm = optimise(optC4,interval=c(0.1,caseq[i]),om=oa,tl=tl,PAR=PAR, p4=p4, gsc=gsc, cs=cs,tol=0.2)$minimum 
      a[i]<-aC4( cm, oa, tl, PAR, p4 )
    }}
  return(a)        
}



