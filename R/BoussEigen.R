#' An eigenvector implementation of the one-dimensional Dupuit-Bousinesq grounwater equation
#'
#'This function takes input data for a catchment about vadose zone recharge, groundwater pumping, soil, the vadose zone and aquifer properties
#'and simulates observed streamflow and well levels.
#'Equations are from:
#' Bidwell, V., Burbery, L., 2011. Groundwater Data Analysis - quantifying aquifer dynamics. Prepared for Envirolink Project 420-NRLC50 (No. 4110/1). Lincoln Ventures Ltd.
#' which in turn cites the following, though note that the symbol labels are different.
#' Bidwell, V.J., Stenger, R., Barkle, G.F., 2008. Dynamic analysis of groundwater discharge and partial-area contribution to Pukemanga Stream, New Zealand. Hydrol. Earth Syst. Sci. 12, 975-987. doi:10.5194/hess-12-975-2008
#' Sloan, W.T., 2000. A physics-based function for modeling transient groundwater discharge at the watershed scale. Water Resour. Res. 36, 225-241. doi:10.1029/1999WR900221
#' @param WellDistance Distance the well is from the upper edge of the groundwater zone in metres
#' @param Storativity the groundwater storativity, i.e. the fraction of space in the groundwater that is available for water
#' @param Transmisivity the two dimensional flow rate of water through the groundwater m2 per day
#' @param ZoneLengths the lengths (in metres) of each of the zones
#' @param DischargeScaleFactor discharge to groundwater level response gain factor
#' @param RechargeData The daily timeseries of vadose recharge for each zone in metres. Output of Vadose.R
#' @param GWBypassFlow The amount of discharge that "slips by" the flow recorder, effectively an offset. It could also be considered to be related to the initial depth to the GW, so it can be a negative number. Ideally calibrate for this using the mean error.
#' @param InitialGWLevel This is a pseudo offset for the groundwater level, but is not quite......
#' @param RiverH This is the constant loss from the groundwater level, in metres, to rivers. Effectively a groundwater level offset
#' @keywords groundwater, hydrology
#' @export

bouss.eigen <- function(WellDistance=65340,Storativity=0.0075,Transmisivity=71940000,
                        ZoneLengths=c(35600,5700,8900,15800),DischargeScaleFactor=500,
                        RechargeData=ZoneVadoseRecharge,GWBypassFlow=1,InitialGWLevel=0,RiverH=0)

{
#Load specific libraries
library(hydroTSM) #this is used to generate the flow duration curve from the discharge timeseries
library(TTR)      #this provides the MovingAverage functions (in particular SMA, Simple Moving Average and EMA exponential weighted moving average)
library(zoo)

NumberEigenvalues          <-  66                                                                #The greater the number the greater the convergence to a true Boussinesq estimate. Cost is time and memeory

NumberOfZones               <-  length(ZoneLengths)
Zone_attribute_names        <-  c("length","recharge_factors","river_recharge","start_length_fraction","end_length_fraction")  #the length fractions refers to the length, as a fraction of the the total groundwater length, where the individual zones start and finish
zone_attributes             <-  array(0, dim=c(NumberOfZones,length(Zone_attribute_names)),
                                      dimnames=list(NULL,Zone_attribute_names))                    #This is simply an array to hold all the zone attributes
zone_attributes[,"length"]  <-  ZoneLengths

# Two sets of eigenvalues are calculated. One for the surface recharge into the groundwater which may be broken up into zones with each having their own eigenvalues,
# and one for the groundwater itself which is just one zone.

zone_eigenfunction_parameters     <-  c("c_ij")                                   #c_ij is taken from Bidwell 2008 equation A6, except it is calculated for each zone
zone_eigenvalues                  <-  array(0, dim=c(NumberOfZones,NumberEigenvalues,length(zone_eigenfunction_parameters)),
                                      dimnames=list(NULL,NULL,zone_eigenfunction_parameters))
general_eigenfunction_parameters  <-  c("k_i","P_i","F_i")                                          #k_i is from Bidwell 2011 equation 8
GeneralEigenvalues               <-  array(0, dim=c(length(general_eigenfunction_parameters),NumberEigenvalues),
                                      dimnames=list(general_eigenfunction_parameters,NULL))

#derived parameters
InitialEigenState                 <- InitialGWLevel - RiverH
groundwater_length                <-  sum(zone_attributes[,"length"])                                   #this is the total length of the groundwater zones
T_over_SL2                        <-  Transmisivity /(Storativity * groundwater_length^2)               #this is used in the model
zone_attributes[,"start_length_fraction"]<-(cumsum(zone_attributes[,"length"])-
                                         zone_attributes[,"length"])/groundwater_length           #this is the distance (as a fraction of the total groundwater reservoir) to the start of each zone
zone_attributes[,"end_length_fraction"]<-cumsum(zone_attributes[,"length"])/groundwater_length         #this is the distance (as a fraction of the total groundwater reservoir) to the end of each zone
GeneralEigenvalues["k_i",]       <-  ((2*col(GeneralEigenvalues)[1,]-1)*pi/2)^2*T_over_SL2

GeneralEigenvalues["P_i",]        <-  cos((2*col(GeneralEigenvalues)[1,]-1)*
                                            pi*WellDistance/groundwater_length/2)/Storativity

GeneralEigenvalues["F_i",]        <- 2*(-1)^(col(GeneralEigenvalues)[1,]+1)/((2*col(GeneralEigenvalues)[1,]-1)*pi)  #note that this Bidwell 2011 F_i is different to Bidwell 2008 Fi

eigenindex                       <-  slice.index(zone_eigenvalues,2)[,,1]                               #This is simply the index of the eigen values
end_length_fraction_index        <-  zone_attributes[slice.index(zone_eigenvalues,1)[,,1],"end_length_fraction"]
dim(end_length_fraction_index)   <-  dim(eigenindex)
start_length_fraction_index      <-  zone_attributes[slice.index(zone_eigenvalues,1)[,,1],"start_length_fraction"]
dim(start_length_fraction_index) <-  dim(eigenindex)

zone_eigenvalues[,,"c_ij"]      <-  (4/(pi*(2*eigenindex-1)))*(sin((2*eigenindex-1)*
                                      pi*end_length_fraction_index/2)-sin((2*eigenindex-1)*pi*start_length_fraction_index/2))

#create the zone timeseries array and populate as much as possible
ZoneTimeseriesNames       <-  c("vadose_recharge","total_recharge")
ZoneTimeseries             <-  array(0, dim=c(nrow(RechargeData),length(ZoneTimeseriesNames),NumberOfZones),
                                      dimnames=list(rownames(RechargeData),ZoneTimeseriesNames,NULL))
ZoneTimeseries[,"vadose_recharge",]<-data.matrix(RechargeData[,1:(NumberOfZones)])
ZoneTimeseries[,"total_recharge",]<-(ZoneTimeseries[,"vadose_recharge",])

#create an empty array ready for catchment timeseries
CatchmentTimeseriesNames  <-  c("estimated_groundwater_depth")
CatchmentTimeseries        <-  array(0, dim=c(nrow(RechargeData),length(CatchmentTimeseriesNames)),
                                      dimnames=list(rownames(RechargeData),CatchmentTimeseriesNames))

#create an empty array, ready for the daily eigenvalue components of the predicted well depth
EigenTimeseries            <-  array(0, dim=c(nrow(RechargeData),NumberEigenvalues,3),dimnames=list(rownames(RechargeData),NULL,c("groundwater_storage","groundwater_level","groundwater_discharge")))

#This function multiplies a dataframe of zone timeseries by a zone parameter for each eigen component and combines the zone results
rechargeByParameter <- function(ZoneSeriesDataframe = ZoneTimeseries[1,"total_recharge",],ZoneParameters=zone_eigenvalues[,,"c_ij"]){
   NumberZones <- length(ZoneSeriesDataframe)
   NumberEigens <- length(ZoneParameters[1,])
   zone_recharge_by_parameter             <- matrix(t(ZoneSeriesDataframe),nrow=NumberZones,ncol=NumberEigens)*ZoneParameters
   rechargeByParameter <- colSums(zone_recharge_by_parameter)
   return(rechargeByParameter)
}

# apply the above function to each row of the recharge data
#by turning each row into an element in a list
totalRechargeList <- split(ZoneTimeseries[,"total_recharge",],row(ZoneTimeseries[,"total_recharge",]))

#calculate gains for each zone by multiplying the zone eigen values by the Pi parmaeter and divide by the Cij parameter
Zone_EigenGains <- zone_eigenvalues[,,1] * t(replicate(3,GeneralEigenvalues["P_i",])) / t(replicate(3,GeneralEigenvalues["k_i",]))

#Now combine the recharge time series with the gains
RechargeByGain <- t(sapply(totalRechargeList,rechargeByParameter,ZoneParameters=Zone_EigenGains))

# Set the first value according to the input argument (this is just helps to reduce the spin-up time)
RechargeByGain[1,1] <- InitialEigenState

#I need, for each eigen component a list that has a list of the daily totalRechargeByCijByEigen values and the Eigenvalue (i.e. a list of lists)
EigenList <- list()
 for (EigenNo in 1:(NumberEigenvalues)){
      EigenList[[length(EigenList)+1]]<-list(series=RechargeByGain[,EigenNo],EigenValue=GeneralEigenvalues["k_i",EigenNo])
 }

#Now I need to create the weighted average of the current and previous eigenvalue
EigenSeries <- sapply(seq_along(EigenList), function(EigenListIndex) {
  #browser()
  GainedRechargeSeries <- EigenList[[EigenListIndex]]$series
  EigenValue <- EigenList[[EigenListIndex]]$EigenValue

  CurrentDaysWeightedValue <- GainedRechargeSeries*(1-exp(-1*EigenValue))
  if(EigenListIndex == 1){CurrentDaysWeightedValue[1] <- InitialEigenState}

  WeightedSeries <- rep(0,length(GainedRechargeSeries))
  for (ValueIndex in 1:length(GainedRechargeSeries)) {
    if(ValueIndex == 1) WeightedSeries[ValueIndex] <- CurrentDaysWeightedValue[ValueIndex]
    else {WeightedSeries[ValueIndex] <- WeightedSeries[ValueIndex-1]*exp(-1*EigenValue)+CurrentDaysWeightedValue[ValueIndex]}
  }
  return(WeightedSeries)
})

EigenTimeseries[,,"groundwater_level"] <- EigenSeries

#Sum the eigenvalues to give the estimated groundwater depth
CatchmentTimeseries[,"estimated_groundwater_depth"]   <-  apply(EigenTimeseries[,,"groundwater_level"],MARGIN=1,sum) + RiverH

#Convert to a zoo data type
CatchmentTimeseries <- as.zoo(CatchmentTimeseries,as.Date(row.names(CatchmentTimeseries),format="%d/%m/%Y"))
} #end of function


