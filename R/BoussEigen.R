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
#' @param Storativity                #the groundwater storativity, i.e. the fraction of space in the groundwater that is available for water
#' @param Transmisivity              #the two dimensional flow rate of water through the groundwater m2 per day
#' @param ZoneLengths                #the lengths (in metres) of each of the zones
#' @param DischargeScaleFactor       #discharge to groundwater level response gain factor
#' @param RechargeFileName           #The csv file with all the data in it. This is the daily timeseries of observed discharge and groundwater level, vadose recharge for each zone, river recharge to groundwater and groundwater pumping in mm
#' @param GWBypassFlow               #The amount of discharge that "slips by" the flow recorder, effectively an offset. It could also be considered to be related to the initial depth to the GW, so it can be a negative number. Ideally calibrate for this using the mean error.                                       
#' @param InitialEigenState          #This is a pseudo offset for the groundwater level, but is not quite......
#' @keywords groundwater, hydrology
#' @export
#' @examples
#' fishRecharge<-vadose.recharge(ZoneStorageTime=c(4,3,3,3),AquiferZoneDryFractions = list(c(0.8,0.8,0.8,0.8),c(0,0,0,0),c(0,0,0.193,0.8263)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.0069,0.1737)),RiverRechargeFractions=c(1,1,0.958,1),RechargeFileName=system.file("extdata","GoldenBayLandSurfaceRechargeData.csv",package="GDATools"),PumpingFileName=system.file("extdata","GoldenBayGWPumpingData.csv",package="GDATools"))
#' fishDischarge <- bouss.eigen(WellDistance=66000,ZoneLengths=c(35600,5700,8900,15800),Storativity=0.01,Transmisivity=313632,DischargeScaleFactor=950,RechargeData=fishRecharge,GWBypassFlow=1)



bouss.eigen <- function(WellDistance=65340,Storativity=0.0075,Transmisivity=71940000,
                        ZoneLengths=c(35600,5700,8900,15800),DischargeScaleFactor=500,
                        RechargeData=ZoneVadoseRecharge,GWBypassFlow=1,InitialEigenState=0)

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

#create recharge data file column headings
ZoneNames                  <-  paste0("Z",seq(1:NumberOfZones))
RechargeTypeNames         <-  c("dryland","irrigated","river","pumped")

#build the column headings to be read from the input data
data_headings<-c(t(outer(ZoneNames,RechargeTypeNames, paste,sep="_")))

ZoneTimeseriesNames       <-  c("vadose_recharge","total_recharge")
ZoneTimeseries             <-  array(0, dim=c(nrow(RechargeData),length(ZoneTimeseriesNames),NumberOfZones),
                                      dimnames=list(rownames(RechargeData),ZoneTimeseriesNames,NULL))

#ZoneTimeseries[,"river_recharge",]<- data.matrix(RechargeData[,(2*NumberOfZones):(3*NumberOfZones-1)])
#ZoneTimeseries[,"pumped_removal",]<- data.matrix(RechargeData[,(3*NumberOfZones):(4*NumberOfZones-1)])
ZoneTimeseries[,"vadose_recharge",]<-data.matrix(RechargeData[,1:(NumberOfZones)])

#Calculate the total groundwater recharge
ZoneTimeseries[,"total_recharge",]<-(ZoneTimeseries[,"vadose_recharge",])

#create an array for catchment timeseries
CatchmentTimeseriesNames  <-  c("estimated_groundwater_depth","scenario_groundwater_depth","estimated_groundwater_discharge","scaled_est_discharge","scaled_est_discharge_less_SW_takes","scenario_groundwater_discharge")
CatchmentTimeseries        <-  array(0, dim=c(nrow(RechargeData),length(CatchmentTimeseriesNames)),
                                      dimnames=list(rownames(RechargeData),CatchmentTimeseriesNames))

#create the array for the daily eigenvalue components of the predicted well depth
EigenTimeseries            <-  array(0, dim=c(nrow(RechargeData),NumberEigenvalues,3),dimnames=list(rownames(RechargeData),NULL,c("groundwater_storage","groundwater_level","groundwater_discharge")))

#This function multiplies the recharge totals for each zone by their respective Cij parameter for each eigen component
rechargeByCij <- function(totalRecharge = ZoneTimeseries[1,"total_recharge",],Eigenvalues=zone_eigenvalues[,,"c_ij"]){
   NumberZones <- length(totalRecharge)
   NumberEigens <- length(Eigenvalues[1,])
   zone_recharge_by_c_ij             <- matrix(t(totalRecharge),nrow=NumberZones,ncol=NumberEigens)*Eigenvalues
   rechargeByCij <- colSums(zone_recharge_by_c_ij)
   return(rechargeByCij)
}

# apply the above function to each row of the recharge data
#by turning each row into an element in a list
totalRechargeByCijList <- split(ZoneTimeseries[,"total_recharge",],row(ZoneTimeseries[,"total_recharge",]))
#then apply the function to each list
totalRechargeByCij <- t(sapply(totalRechargeByCijList,rechargeByCij,Eigenvalues=zone_eigenvalues[,,"c_ij"]))  #note that I have transformed it back to have each row a different day

# Divide these "Eigen" totals by the Eigenvalue.
totalRechargeByCijByEigen <- sweep(totalRechargeByCij,MARGIN=2, GeneralEigenvalues["k_i",], FUN='/')

# Set the first value according to the input argument (this is just helps to reduce the spin-up time)
totalRechargeByCijByEigen[1,1] <- InitialEigenState

#Which enables the use of the EMA function to calculate the daily eigentotals.
#so I need, for each eigen component a list that has a list of the daily totalRechargeByCijByEigen values and the Eigenvalue (i.e. a list of lists)
EigenList <- list()
for (EigenNo in 1:(NumberEigenvalues)){
  EigenList[[length(EigenList)+1]]<-list(series=totalRechargeByCijByEigen[,EigenNo],EigenValue=GeneralEigenvalues["k_i",EigenNo])
}

#And then actually apply the weighted averaging
EigenTimeseriesEMA <- sapply(EigenList, function(x) EMA(x$series,n=1,ratio=1-exp(-1*x$EigenValue)))
EigenTimeseries[,,"groundwater_storage"] <- EigenTimeseriesEMA

#calculate the eigenvalue components of the predicted well depth
EigenTimeseries[,,"groundwater_level"]            <-  EigenTimeseries[,,"groundwater_storage"] * rep(GeneralEigenvalues["P_i",],each=nrow(EigenTimeseries[,,"groundwater_storage"]))

#I think the following line is equivalent
#EigenTimeseries[,,"groundwater_level"]            <-  EigenTimeseries[,,"groundwater_storage"] *t(replicate(nrow(EigenTimeseries[,,"groundwater_storage"]),GeneralEigenvalues["P_i",]))

#if (i==1) {
#    if (!is.na(CatchmentTimeseries[1,"observed_groundwater_depth"])){
#      EigenTimeseries[i,1,"groundwater_level"]<-  CatchmentTimeseries[1,"observed_groundwater_depth"]}                           #Initial principal eigenvalues set to match observed depth if observed depth is available
#}

#calculate the eigenvalue components of the predicted groundwater discharge
#for (i in 1:2) {             #used in debugging

  EigenTimeseries[,,"groundwater_discharge"]            <-  EigenTimeseries[,,"groundwater_storage"] * rep(GeneralEigenvalues["F_i",]*GeneralEigenvalues["k_i",],each=nrow(EigenTimeseries[,,"groundwater_storage"]))
#  if (i==1) {
#    if (!is.na(CatchmentTimeseries[1,"observed_groundwater_discharge"])){
#      EigenTimeseries[i,1,"groundwater_discharge"]<-  CatchmentTimeseries[1,"observed_groundwater_discharge"]}                           #Initial principal eigenvalues set to match observed flow if observed is available
#}

#Sum the eigenvalues to give the estimated groundwater depth
CatchmentTimeseries[,"estimated_groundwater_depth"]   <-  apply(EigenTimeseries[,,"groundwater_level"],MARGIN=1,sum)
#if (!is.na(CatchmentTimeseries[1,"observed_groundwater_depth"])){
#  CatchmentTimeseries[1,"estimated_groundwater_depth"]    <-  CatchmentTimeseries[1,"observed_groundwater_depth"]    #Initial estimated set to match observed depth
#}

#Sum the eigenvalues to give the estimated groundwater discharge
CatchmentTimeseries[,"estimated_groundwater_discharge"]   <-  apply(EigenTimeseries[,,"groundwater_discharge"],MARGIN=1,sum)
#if (!is.na(CatchmentTimeseries[1,"observed_groundwater_discharge"])){
#  CatchmentTimeseries[1,"estimated_groundwater_discharge"]    <-  CatchmentTimeseries[1,"observed_groundwater_discharge"]    #Initial estimated set to match observed depth
#}
CatchmentTimeseries[,"scaled_est_discharge"] <- CatchmentTimeseries[,"estimated_groundwater_discharge"] * DischargeScaleFactor

#Add an offset
CatchmentTimeseries[,"scaled_est_discharge"] <- pmax(CatchmentTimeseries[,"scaled_est_discharge"] - GWBypassFlow,0)


#Convert to a zoo data type
CatchmentTimeseries <- as.zoo(CatchmentTimeseries,as.Date(row.names(CatchmentTimeseries),format="%d/%m/%Y"))
} #end of function


