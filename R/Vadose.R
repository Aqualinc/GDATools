#' A Vadose recharge time series generator.
#'
#'This function reads in a timeseries of surface recharge for zones in a catchment and estimates the timeseries of vadose recharge.
# Each zone has recharge for dryland, irrigated land and from the river.
# The dryland and irrigated land can be recharging to different aquifers which contribute different amounts to the groundwater that is being modelled.
# The vadose recharge is assumed to be well estimated by an exponential weighted moving average.
#'Equations are from:
#' Bidwell, V., Burbery, L., 2011. Groundwater Data Analysis - quantifying aquifer dynamics. Prepared for Envirolink Project 420-NRLC50 (No. 4110/1). Lincoln Ventures Ltd.
#' which in turn cites the following, though note that the symbol labels are different.
#' Bidwell, V.J., Stenger, R., Barkle, G.F., 2008. Dynamic analysis of groundwater discharge and partial-area contribution to Pukemanga Stream, New Zealand. Hydrol. Earth Syst. Sci. 12, 975-987. doi:10.5194/hess-12-975-2008
#' @param ZoneStorageTime The time (in days) that it takes for water to get through the vadose zone into the groundwater, one value for each zone.
#' @param AquiferZoneDryFractions a list of vectors of the proportion of each aquifer's dryland recharge that is contributing to each zones vadose recharge
#' @param AquiferZoneIrrigFractions a list of vectors of the proportion of each aquifer's irrigation recharge that is contributing to each zones vadose recharge
#' @param RiverRechargeFractions The fraction of the river discahrge that contributes to the groundwater for each zone.  
#' @param RechargeFilename The csv file with all the data in it. This is the daily timeseries of vadose recharge for each zone, and river recharge in mm
#' @param PumpingFileName The csv file with the pumping data in it. In mm. One series per aquifer per zone, ordered as Z1A1, Z1A2..Z1An,Z2A1,Z2A2..Z2An..ZnA1,ZnA2...ZnAn
#' @keywords groundwater, hydrology
#' @export
#' @examples
#' fishRecharge<-vadose.recharge(ZoneStorageTime=c(4,3,3,3),AquiferZoneDryFractions = list(c(0.8,0.8,0.8,0.8),c(0,0,0,0),c(0,0,0.193,0.8263)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.0069,0.1737)),RiverRechargeFractions=c(1,1,0.958,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")


vadose.recharge <- function(ZoneStorageTime = c(30,30,30,20),
                            AquiferZoneDryFractions = list(c(0.3,0.3,0.3,0.3),
                                                              c(0,0,0.1,0.1),
                                                              c(0.991,0.9398,0.7,0.7)),
                            AquiferZoneIrrigFractions = list(c(0,0,0,0),
                                                              c(0,0,0,0),
                                                              c(0.009,0.0602,0.0069,0.1737)),
                            RiverRechargeFractions = c(0,0,0,0),
                            RechargeFileName ="GoldenBayLandSurfaceRechargeData.csv",
                            PumpingFileName = "GoldenBayGWPumpingData.csv")  
  
{  #Start of the function
  
#Load libraries
library(TTR)          #This library includes the EMA exponentially weighted moving average

#read in the data
ZoneTimeseries        <-  read.csv(RechargeFileName)                              #this is the daily timeseries of observed discharge and groundwater level, vadose recharge for each zone, river recharge to groundwater and groundwater pumping in mm
PumpingTimeseries     <-  read.csv(PumpingFileName)                               #This is the daily timeseries of pumping for each aquifer in each zone. Same date range as the Recharge data.

#Calculate the number of zones and aquifers based on the 
NumberOfZones         <- length(AquiferZoneDryFractions[[1]])
NumberOfAquifers      <- length(AquiferZoneDryFractions)

#Multiply the aquifer recharge fractions by the recharge timeseries to calculate a total recharge for each zone
ZoneFractions  <- c()                                               #Initialise ZoneFractions to an empty set
for (ZoneNo in 1:NumberOfZones){
  
  for (AquiferNo in 1:NumberOfAquifers){
    ZoneFractions <- c(ZoneFractions,AquiferZoneDryFractions[[AquiferNo]][ZoneNo],AquiferZoneIrrigFractions[[AquiferNo]][ZoneNo])
  }
}

#Multiply the recharge timeseries by their respective fractions
ZoneSurfaceRecharge  <- sweep(ZoneTimeseries[2:(ncol(ZoneTimeseries)-2*NumberOfZones)],MARGIN=2,ZoneFractions,'*')
ZoneRiverRecharge <- sweep(ZoneTimeseries[(ncol(ZoneTimeseries)-2*NumberOfZones+1):(ncol(ZoneTimeseries)-NumberOfZones)],MARGIN=2,RiverRechargeFractions,'*')

#Get the sums for each zone so that we are left with just one timeseries for each zone
ZoneSurfaceRechargeTotals    <- c()     #Initialise
for (ZoneNo in 1:(NumberOfZones)){
  ZoneSurfaceRechargeTotals <- cbind(ZoneSurfaceRechargeTotals,rowSums(ZoneSurfaceRecharge[,c(((ZoneNo-1)*(2*NumberOfAquifers)+1):(ZoneNo*2*NumberOfAquifers))]))
}

#Uncomment below and comment out the line below below if the river recharge is to have a delay into the vadose zone
##Add in the river recharge
#ZoneRechargeTotals <- ZoneSurfaceRechargeTotals + ZoneRiverRecharge
ZoneRechargeTotals <- ZoneSurfaceRechargeTotals

#build a list of lists of the recharge data with the weighting time constant, needed to apply the exp weighted moving average
RechargeList <- list()
for (ZoneNo in 1:(NumberOfZones)){
  RechargeList[[length(RechargeList)+1]]<-list(series=ZoneRechargeTotals[,ZoneNo],timeConstant=ZoneStorageTime[ZoneNo])
}

#Apply the exponential weighted moving average
ZoneVadoseRecharge <- sapply(RechargeList, function(x) EMA(x$series,n=1,ratio=1-exp(-1/max(x$timeConstant,0.000001))))
names(ZoneVadoseRecharge) <- paste0("VadoseZ",seq(1,length.out=NumberOfZones))

#This recharge is for when the stream is directly connected to the groundwater, so no delay.
#Add in the river recharge
ZoneVadoseRecharge <- ZoneVadoseRecharge + ZoneRiverRecharge

#Take out any pumping
#Somehow get the proportions of each aquifer for each zone for pumping.
#Add together the dry land and irrigated land fractions to get aquifer fractions
#Use these aquifer fractions to get total pumping
AquiferFractions <- sapply(as.matrix(AquiferZoneDryFractions), unlist) + sapply(as.matrix(AquiferZoneIrrigFractions),unlist)

PumpingScaled <- sweep(PumpingTimeseries[,2:ncol(PumpingTimeseries)],MARGIN=2, as.vector(t(AquiferFractions)),'*')
#Now I need to add together the aquifers
PumpingTotals    <- c()     #Initialise
for (ZoneNo in 1:(NumberOfZones)){
  PumpingTotals <- cbind(PumpingTotals,rowSums(PumpingScaled[,c(((ZoneNo-1)*(NumberOfAquifers)+1):(ZoneNo*NumberOfAquifers))]))
}

ZoneVadoseRecharge <- ZoneVadoseRecharge - PumpingTotals

#Convert from millimetres to metres
ZoneVadoseRecharge<- ZoneVadoseRecharge / 1000

#Put the dates back on
rownames(ZoneVadoseRecharge) <- ZoneTimeseries[,1]

#Rename the columns
colnames(ZoneVadoseRecharge) <- paste0("Zone",c(1:NumberOfZones))

#Put the observed data back on
#ZoneVadoseRecharge <- cbind(ZoneTimeseries[,2:3],ZoneVadoseRecharge) 

return(ZoneVadoseRecharge)
} #End of function