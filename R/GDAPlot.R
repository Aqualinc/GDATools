#' A simple plot function
#'
#' This function lets you plot the estimated groundwater response variable (water depth or discharge) against the observed
#' and provides a bunch of objective function statistics.
#' It includes the option for start and end times for the plot.
#' @param EstimatedTimeseries
#' @param ObservedTimeseries
#' @param StartDate A date in YYYY-MM-dd format. Defaults to the minimum of the input timeseries'
#' @param EndtDate A date in YYYY-MM-dd format. Defaults to the maximum of the input time series'
#' @keywords groundwater, hydrology
#' @export
#' @examples
#' fishRecharge<-vadose.recharge(ZoneStorageTime=c(4,3,3,3),AquiferZoneDryFractions = list(c(0.8,0.8,0.8,0.8),c(0,0,0,0),c(0,0,0.193,0.8263)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.0069,0.1737)),RiverRechargeFractions=c(1,1,0.958,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
#' fishDischarge <- bouss.eigen(WellDistance=66000,ZoneLengths=c(35600,5700,8900,15800),Storativity=0.01,Transmisivity=313632,DischargeScaleFactor=950,RechargeData=fishRecharge,GWBypassFlow=1)
#' Observed <- read.csv("Observeddata.csv")
#' library(zoo)
#' ObservedZoo <- read.zoo(Observeddata,format="%d/%m/%Y")
#' bouss.plot(fishDischarge[,4],ObserveddataZoo$FishCreekDischarge)

bouss.plot <- function(EstimatedTimeseries,ObservedTimeseries,startDate=NULL,endDate=NULL) {
  library(hydroGOF)
  
  #  observed   <- timeseries[,4]
  #  estimated  <- timeseries[,4]
  Estimated <- window(EstimatedTimeseries,start=startDate,end=endDate)
  Observed <- window(ObservedTimeseries,start=startDate,end=endDate)
  ggof(Estimated,Observed,lwd=c(1.5,1.5),pch=c(".","."),lty=c(1,1),col=c("red","black"))
  
}