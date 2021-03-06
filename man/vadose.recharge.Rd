% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Vadose.R
\name{vadose.recharge}
\alias{vadose.recharge}
\title{A Vadose recharge time series generator.}
\usage{
vadose.recharge(ZoneStorageTime = c(30, 30, 30, 20),
  AquiferZoneDryFractions = list(c(0.3, 0.3, 0.3, 0.3), c(0, 0, 0.1, 0.1),
  c(0.991, 0.9398, 0.7, 0.7)), AquiferZoneIrrigFractions = list(c(0, 0, 0, 0),
  c(0, 0, 0, 0), c(0.009, 0.0602, 0.0069, 0.1737)),
  RiverRechargeFractions = c(0, 0, 0, 0),
  RechargeFileName = "GoldenBayLandSurfaceRechargeData.csv",
  PumpingFileName = "GoldenBayGWPumpingData.csv")
}
\arguments{
\item{ZoneStorageTime}{The time (in days) that it takes for water to get through the vadose zone into the groundwater, one value for each zone.}

\item{AquiferZoneDryFractions}{a list of vectors of the proportion of each aquifer's dryland recharge that is contributing to each zones vadose recharge}

\item{AquiferZoneIrrigFractions}{a list of vectors of the proportion of each aquifer's irrigation recharge that is contributing to each zones vadose recharge}

\item{RiverRechargeFractions}{The fraction of the river discahrge that contributes to the groundwater for each zone.}

\item{PumpingFileName}{The csv file with the pumping data in it. In mm. One series per aquifer per zone, ordered as Z1A1, Z1A2..Z1An,Z2A1,Z2A2..Z2An..ZnA1,ZnA2...ZnAn}

\item{RechargeFilename}{The csv file with all the data in it. This is the daily timeseries of vadose recharge for each zone, and river recharge in mm}
}
\description{
This function reads in a timeseries of surface recharge for zones in a catchment and estimates the timeseries of vadose recharge.
Equations are from:
Bidwell, V., Burbery, L., 2011. Groundwater Data Analysis - quantifying aquifer dynamics. Prepared for Envirolink Project 420-NRLC50 (No. 4110/1). Lincoln Ventures Ltd.
which in turn cites the following, though note that the symbol labels are different.
Bidwell, V.J., Stenger, R., Barkle, G.F., 2008. Dynamic analysis of groundwater discharge and partial-area contribution to Pukemanga Stream, New Zealand. Hydrol. Earth Syst. Sci. 12, 975-987. doi:10.5194/hess-12-975-2008
}
\examples{
fishRecharge<-vadose.recharge(ZoneStorageTime=c(4,3,3,3),AquiferZoneDryFractions = list(c(0.8,0.8,0.8,0.8),c(0,0,0,0),c(0,0,0.193,0.8263)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.0069,0.1737)),RiverRechargeFractions=c(1,1,0.958,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
}
\keyword{groundwater,}
\keyword{hydrology}

