#This file is used to scrape USGS Stream Gauge data for the Excel model if you want a new date range of USGS Gauge data or if you go to a new location.
#The script scrapes the Gauge data and saves it as a csv. You open the csv manually and copy and paste the data from there into the Excel model you are working with. 
#You put it in the Flow_data sheet in column M starting M10 without the column header.

library(dplyr)
library(EGRET)
library(lubridate)

path<- "C:/Users/jcrane/repos/combined-water-balance/WaterBalance"
codePath <- file.path(path, 'Code')
outPath <- file.path(path,"Output")
dataPath <- file.path(path, "Data")

startDate <- "1979-01-01"
endDate <- "2023-06-24"
siteID <- "09505200"
parameterCd = "00060"#denotes pulling discharge
siteName<- "WET BEAVER CREEK NEAR RIMROCK, AZ"

Daily <- EGRET::readNWISDaily(siteNumber = siteID, parameterCd = parameterCd, 
                              startDate = startDate, endDate = endDate) |>
  dplyr::filter(grepl('A', Qualifier)) |> #this filters for any Qualifier that has an A. It will return A and A:E
  dplyr::mutate(CFS = Q*35.314666212661) #converts Q to flow cfs
Meas<- Daily[,c("Date", "CFS")]

#Fill Missing values of Gauge data with NAs by creating a date sequence to merge the gauge data with
DateSeq<- data.frame(Date = seq.Date(as.Date(startDate), as.Date(endDate), by="day"))
MeasFullDates = dplyr::full_join(Meas, DateSeq)
MeasFullDates = MeasFullDates[order(MeasFullDates$Date), ]

#write to Excel
write.csv(MeasFullDates, file.path(dataPath, paste(siteName, "USGS.csv")))






