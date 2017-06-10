### Reading in and preparing the rainfall data from the stations in WA

Bencubbin <- read.csv("M:/Master thesis/Rainfall Data Bencubbin WA/IDCJAC0009_10007_1800_Data.csv", header= T, sep = ",")
Bencubbin <- read.csv("M:/Master thesis/Rainfall Data Margaret River WA/IDCJAC0009_2019_1800_Data.csv", header= T, sep = ",")


millimeters <- Bencubbin$Rainfall.amount..millimetres.
barplot(millimeters)
rain_data <- ts(Bencubbin[3624:37406,])
barplot(rain_data)
plot.ts(rain_data)
components <- decompose(rain_data)
plot(components)


rain_data$Year <- as.Date( paste( rain_data$Year,  rain_data$Month , rain_data$Day, sep = ", " )  , format = "%y,%m,%d" )
head(rain_data)

millimeters[3624:37406]



