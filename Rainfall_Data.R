### Reading in and preparing the rainfall data from the stations in WA

# Bencubbin <- read.csv("M:/Master thesis/Rainfall Data Bencubbin WA/IDCJAC0009_10007_1800_Data.csv", header= T, sep = ",")
# Margaret <- read.csv("M:/Master thesis/Rainfall Data Margaret River WA/IDCJAC0009_2019_1800_Data.csv", header= T, sep = ",")
# 
# TennantPost <- read.csv("M:/Master thesis/tennant creek post/IDCJAC0009_015087_1800_Data.csv", header= T, sep = ",")

TennantAirport <- read.csv("M:/Master thesis/tennant creek airport/IDCJAC0009_015135_1800_Data.csv", header= T, sep = ",")


# Bencubbin <- na.omit(Bencubbin)
# rainfall <- c(Bencubbin$Rainfall.amount..millimetres.)
# 
# millimeters <- Bencubbin$Rainfall.amount..millimetres.
# barplot(millimeters)

# TennantPost <- na.omit(TennantPost)
# rainfall <- c(TennantPost$Rainfall.amount..millimetres.)
# 
# which(is.na(TennantPost))
# 
# millimeters <- TennantPost$Rainfall.amount..millimetres.
# barplot(millimeters)


TennantAirport <- na.omit(TennantAirport)
rainfall <- c(TennantAirport$Rainfall.amount..millimetres.)
rainfall

rainfallwhich(is.na(TennantAirport))

millimeters <- TennantAirport$Rainfall.amount..millimetres.
barplot(millimeters)

# # Convert it to a time series object.
rainfall.timeseries <- ts(rainfall, frequency = (365)) #start = c(2012,1),frequency = 12
# 
# # Print the timeseries data.
# print(rainfall.timeseries)
# 
# # Give the chart file a name.
# png(file = "rainfall.png")
# 
# # Plot a graph of the time series.
# plot(rainfall.timeseries)
# 
# # Save the file.
# dev.off()
acf(rainfall.timeseries)

components <- decompose(rainfall.timeseries)
plot(components)
plot(stl(rainfall.timeseries, "per"))

install.packages("forecast")
library(forecast)
trend = ma(rainfall.timeseries, order = 4, centre = F)
plot(as.ts(rainfall.timeseries))
lines(trend, col="red")
plot(as.ts(trend))

detrended <- as.ts(rainfall.timeseries - trend)
plot(detrended)

m_beer = t(matrix(data = detrend_beer, nrow = 4))
seasonal_beer = colMeans(m_beer, na.rm = T)
plot(as.ts(rep(seasonal_beer,16)))

# Bencubbin_rain <- data.frame(Bencubbin$Year, Bencubbin$Month, Bencubbin$Day, Bencubbin$Rainfall.amount..millimetres.)
# Margaret_rain <- data.frame(Margaret$Year, Margaret$Month, Margaret$Day, Margaret$Rainfall.amount..millimetres.)

# barplot(Bencubbin_rain, main = "Bencubbin, WA")
# barplot(Margaret_rain, main = "Margaret, WA")
# attach(Margaret)
# Margaret$date <- as.Date(with(Margaret, paste(Year, Month, Day,sep="-")), "%Y-%m-%d")
# 
# ts_Margaret_rain <- data.frame(Margaret$Rainfall.amount..millimetres.)

ts_Margaret_rain <-ts(Margaret_rain)
plot.ts(ts_Margaret_rain)
components <- decompose(ts_Margaret_rain)
plot(components)
head(Margaret)
stl(as.ts(Margaret$Rainfall.amount..millimetres.,"per"))
plot(stl(as.ts(ts_Margaret_rain,"per")))

ts_Bencubbin_rain <- ts(Bencubbin$Rainfall.amount..millimetres.)
plot.ts(ts_Bencubbin_rain)

require(graphics)

plot(stl(nottem, "per"))
plot(stl(nottem, s.window = 7, t.window = 50, t.jump = 1))

plot(stllc <- stl(log(co2), s.window = 21))
summary(stllc)
## linear trend, strict period.
plot(stl(log(co2), s.window = "per", t.window = 1000))

## Two STL plotted side by side :
stmd <- stl(mdeaths, s.window = "per") # non-robust
summary(stmR <- stl(mdeaths, s.window = "per", robust = TRUE))
op <- par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0), mfcol = c(4, 2))
plot(stmd, set.pars = NULL, labels  =  NULL,
     main = "stl(mdeaths, s.w = \"per\",  robust = FALSE / TRUE )")
plot(stmR, set.pars = NULL)
# mark the 'outliers' :
(iO <- which(stmR $ weights  < 1e-8)) # 10 were considered outliers
sts <- stmR$time.series
points(time(sts)[iO], 0.8* sts[,"remainder"][iO], pch = 4, col = "red")
par(op)   # reset


Margaret$Month <- month.abb[Margaret$Month]
Margaret$Year <- as.character(Margaret$Year)
Margaret$Day <- as.character(Margaret$Day)
Margaret$Date <- as.Date( paste( Margaret$Year,  Margaret$Month , Margaret$Day, sep = "," ) , format = "%y,%m,%d")
# head(Margaret)
# 
# millimeters[3624:37406]
# month.abb[Margaret$Month]
# 


