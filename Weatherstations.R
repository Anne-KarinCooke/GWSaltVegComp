setwd("M:/Master thesis")

##CSV1
stations1 <-read.table("DC02D_StnDet_999999999142658.txt", sep = ",", header=T)
#stations1
head(stations1)
attach(stations1)
complete1 <- subset(stations1, stations1$Percentage.complete.between.first.and.last.records > 90)
head(complete1)

duration1 <- subset(complete1, complete1$Last.year.of.data.supplied.in.data.file - complete1$First.year.of.data.supplied.in.data.file >= 100)
head(duration1)
duration1

##CSV2
stations2 <-read.table("DC02D_StnDet_999999999142659.txt", sep = ",", header=T)
#stations2
head(stations2)
attach(stations2)
complete2 <- subset(stations2, stations1$Percentage.complete.between.first.and.last.records > 90)
head(complete2)

duration2 <- subset(complete2, complete2$Last.year.of.data.supplied.in.data.file - complete2$First.year.of.data.supplied.in.data.file >= 100)
head(duration2)
duration2


##CSV3
stations3 <-read.table("DC02D_StnDet_999999999142663.txt", sep = ",", header=T)
#stations3
head(stations3)
attach(stations3)
complete3 <- subset(stations3, stations1$Percentage.complete.between.first.and.last.records > 90)
head(complete3)

duration3 <- subset(complete3, complete3$Last.year.of.data.supplied.in.data.file - complete3$First.year.of.data.supplied.in.data.file >= 100)
head(duration3)
duration3


##CSV4
stations4 <-read.table("DC02D_StnDet_999999999142664.txt", sep = ",", header=T)
#stations4
head(stations4)
attach(stations4)
complete4 <- subset(stations4, stations1$Percentage.complete.between.first.and.last.records > 90)
head(complete4)

duration4 <- subset(complete4, complete4$Last.year.of.data.supplied.in.data.file - complete4$First.year.of.data.supplied.in.data.file >= 100)
head(duration4)
duration4

mapCountryData("Australia")
data(countryExData)
install.packages("ggmap")
library(ggmap)
map <- get_map(location = 'Australia')

?get_map()

plot(duration4$Latitude.to.4.decimal.places.in.decimal.degrees, duration4$Longitude.to.4.decimal.places.in.decimal.degrees)
plot(duration3$Latitude.to.4.decimal.places.in.decimal.degrees, duration3$Longitude.to.4.decimal.places.in.decimal.degrees)
plot(duration2$Latitude.to.4.decimal.places.in.decimal.degrees, duration2$Longitude.to.4.decimal.places.in.decimal.degrees)
plot(duration1$Latitude.to.4.decimal.places.in.decimal.degrees, duration1$Longitude.to.4.decimal.places.in.decimal.degrees)


library(maps)       # Provides functions that let us plot the maps
library(mapdata)    # Contains the hi-resolution points that mark out the countries
#points(-1.615672,54.977768,col=2,pch=18)

map('worldHires','Australia')
rasterize(duration4$Latitude.to.4.decimal.places.in.decimal.degrees)
duration4$Longitude.to.4.decimal.places.in.decimal.degrees
rasteriz


# seasonalRain <-read.table("seasrain.txt", sep = ",", header=T)
# plot(raster(seasonalRain))

head(seasonalRain)
plot(seasonalRain)


df1<- as.data.frame(as.numeric(duration4$Latitude.to.4.decimal.places.in.decimal.degrees), as.numeric(duration4$Longitude.to.4.decimal.places.in.decimal.degrees))
df1

# loading the required packages
library(ggplot2)
library(ggmap)

# creating a sample data.frame with your lat/lon points
lon <- c(-38.31,-35.5)
lat <- c(40.96, 37.5)
df <- as.data.frame(cbind(lon,lat))

# getting the map
mapgilbert <- get_map(location = c(lon = mean(df$lon), lat = mean(df$lat)), zoom = 4,
                      maptype = "satellite", scale = 2)

# plotting the map with some points on it
ggmap(mapgilbert) +
  geom_point(data = df, aes(x = lon, y = lat, fill = "red", alpha = 0.8), size = 5, shape = 21) +
  guides(fill=FALSE, alpha=FALSE, size=FALSE)
