library(RCurl)
library(RJSONIO)
library(plyr)

getDistance <- function(origin, destination, key, return.call = "json") {
	# Input:  origin, destination city 
	# Output: duration (hours), distance of travel (kilometers)
	# inspired by http://www.r-bloggers.com/using-google-maps-api-and-r/
	
	root <- "https://maps.googleapis.com/maps/api/distancematrix/"
	address <- paste(root, return.call, "?", "origins=", origin, "&destinations=", destination,"&key=", key, sep = "")
	doc <- getURL(address)
	data <- fromJSON(doc, simplify = FALSE)
	if(data$status=="OK") {
		distance <- data$rows[[1]]$elements[[1]]$distance$value/(1000) # from meters to kilometers
		duration <- data$rows[[1]]$elements[[1]]$duration$value/(360)  # from seconds to hours
		return(c(distance, duration))
	} else {
		warning(paste('Didn\'t get data for ',origin,' => ',destination,'. Returning NAs!', sep = ''))
		return(c(NA,NA))
	}
}

source('googleKey.r',echo = F)