library(rgdal)
library(sp)
library(maptools)
library(tlocoh)
library(foreach)
library(doParallel)
library(forecast)

#Create a list of .csv file names to import for parallel processing
name.list <- c("AAA", "BBB", "CCC", "DDD")

# Assign n-1 of the total cores to the cluster operation
num_cores <- detectCores()
cl<-makeCluster(num_cores)
registerDoParallel(cl)

#Create object to store elements made during parallel processing (using foreach loop)
track <- foreach(p = 1:length(name.list), .packages=c("tlocoh", "sp", "maptools", "forecast"), .errorhandling = 'remove', .combine='rbind') %dopar% {

#Import data from appropriate file path, exluding the actual file name
data <- read.csv(paste("/global/home/users/", name.list[p], ".csv", sep=''))
  
#Redefine the time signature
data$Datetime <- as.POSIXct(data$date, tz="GMT", origin = "1970-01-01 00:00:00")
full <- data.frame(cbind(data$x, data$y, data$Datetime))
colnames(full) <- c('long', 'lat', 'datetime')
full$datetime <- as.POSIXct(full$datetime, tz="GMT", origin = "1970-01-01 00:00:00")

# Two alternatives here: the first is to remove all NAs entirely from the dataset
#full1 <- full[!is.na(full[,1]),]

# The second is to use a Kalman smoother to fill in the NAs
long <- full$long
z <- long
fit <- auto.arima(long)
kr <- KalmanSmooth(long, fit$model)
id.na <- which(is.na(long))
num <- ncol(kr$smooth)
for (j in id.na) {
  z[j] <- kr$smooth[j,num]
}
long <- z

lat <- full$lat
z <- lat
fit <- auto.arima(lat)
kr <- KalmanSmooth(lat, fit$model)
id.na <- which(is.na(lat))
num <- ncol(kr$smooth)
for (j in id.na) {
  z[j] <- kr$smooth[j,num]
}
lat <- z

full2 <- cbind(long, lat)
full2 <- cbind(full2, data$Datetime)
colnames(full2) <- c('long', 'lat', 'datetime')
full2 <- as.data.frame(full2)
full2$datetime <- as.POSIXct(full2$datetime, tz="GMT", origin = "1970-01-01 00:00:00")

# Convert to UTM (with the Kalman smoothed dataset), adjust UTM Zone as necessary
full.sp.latlong <- SpatialPoints(full2[ , c("long","lat")], proj4string=CRS("+proj=longlat +ellps=WGS84"))
full.sp.utm <- spTransform(full.sp.latlong, CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
full.mat.utm <- coordinates(full.sp.utm)
colnames(full.mat.utm) <- c("x","y")
full.lxy <- xyt.lxy(xy=full.mat.utm, dt=full$datetime, id=as.character(name.list[p]), proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"), dup.dt.check=FALSE)

# Create 100 training/testing splits
df <- data.frame(matrix(TRUE,nrow(data),100))

for (i in 1:ncol(df)) {
  samp <- sample(seq(1,nrow(data),1), round(0.002222*(nrow(data))))
  for (j in 1:length(samp)) {
    for (k in 1:nrow(df)) {
      if (k == samp[j]) {
        df[k,i] <- FALSE
      }
    }
  }
}

count <- 0
for (i in 1:ncol(df)) {
  new <- table(df[,i])[1]
  count <- count + new
}

for (i in 3:nrow(df)) {
  for (j in 1:ncol(df)) {
    if (df[i-1,j] == FALSE && df[i-2,j] != FALSE) {
      for (k in 0:98) {
        df[i+k,j] <- FALSE
      }
    }
  }
}

df <- df[1:nrow(data),]

# Create list for output of cluster operation
s.loops <- data.frame(matrix(0,1189,9))

for (k in 2:30) {
  for (z in 1:41) { 
    
    #Define s and k values to perform grid search
    current.s_val <- (z-1) * 0.025 
    current.k_val <- k
    
    # Calculate the nearest neighbors and create lhs object for full dataset
    full.lxy <- lxy.nn.add(full.lxy, k=current.k_val, s=current.s_val, status=F)
    full.lhs <- lxy.lhs(full.lxy, k=current.k_val, s=current.s_val, status=F)
    coords <- full.lhs[[1]]$pts@coords
    
    # Create list for the negative log likelihood values for each test/train split in df
    likelihood <- list()
    
    for (j in 1:ncol(df)) {
      
      # Create a one-column data frame from df
      df1 <- df[1:nrow(data),j]
      
      # Create selection of hulls based on Boolean
      hulls.sel.idx <- which(df1)
      full.hulls <- hulls(full.lhs)
      selected.hulls <- full.hulls[[1]] [ full.hulls[[1]]@data$pts.idx %in% hulls.sel.idx , ]
      
      # Determine the number of points in training dataset
      df.temp <- data.frame(as.numeric(df[,j]))
      colnames(df.temp) <- "subset"
      total.pts <- sum(df.temp)
      
      # Find middle points of testing data and define as -1
      for (i in 2:nrow(df.temp)) {
        if (df.temp[i,1] == 0 && df.temp[i-1,1] != 0 && df.temp[i-1,1] != -1) {
          df.temp[i+50,1] <- -1
        }
      }
      
      df.temp <- df.temp[1:nrow(data),]
      df.temp <- cbind(coords, df.temp)
      
      # Extract the middle points for testing and record associated coordinates
      test.pts <- data.frame()
      q = 1
      for (i in 1:nrow(df.temp)) {
        if (df.temp[i,3] == -1) {
          test.pts[q,1] <- df.temp[i,1]
          test.pts[q,2] <- df.temp[i,2]
          q = q + 1
        }
      }
      
      colnames(test.pts) <- c("x", "y")
      test.pts <- SpatialPoints(test.pts[ , c("x","y")], proj4string=CRS("+proj=utm +south +zone=35 +ellps=WGS84"))
      
      poly <- SpatialPolygons(selected.hulls@polygons, proj4string = CRS('+proj=utm +south +zone=35 +ellps=WGS84'))
      
      # Determine the number of hulls under each test point,
      overlay <- data.frame(matrix(0,length(test.pts@coords[,1]),1))
      
      for (i in 1:length(test.pts@coords[,1])) {
        overlay.list <- over(test.pts[i], poly, returnList=TRUE)
        overlay[i,1] <- length(overlay.list[[1]])
      }
      
      hull.mean <- mean(full.lhs[[1]]$hulls@data$area)
      
      #Calculate likelihood
      for (i in 1:nrow(overlay)) {
        overlay[i,2] <- overlay[i,1]/length(selected.hulls[[1]])
        overlay[i,3] <- -log(overlay[i,2])
        overlay[i,4] <- log(overlay[i,2],exp(1))
        if (overlay[i,1] == 0) {
          overlay[i,3] <- -log((1/nrow(data))/100)
          overlay[i,4] <- log((1/nrow(data))/100, exp(1))
          overlay[i,4] <- log((1/nrow(data))/100, exp(1))
        }
      }
      
      colnames(overlay) <- c("over", "prob", "loglike", "ln.like")
      
      # Add values to likelihood list
      likelihood[[j]] <- as.list(overlay)
    }
    
    #Likelihood calculation across subsets, column 3 will represent the AIC equivalent and column 4 is akin to the BIC
    log.like <- data.frame(matrix(0,length(likelihood),4))
    for (i in 1:length(likelihood)) {
      log.like[i,1] <- sum(likelihood[[i]]$loglike, na.rm=TRUE)
      log.like[i,2] <- sum(likelihood[[i]]$ln.like, na.rm=TRUE)
      log.like[i,3] <- -2*(log.like[i,2]) + 2*k
      log.like[i,4] <- -2*(log.like[i,2]) + k*log(as.numeric(count),exp(1))
    }
    
    #Sum across each form of likelihood calculation
    new.postLike <- sum(log.like[,1])
    new.lnLike <- sum(log.like[,2])
    new.AIC <- sum(log.like[,3])
    new.BIC <- sum(log.like[,4])
    
    s.loops[((k-2)*41)+z,1] <- z
    s.loops[((k-2)*41)+z,2] <- current.k_val
    s.loops[((k-2)*41)+z,3] <- current.s_val
    s.loops[((k-2)*41)+z,4] <- hull.mean
    s.loops[((k-2)*41)+z,5] <- new.postLike
    s.loops[((k-2)*41)+z,6] <- new.lnLike
    s.loops[((k-2)*41)+z,7] <- new.AIC
    s.loops[((k-2)*41)+z,8] <- new.BIC
    s.loops[((k-2)*41)+z,9] <- as.character(name.list[p])
    
  } # End of s loop
  
} # End of k loop

print(s.loops)

} # End of parallel p loop

write.csv(track, "Excess_Track_All.csv")

w <- 1
z <- 1189
for (i in 1:length(name.list)) {
  temp.track <- as.data.frame(track[w:z,])
  write.csv(temp.track, file = paste(name.list[i], "_Trace.csv", sep=''))
  w <- w + 1189
  z <- z + 1189
}

# End cluster operation
stopCluster(cl)

