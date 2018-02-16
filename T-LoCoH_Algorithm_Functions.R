# Function for creating training-testing datasets (with approximately 80% training)
# seed argument allows for recreation of the same dataset for the sake of replication
train.test <- function (data, seed=1) {
  df <- data.frame(matrix(TRUE,nrow(data),100))
  set.seed(seed)
  
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
  return(df)
}

# Function for performing the actual algorithm based on:
# a) train.test - the result of the previous function
# b) k.vals - a vector of k.vals over which you would like to search
# c) s.max - the maximum value of s that the grid-based search will cover
# d) num.s.vals - number of equal intervals to insert between 0 and s.max
# e) data.lxy - a properly formatted lxy-object from tlocoh
# f) data - the original dataset (xy-coords are sufficient)
algo <- function(train.test, k.vals, s.max, num.s.vals, data.lxy, data) {
  
  # Determine number of test points across all 100 train/test splits for BIC calculation
  count <- 0
  for (i in 1:ncol(train.test)) {
    new <- table(train.test[,i])[1]
    count <- count + new
  }
  
  trace <- data.frame(matrix(0,length(k.vals),7))
  
  for (k in 1:length(k.vals)) {
    for (z in 1:(num.s.vals + 1)) {
      
      current.k_val <- k.vals[k]
      current.s_val <- (z-1) * (s.max/num.s.vals)
      
      # Calculate the nearest neighbors and create lhs object for full dataset
      full.lxy <- lxy.nn.add(data.lxy, k=current.k_val, s=current.s_val, status=F)
      full.lhs <- lxy.lhs(full.lxy, k=current.k_val, s=current.s_val, status=F)
      coords <- full.lhs[[1]]$pts@coords
      
      # Create list for the negative log likelihood values for each test/train split in df
      likelihood <- list()
      
      for (j in 1:ncol(train.test)) {
        
        # Create a one-column data frame from df
        df1 <- train.test[1:nrow(train.test),j]
        
        # Create selection of hulls based on Boolean
        hulls.sel.idx <- which(df1)
        full.hulls <- hulls(full.lhs)
        selected.hulls <- full.hulls[[1]] [ full.hulls[[1]]@data$pts.idx %in% hulls.sel.idx , ]
        
        # Determine the number of points in training dataset
        df.temp <- data.frame(as.numeric(train.test[,j]))
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
        
        # Calculate likelihood
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
      
      log.like <- data.frame(matrix(0,length(likelihood),4))
      for (i in 1:length(likelihood)) {
        log.like[i,1] <- sum(likelihood[[i]]$loglike, na.rm=TRUE)
        log.like[i,2] <- sum(likelihood[[i]]$ln.like, na.rm=TRUE)
        log.like[i,3] <- -2*(log.like[i,2]) + 2*k
        log.like[i,4] <- -2*(log.like[i,2]) + k*log(as.numeric(count),exp(1))
      }
      
      # Assign mean of the means of negative log likelihoods as the current posterior probability for Metropolis-Hastings
      new.postLike <- sum(log.like[,1])
      new.lnLike <- sum(log.like[,2])
      new.AIC <- sum(log.like[,3])
      new.BIC <- sum(log.like[,4])
      
      #cat("iteration:", counter, "chain:", current.k_val, current.s_val, "likelihood:", new.postLike, "\n")
      
      trace[((k-1)*(num.s.vals+1))+z,1] <- current.s_val
      trace[((k-1)*(num.s.vals+1))+z,2] <- current.k_val
      trace[((k-1)*(num.s.vals+1))+z,3] <- hull.mean
      trace[((k-1)*(num.s.vals+1))+z,4] <- new.postLike
      trace[((k-1)*(num.s.vals+1))+z,5] <- new.lnLike
      trace[((k-1)*(num.s.vals+1))+z,6] <- new.AIC
      trace[((k-1)*(num.s.vals+1))+z,7] <- new.BIC
      colnames(trace) <- c("s.val", "k.val", "hull.mean", "postLike", "lnLike", "AIC", "BIC")
      
    } # End of s loop
    
  } # End of k loop
  
  return(trace) 
}

###############################

#Code for visualizing the results (which are saved in an object called 'trace')
library(lattice)

my.palette <- c("#FFF5F0", "#FEE0D2", "#FEE0D2", "#FCBBA1", "#FCBBA1", "#FC9272", "#FC9272", "#FB6A4A", "#FB6A4A", "#EF3B2C", "#EF3B2C", "#CB181D", "#CB181D", "#A50F15", "#A50F15", "#67000D")

levelplot(-trace$X7 ~ trace$X2 * trace$X3, xlab="k value", ylab="s value", zlab="IC", screen = list(z = -30, x=-60), regions=TRUE, cuts=15, col.regions=my.palette)

wireframe(-trace$X7 ~ trace$X2 * trace$X3, xlab="k value", ylab="s value", zlab="IC", screen = list(z = -30, x=-60), drape=TRUE, cuts=15, col.regions=my.palette)
