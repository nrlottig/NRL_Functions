#Linearly interpolate values for strings (up to specified length) of missing data
#CTS 31 Jul 2009

#Inputs
# dataIn:     A data.frame with two columns, first is "datetime", second is data
# maxLength:  Maximum length of NA string that you are willing to interpolate across.
#             NOTE that this is given in minutes
# timeStep:   The time step of the data


fillHoles <- function(dataIn,maxLength,timeStep)
{
  
  #Number of rows in dataIn
  nObs <- dim(dataIn)[1]
  
  #Express maxLength as number of time steps instead of number of minutes
  maxLength <- maxLength/timeStep
  
  #Temporarily replace NAs in data with -99999
  whichNA <- which(is.na(dataIn[,2]))
  dataIn[whichNA,2] <- -99999
  #Identify strings of NA (-99999) values
  rleOut <- rle(dataIn[,2])
  which9 <- which(rleOut$values==-99999)
  
  #If no NA strings in data, return dataIn
  if (length(which9)==0)
  {
    return(dataIn)
  } else
    
    #Otherwise, continue
  {
    
    #Interpolate valus for each string of NA values
    for (i in 1:length(which9))
    {
      
      #Determine start and end index of string i, and calculate length of string
      if (which9[i]==1)
      {
        stringStart <- 1
      } else
        
      {
        stringStart <- 1 + sum(rleOut$lengths[1:(which9[i]-1)])
      }
      
      stringEnd <- sum(rleOut$lengths[1:which9[i]])
      stringLength <- stringEnd-stringStart+1
      
      #Skip to next i if:
      #  -length of string exceeds maxLength,
      #  -stringStart is the first obs in dataIn
      #  -stringEnd is the last obs in dataIn
      if (stringLength > maxLength | stringStart==1 | stringEnd==nObs) next else
      {
        
        #Linearly interpolate missing values
        interp <- approx(x=c(dataIn[stringStart-1,"datetime"],dataIn[stringEnd+1,"datetime"]),
                         y=c(dataIn[stringStart-1,2],dataIn[stringEnd+1,2]),
                         xout=dataIn[stringStart:stringEnd,"datetime"],
                         method="linear")
        dataIn[stringStart:stringEnd,2] <- interp$y
      }
    }
    
    #Replace any remaining -99999 with NA
    dataIn[which(dataIn[,2]==-99999),2] <- NA
    
    #Return result
    return(dataIn)
  }
  
}

####Code for aggregating data
# Inputs
# x: vector of data that you want to aggregate over a rolling time period
# width: number of observations to aggregate over
aggregate_metab <- function(x,width=3){
  packages = c("zoo")
  package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })
  aggregate_odd <- function(x,width) {
    rollapply(x, 
              width=width,
              FUN=function(x) mean(x,na.rm=T), 
              by=1, 
              by.column=TRUE, 
              partial=TRUE, 
              fill=NA, 
              align="center")
  }
#aggregate data if width is an even number  
  aggregate_even <- function(x,width) {
    rollapply(x, 
              width=width, 
              FUN=function(x) mean(x,na.rm=T), 
              by=1, 
              by.column=TRUE, 
              partial=TRUE, 
              fill=NA, 
              align="center")
  }
#shift the data to center aggregation on an even width  
  shift <- function(x, lag) {
    n <- length(x)
    xnew <- rep(NA, n)
    if (lag < 0) {
      xnew[1:(n-abs(lag))] <- x[(abs(lag)+1):n]
    } else if (lag > 0) {
      xnew[(lag+1):n] <- x[1:(n-lag)]
    } else {
      xnew <- x
    }
    return(xnew)
  }
#calculate centered aggregation on even number 
  aggregate_even_center <- function(x,width) {
    rowMeans(cbind(aggregate_even(x,width),shift(aggregate_even(x,width),1)),na.rm=T)
  } 
  
  
  is_odd = width %% 2 # check if the rolling width is an odd number
  if(is_odd == 1) {
    x_aggregated <- aggregate_odd(x=x,width=width)
  } else {
    x_aggregated <- aggregate_even_center(x = x,width = width)
  }
  return(x_aggregated)
}




######This code is for centering aggregated data when period of aggregation is a evan number
######Essentially it weights tail values by 1/2 the other values in the aggregation


