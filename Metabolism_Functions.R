fillHoles <- function(dataIn,maxLength,timeStep)
  #Linearly interpolate values for strings (up to specified length) of missing data
  #CTS 31 Jul 2009
  
  #Inputs
  # dataIn:     A data.frame with two columns, first is "datetime", second is data
  # maxLength:  Maximum length of NA string that you are willing to interpolate across.
  #             NOTE that this is given in minutes
  # timeStep:   The time step of the data
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


aggregate_metab <- function(x,width=3){
  ####Code for aggregating data
  # Inputs
  # x: vector of data that you want to aggregate over a rolling time period
  # width: number of observations to aggregate over
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
merge_met_data <- function(buoy,airport,time_step, startdate, enddate){
  #buoy is a two column data frame of buoy data. column 1 is datetime,
  #column 2 is met data
  #airport is a two column data frame of airport met data. Column 1 is datetime,
  #   column 2 is met data to match buoy
  #time_step is the timestep of buoy data in minutes
  #start and enddates are the start and end date for the data being merged
  # defaults to range of buoy data
  ###Fix for requiring lubridate and tidyverse
  
  packages = c("tidyverse","lubridate")
  package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })
  
  #get start and end dates if they are not provided. defaults to range of buoy
  #data.
  if(is.null(startdate)) startdate <- as_date(min(buoy$datetime,na.rm=TRUE))
  if(is.null(enddate)) startdate <- as_date(min(buoy$datetime,na.rm=TRUE))
  
  
  dt_seq = data.frame(datetime = seq(
    from=as_datetime(paste(startdate," 00:00:00")),
    to=as_datetime(paste(enddate," 23:59:00")),
    by=time_step*60))
  
  dat = dt_seq %>% 
    left_join(buoy) %>% 
    left_join(airport,by = "datetime")
  
  dat[[2]][which(is.na(dat[[2]]))] <- dat[[3]][which(is.na(dat[[2]]))]
  
  dat <- dat[,c(1:2)]
  names(dat) = names(buoy)
  return(dat)
}

######
# Determine Time Frequence of data
######
time.freq <- function(x,obs.freq=c(1,10,30)){
  x_diff <- abs(as.numeric(difftime(x[-length(x)],x[-1],units="mins")))
  freq.dat <- data.frame(datetime=x,time_freq = c(x_diff,NA))
  freq.dat <- freq.dat %>% 
    mutate(time_freq=replace(time_freq,which(!(time_freq %in% obs.freq)),NA)) %>%
    drop_na() %>% 
    group_by(sampledate=as_date(datetime)) %>% 
    summarize(time.freq = paste(unique(time_freq,na.rm=TRUE),collapse = ","))
  print(freq.dat[which(!(freq.dat$time.freq %in% obs.freq)),])
  return(freq.dat)
}

######
# Exract aggregated data
######
extract.metab.data <- function(dataIn, time.freq){
  #DataIn: output from aggregate.metab function
  #time.freq: time frequency in minutes of data aggregation from aggregate.metab
  start.datetime <- min(floor_date(dataIn$datetime,unit="hours"))
  end.datetime <- max(ceiling_date(dataIn$datetime,unit="hours"))
  var.names <- names(dataIn)
  dt_seq = data.frame(datetime = seq(
    from=as_datetime(start.datetime),
    to=as_datetime(end.datetime),
    by=time.freq*60))
  dataIn <- dataIn %>% right_join(dt_seq) %>% select(datetime,val.30)
  names(dataIn)[2] <- var.names[2]
  return(dataIn)
}

#function to remove points 
identifyPch <- function(x, y = NULL, n = 2, pch = 19, ...)
{
  xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
  sel <- rep(FALSE, length(x)); res <- integer(0)
  while(sum(sel) < n) {
    ans <- identify(x[!sel], y[!sel], n = 1, plot = FALSE, ...)
    if(!length(ans)) break
    ans <- which(!sel)[ans]
    points(x[ans], y[ans], pch = pch)
    sel[ans] <- TRUE
    res <- c(res, ans)
  }
  res
}
complete.ts <- function(dataIn,time.freq) {
  #Function creates complete ts without missing time periods
  #DataIn: two colomn dataframe. 1st column datetime, 2nd column data
  #time.freq: time frequency of data
  start.datetime <- min(floor_date(dataIn$datetime,unit="hours"))
  end.datetime <- max(ceiling_date(dataIn$datetime,unit="hours"))
  dt_seq = data.frame(datetime = seq(
    from=as_datetime(start.datetime),
    to=as_datetime(end.datetime),
    by=time.freq*60))
  dataIn <- dataIn %>% right_join(dt_seq)
  return(dataIn)
}
