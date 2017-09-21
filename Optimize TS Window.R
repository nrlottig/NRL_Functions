#Script to identify time windows and missing data rates for time series analysis of multiple different systems.

#Check to make sure you have all the needed packages installed
packages <- c("RCurl", "reshape2", "zoo")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

require(RCurl)

#Example data file take from Oliver et al. 2017 (chla data)
dat = read.csv(text=getURL("https://raw.githubusercontent.com/nrlottig/NRL_Functions/master/datafiles/data_optimize_ts_window.csv"), header=T)

#Function for identifying optimal/maximum time series and time periods of record based on two input criteria
#record.length is length or long-term record to be analyzed or vector of potential long-term record lengths
#missing.freq is the extent of missing data allowed (single value or vector of values) as percent of data (e.g., 10%)
TS.Length = function(data,record.length,missing.freq){
  require(reshape2)
  require(zoo)
  #Misc functions needed
  RL <- function(data){sum(!is.na(data))}
  count.records = function(data){sum(data>=missing.freq.t[j])}

  out.data = data.frame(matrix(NA, nrow=length(record.length), ncol=length(missing.freq)))
  row.names(out.data) = record.length
  colnames(out.data) = missing.freq

  for(x in 1:length(record.length)){ #loop through record lengths
    new.data = data.frame(matrix(NA, nrow=nrow(data), ncol=ncol(data)-record.length[x]+1))

    for(i in 1:nrow(data)){ #calculate number of observations within each time period
      new.data[i,1] = data[i,1]
      x.temp = as.numeric(data[i,2:ncol(data)])
      new.data[i,2:ncol(new.data)] = rollapply(x.temp,width = record.length[x],FUN=RL)
    }

    missing.freq.t = ceiling(record.length[x] - (missing.freq/100*record.length[x]))
    summary.data = data.frame(matrix(NA, nrow=length(missing.freq.t), ncol=ncol(new.data)))

    for(j in 1:length(missing.freq.t)){
      summary.data[j,1] = missing.freq.t[j]
      if(ncol(new.data)>2) summary.data[j,2:ncol(summary.data)]=apply(new.data[,2:ncol(new.data)],2,count.records) else {
        summary.data[j,2:ncol(summary.data)] = count.records(new.data[,2])
      }
    }

    names(summary.data)[2:ncol(summary.data)] = names(data[2:ncol(summary.data)])
    names(summary.data)[1] = "min_rec_len"

    for(z in 1:nrow(summary.data)){
      out.data[x,z] = max(summary.data[z,2:ncol(summary.data)])
    }
  }
  
  if(length(record.length)==1){
    return(summary.data)
    break()
  }
  matplot(out.data,type="l",xlab="Maximum Years of Data",xaxt="n",ylab="Number of Lakes")
  axis(side=1,at=c(1:nrow(out.data)),labels=row.names(out.data))
  legend('topright',legend=missing.freq,lwd=2,col=c(1:length(missing.freq)),title="Freq. of Missing Data (%)")
  return(out.data)
}

#Example of running code
temp=TS.Length(data = dat,record.length = c(10:24),missing.freq = c(10,20,30,40))
#Run with a single length of time to extract table of when to start time series
temp=TS.Length(data = dat,record.length = 20,missing.freq = c(10,20,30,40))
#Look at output file. 530 lakes have 14-20 years of data (30% max missing data) starting with 1994.