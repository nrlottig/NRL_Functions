#Check to make sure you have all the needed packages installed
packages <- c("RCurl", "reshape2", "zoo","maps","maptools","reshape2")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

#platLAGOSLT.data extracts timeseries from LAGOS data that meet criteria set in function, plots the data, and returns a dataframe of the 
#mean summer values for the lakes with associated time series. The algorithm will automatically select the time period that has the most lakes 
#and if more than one time period has the same (max) number of lakes, it will select the most recent time period.
plotLAGOSLT.data = function(data,TSlength,miss.freq,locos,param.name){
  #Data columns (1) lagoslakeid, (2) sampledate [format="%m/%d/%Y"], (3) parameter value 
  #TSlength is length of time series desired
  #missing.freq  is amount of missing data within each TS as a percent (e.g., 20)
  #locos is LAGOS locus data table
  #param.name is string with parameter name for plotting purposes
  
  require(maps)
  require(maptools)
  require(reshape2)
  data.complete = data[complete.cases(data[,1:3]),]
  names(data.complete)[3] = "value"
  data.complete$DoY = as.numeric(strftime(as.Date(data.complete$sampledate,format="%m/%d/%Y"), format = "%j"))
  data.complete$Year = as.numeric(strftime(as.Date(data.complete$sampledate,format="%m/%d/%Y"), format = "%Y"))
  data.complete = data.complete[which(data.complete$DoY>=166 & data.complete$DoY<=258),] #subset to samples collected in summer
  data.complete = aggregate(data.complete$value,by=list(data.complete$lagoslakeid,data.complete$Year),FUN=mean) #estimate mean summer value
  names(data.complete) = c("lagoslakeid","year","value")
  years.data = range(data.complete$year)
  years.data = data.frame(year=seq(from=years.data[1],to=years.data[2],by=1),temp=NA)
  lakes = unique(data.complete$lagoslakeid)
  for(i in 1:length(lakes)){
    if(i == 1){
      temp.data = data.complete[which(data.complete$lagoslakeid==lakes[i]),]
      data.full = merge(temp.data,years.data,all=TRUE)
      data.full$lagoslakeid = lakes[i]
      data.full = data.full[,-4]
    } else {
      temp.data = data.complete[which(data.complete$lagoslakeid==lakes[i]),]
      temp.data = merge(temp.data,years.data,all=TRUE)
      temp.data$lagoslakeid = lakes[i]
      temp.data = temp.data[,-4]
      data.full = rbind(data.full,temp.data)
    }
    
  }
  data.wide <- dcast(data.full, lagoslakeid ~ year,mean)
  ts.out = TS.Length(data = data.wide,record.length = TSlength,missing.freq = miss.freq)
  start.year = as.numeric(names(ts.out)[which(ts.out==max(ts.out))[length(which(ts.out==max(ts.out)))]])
  years.data = seq(from=start.year,to=start.year+TSlength-1,by=1)
  
  #subset data to years of interest with correct record length
  data.full = data.full[which(data.full$year %in% years.data),]
  data.wide <- dcast(data.full, lagoslakeid ~ year,mean)
  data.wide$recordlength = rowSums(!is.na(data.wide[,2:(TSlength+1)]))
  data.wide = data.wide[which(data.wide$recordlength>=(TSlength-(TSlength*miss.freq/100))),]
  data.wide = merge(data.wide,locos[,c("lagoslakeid","nhd_lat","nhd_long")])
  map.regions = c('wisconsin','minnesota','vermont','maine','michigan','missouri','rhode island','new york','iowa','illinois','indiana','ohio','new hampshire','pennsylvania','connecticut','massachusetts','new jersey')
  map('state',region=map.regions,col=grey(.98),fill=TRUE,resolution = 0,mar=c(0,0,0,0),border=grey(.5),lty=5)
  points(data.wide$nhd_long,data.wide$nhd_lat,pch=16,col=rgb(0,0,0,.65),cex=.75)
  mtext(side=3,line=-1,paste(param.name," LT Data (",range(years.data)[1],":",range(years.data)[2],")",sep=""))
  mtext(side=3,line=-2,paste(nrow(data.wide)," Lakes w/ ",(TSlength-(TSlength*miss.freq/100))," - ",TSlength," yrs of data",sep=""))
  mtext(side=3,line=-3,paste(TSlength," yr record w/ max of ",miss.freq,"% missing data",sep=""))
  return(data.wide)
}
#Example
#out = plotLAGOSLT.data(data=LAGOS_Limno[,c("lagoslakeid","sampledate","tp")],TSlength = 25,miss.freq = 20,locos = LAGOS_Lakes,param.name = "Total P")

#TS.Length function for identifying optimal/maximum time series and time periods of record based on two input criteria
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
#Example data file take from Oliver et al. 2017 (chla data)
#require(Rcurl)
#dat = read.csv(text=getURL("https://raw.githubusercontent.com/nrlottig/NRL_Functions/master/datafiles/data_optimize_ts_window.csv"), header=T)
#temp=TS.Length(data = dat,record.length = c(10:24),missing.freq = c(10,20,30,40))
#Run with a single length of time to extract table of when to start time series
#temp=TS.Length(data = dat,record.length = 20,missing.freq = c(10,20,30,40))
#Look at output file. 530 lakes have 14-20 years of data (30% max missing data) starting with 1994.