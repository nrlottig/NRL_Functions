#Script to identify time windows and missing data rates for time series analysis of multiple different systems.
data = subset(LAGOS_summer_meanvals,LAGOS_summer_meanvals$variable=="chl_ugL")
data = data[,c(1:3)]
data3 = dcast(data,lagoslakeid~sampleyear)


TS.Length = function(data,record.length,missing.freq){
  library(reshape2)
  library(zoo)
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

    missing.freq.t = record.length[x] - (missing.freq/100*record.length[x])
    summary.data = data.frame(matrix(NA, nrow=length(missing.freq.t), ncol=ncol(new.data)))

    for(j in 1:length(missing.freq.t)){
      summary.data[j,1] = missing.freq.t[j]
      if(ncol(new.data)>2) summary.data[j,2:ncol(summary.data)]=apply(new.data[,2:ncol(new.data)],2,count.records) else {
        summary.data[j,2:ncol(summary.data)] = count.records(new.data[,2])
      }
    }

    names(summary.data)[2:ncol(summary.data)] = names(data[2:ncol(summary.data)])

    for(z in 1:nrow(summary.data)){
      out.data[x,z] = max(summary.data[z,2:ncol(summary.data)])
    }
  }
  matplot(out.data,type="l",xlab="Years of Data",xaxt="n",ylab="Number of Lakes")
  axis(side=1,at=c(1:nrow(out.data)),labels=row.names(out.data))
  legend('topright',legend=missing.freq,lwd=2,col=c(1:length(missing.freq)),title="Freq. of Missing Data (%)")
  if(length(record.length==1)){
    return(summary.data)
    break()
  }
  return(out.data)
}

temp=TS.Length(data = data3,record.length = c(10:24),missing.freq = c(10,20,30,40))
temp=TS.Length(data = data3,record.length = 20,missing.freq = c(10,20,30,40))
