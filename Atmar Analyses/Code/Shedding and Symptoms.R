create.average.df <- function(data) {
  avdata= data.frame("strains"=character(),
                     "virus.type"=character(),
                     
                     "cum.symp.onset"= integer(),
                     "symp.onset.num"= integer(),
                     
                     "cum.symp.dur"= integer(),
                     "symp.dur.num"= integer(),
                     
                     "cum.shed.onset"= integer(),
                     "shed.onset.num"= integer(),
                     
                     "cum.shed.dur"= integer(),
                     "shed.dur.num"= integer(),
                     
                     "cum.peak.vl.day"= integer(),
                     "peak.vl.day.num"=integer(),
                     
                     "size"=integer())
  
  cumul.data= cumData(data, avdata)
  avdata= averageData(cumul.data)
  
  return(avdata)
  
}

cumData <- function(data, avdata) {
  for (i in 1:nrow(data)) {
    
    row= data[i,]
    numRows= nrow(avdata)
    
    shed.onset= row$shed.onset
    shed.dur= row$shed.dur
    symp.onset= row$symptom.onset
    symp.dur= row$symptom.dur
    virus= row$virus
    type= row$virus.type
    peak.vl.day= row$peak.vl.day
    num.participants= row$participants
    
    if (!virus %in% avdata$strains) {
      avdata[numRows + 1, ]= c(virus, type, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    }
    
    rowStrain= which(avdata$strains== virus)
    
    #Symptom Onset
    avdata[rowStrain, ]$cum.symp.onset = 
      sum(c(as.numeric(avdata[rowStrain, ]$cum.symp.onset), num.participants * symp.onset), na.rm=TRUE)
    
    if(!is.na(symp.onset)) {
      avdata[rowStrain, ]$symp.onset.num = as.numeric(avdata[rowStrain, ]$symp.onset.num) + num.participants
    }
    
    #Symptom Duration
    avdata[rowStrain, ]$cum.symp.dur = 
      sum(c(as.numeric(avdata[rowStrain, ]$cum.symp.dur), num.participants * symp.dur), na.rm=TRUE)
    
    if(!is.na(symp.dur)) {
      avdata[rowStrain, ]$symp.dur.num = as.numeric(avdata[rowStrain, ]$symp.dur.num) + num.participants
    }
    
    #Shedding Duration
    avdata[rowStrain, ]$cum.shed.dur = 
      sum(c(as.numeric(avdata[rowStrain, ]$cum.shed.dur), num.participants * shed.dur), na.rm=TRUE)
    
    if(!is.na(shed.dur)) {
      avdata[rowStrain, ]$shed.dur.num = as.numeric(avdata[rowStrain, ]$shed.dur.num) + num.participants
    }
    
    #Shedding Onset
    avdata[rowStrain, ]$cum.shed.onset = 
      sum(c(as.numeric(avdata[rowStrain, ]$cum.shed.onset), num.participants * shed.onset), na.rm=TRUE)
    
    if(!is.na(shed.onset)) {
      avdata[rowStrain, ]$shed.onset.num = as.numeric(avdata[rowStrain, ]$shed.onset.num) + num.participants
    }
    
    #Peak Viral Load Day
    avdata[rowStrain, ]$cum.peak.vl.day = 
      sum(c(as.numeric(avdata[rowStrain, ]$cum.peak.vl.day), num.participants * peak.vl.day), na.rm=TRUE)
    
    if(!is.na(peak.vl.day)) {
      avdata[rowStrain, ]$peak.vl.day.num = as.numeric(avdata[rowStrain, ]$peak.vl.day.num) + num.participants
    }
    
    #Sample Size
  
    avdata[rowStrain, ]$size= as.numeric(avdata[rowStrain, ]$size) + 1
    
  }
  
  
  return(avdata)
}

averageData <- function(data) {
  
  data$symp.onset.avg= 0
  data$symp.dur.avg= 0
  data$shed.onset.avg= 0
  data$shed.dur.avg= 0
  data$peak.vl.day.avg= 0

  for (i in 1:nrow(data)) {
    row= data[i, ]
    
    data[i, ]$symp.onset.avg= as.numeric(row$cum.symp.onset) / as.numeric(row$symp.onset.num)
    data[i, ]$symp.dur.avg= as.numeric(row$cum.symp.dur) / as.numeric(row$symp.dur.num)
    data[i, ]$shed.onset.avg= as.numeric(row$cum.shed.onset) / as.numeric(row$shed.onset.num)
    data[i, ]$shed.dur.avg= as.numeric(row$cum.shed.dur) / as.numeric(row$shed.dur.num)
    data[i, ]$peak.vl.day.avg= as.numeric(row$cum.peak.vl.day) / as.numeric(row$peak.vl.day.num)
  }
  
  return(data)
}

eliminate.incomplete.data <- function(data, aggregated=TRUE) {
  if (aggregated) {
    sheddingonset <- data$shed.onset.avg
    sheddingduration <- data$shed.dur.avg
    symptomonset <- data$symp.onset.avg
    symptomduration <- data$symp.dur.avg
    strains <- data$strains
  } else {
    sheddingonset <- data$shed.onset
    sheddingduration <- data$shed.dur
    symptomonset <- data$symptom.onset
    symptomduration <- data$symptom.dur
    strains <- data$virus
  }
  
  data$complete= FALSE
  
  for (strain in 1:length(strains)) {
    
    plotQualifications= !(is.na(sheddingonset[strain]) || is.na(sheddingduration[strain]) || 
                            is.na(symptomonset[strain]) || is.na(symptomduration[strain]))
    
    if (plotQualifications) {
      
      data$complete[strain]= TRUE
      
    }
  }
  
  return(data[data$complete, ])
}

combineflu <- function(data, flu.name) {
  strains= data$strains
  
  flu.strains= grepl(flu.name, strains, fixed= TRUE)
  
  if(TRUE %in% flu.strains) {
    
    flu= c(flu.name, 'respiratory', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, TRUE)
    
    data= rbind(data, flu)
    
    flu.row= nrow(data)
    
    data[flu.row, ]$symp.onset.avg= weighted.mean(as.numeric(data[flu.strains, ]$symp.onset.avg), 
                                                    as.numeric(data[flu.strains, ]$symp.onset.num))
    data[flu.row, ]$shed.onset.avg= weighted.mean(as.numeric(data[flu.strains, ]$shed.onset.avg), 
                                                    as.numeric(data[flu.strains, ]$shed.onset.num))
    data[flu.row, ]$symp.dur.avg= weighted.mean(as.numeric(data[flu.strains, ]$symp.dur.avg), 
                                                  as.numeric(data[flu.strains, ]$symp.dur.num))
    data[flu.row, ]$shed.dur.avg= weighted.mean(as.numeric(data[flu.strains, ]$shed.dur.avg), 
                                                  as.numeric(data[flu.strains, ]$shed.dur.num))
    data[flu.row, ]$peak.vl.day.avg= weighted.mean(as.numeric(data[flu.strains, ]$peak.vl.day.avg), 
                                             as.numeric(data[flu.strains, ]$peak.vl.day.num))
    data[flu.row, ]$size= sum(as.numeric(data[flu.strains, ]$size))
    
  }
  
  return(data[!flu.strains, ])
}

plot.RespiratoryVsDigestive <- function(data) {
  
  data= data[order(data$virus.type), ]
  
  data= combineflu(data, 'Influenza A')
  data= combineflu(data, 'Influenza B')
  
  plot.SheddingVsSymptoms(data, split=TRUE)
}

plot.SheddingVsSymptoms <- function(data, split=FALSE) {
  
  if(split) {
    par(bty='n', mfrow=c(1, 2), mar=c(5, 10, 2, 2))
    plot.panel(data, 'respiratory')
    mtext('Respiratory Viruses', side=3, font=2)
    
    plot.panel(data, 'digestive')
    mtext('Digestive Viruses', side=3, font=2)
    
  } else {
    par(bty='n', mar=c(5, 14, 2, 2))
    plot.panel(data)
    mtext('Shedding and Symptom Duration', font=2)
  }
  
  legend(15, 2, legend=c("Symptoms", "Shedding"), 
         col=c("#D55E00", "#009E73"), lty=2:1, lwd=2)
  
}

plot.panel <- function(data, type=FALSE) {
  
  size= data$size
  sheddingonset= data$shed.onset.avg
  sheddingduration= data$shed.dur.avg
  symptomonset= data$symp.onset.avg
  symptomduration= data$symp.dur.avg
  strains= data$strains
  
  if(typeof(type) == 'character') {
    complete= data$complete==TRUE & data$virus.type==type
  } else {
    complete= data$complete==TRUE
  }
  
  completeStrains= strains[complete]
  
  plot(0, 0, xlim=c(0, 20), ylim=c(0, length(completeStrains)-1), col = "white", 
       xlab= "Days Post-Inoculation", ylab='', yaxt='n')
  
  axis(2, at=seq(0, length(completeStrains)-1, 1), las=2, labels=completeStrains)
  
  plot.segments(sheddingonset[complete], sheddingduration[complete], FALSE, TRUE, size[complete])
  
  plot.segments(symptomonset[complete], symptomduration[complete], TRUE)
  
}

plot.segments <- function(xstarts, xlengths, aggregated=TRUE, complete=NA, dotted, 
                          includesize=FALSE, samplesize=0) {
  if (dotted) {
    color="#D55E00"
    lty = 2
    offsety= 0.5
  } else {
    lty = 1
    color="#009E73"
    offsety= 0
  }
  
  for (i in 1:length(xstarts)) {
    
    dot <- FALSE
    redcol <- FALSE
    
    xstart <- as.numeric(xstarts[i])
    xlen <- as.numeric(xlengths[i])
    
    if(is.na(xstart)) {
      start <- 0
      redcol <- TRUE
    } else {
      start <- xstart
    }
    
    if(is.na(xlen)) {
      dot <- TRUE
    } else {
      length <- xlen
    }
    
    # If there's missing duration data, just make a point for onset
    if (dot) {
      points(start, i-1 + offsety, col=color, cex=2)
    } else {
      end <- start + length
      segments(start, i-1 + offsety, end, i-1 + offsety, lty=lty, 
               col=ifelse(redcol, "#D41159", color), cex=2)
    }
    
    
    # Star the sets of lines with complete data
    if (!aggregated) {
      if (complete[i] && !dotted) {
        points(end + 3, i-1 + offsety, cex=2, pch=8)
      }
    } 
    
    if (includesize) {
      text(end + 3, i-1, labels=paste('n = ', samplesize[i], sep=""), cex=0.8)
    }
  }
  
}

#Import and extract data from csv/Excel sheet

data <- read.csv('ViralCHIStudies.csv')

#Plot shedding and symptom timelines!
strains <- data$virus

ordered.strains <- rev(order(strains))

ordered.data <- data[ordered.strains,]
norwalk <- ordered.data$virus=="Norwalk virus"
data.norwalkatbottom <- rbind(ordered.data[norwalk, ], ordered.data[!norwalk, ])

avg.data <- eliminate.incomplete.data(create.average.df(ordered.data))

png("Shedding and Symptom Duration.png", width = 700, height = 480)
plot.SheddingVsSymptoms(avg.data)
dev.off()

png("Respiratory vs. Digestive Timing.png", width = 800, height = 550)
plot.RespiratoryVsDigestive(avg.data)
dev.off()

# Making Shedding and Symptom Duration plot with each study's data laid out
numrows <- nrow(data.norwalkatbottom) - 1

# Check for complete dataset with 4 points
data.norwalkatbottom$complete <- !(is.na(data.norwalkatbottom$symptom.onset -
                                           data.norwalkatbottom$symptom.dur -
                                           data.norwalkatbottom$shed.onset -
                                           data.norwalkatbottom$shed.dur))
{
  postscript(file="Complete Shedding Vs. Symptom.eps", paper="special", 
             horizontal=FALSE, height=11, width=14)
  
  par(bty='n', mar=c(5, 15, 2, 2))
  
  plot(0, 0, xlim=c(0, 40), ylim=c(0, numrows), col = "white", 
       xlab= "Days Post-Inoculation", ylab='', yaxt='n', cex.axis=1.5, 
       cex.lab=1.5)
  
  axis(2, at=seq(0, numrows, 1), las=2, labels=c("Norwalk Virus", "", "",
                                                 "Snow Mountain Virus",
                                                 "RSV", "",
                                                 "Rhinovirus", "", "", "", "", "",
                                                 "Influenza B", "",
                                                 "Influenza A",
                                                 "", "", "", "", "", "", "", 
                                                 "", "", "", "", ""
                                                 ), cex.axis=1.5)
  
  plot.segments(data.norwalkatbottom$shed.onset, 
                data.norwalkatbottom$shed.dur, aggregated=FALSE, 
                complete=data.norwalkatbottom$complete,
                dotted=FALSE)
  
  plot.segments(data.norwalkatbottom$symptom.onset, 
                data.norwalkatbottom$symptom.dur, aggregated=FALSE, 
                complete=data.norwalkatbottom$complete,
                dotted=TRUE)
  
  mtext('Shedding and Symptom Duration', font=2, cex=2)
  
  
  legend(18, 10, legend=c("Symptoms", "Shedding", "No Time Of Onset Reported", 
                         "No Duration Reported", "Complete Dataset"), 
         col=c("#E1BE6A", "#009E73", "#D55E00", "black", "black"), 
         lty=c(2, 1, 1, NA, NA), 
         pch=c(NA, NA, NA, 1, 8), lwd=2, cex=1.5)
  
  dev.off()
}
