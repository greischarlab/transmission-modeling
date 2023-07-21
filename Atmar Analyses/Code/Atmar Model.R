# require(metaDigitise)
# 
# atmardata <- metaDigitise("~/Research Project/Summer URE/Process Images",
#                             summary=FALSE)
# 
# atmardata= atmardata$scatterplot$plotc.jpg

setwd("C:/Users/ktopa/OneDrive/Documents/Research Project/Summer URE")

require(ggplot2)

generate.estimates <- function(params, t) {
  
  p1= exp(params[1])
  p2= exp(params[2])
  p3= exp(params[3])
  p4= exp(params[4])
  
  estimate= (2 * p1) / (exp(-p2*(t-p3)) + exp(p4*(t-p3)))
  
  return(log10(estimate))
  
} 

sse.viralLoad <- function(params, data) {
  vl = log10(data$y)
  t = data$x * 24 - 19
  
  estimates= generate.estimates(params, t)
  
  sse= sum((estimates - vl)^2)
  
  if (sse > 1e20) {
    return(1e20)
  } else {
    return(sse)
  }
  
}

get.params <- function(data) {
  
  fits= generate.possible.params(data)
  
  minsse= min(fits[fits$convergence==0, ]$sse)
  
  pars= unlist((fits[fits$sse == minsse, ])[1:4])
  
  return (c(pars, minsse))
}

generate.possible.params <- function(data) {
  numguess= 1e4
  
  fits= data.frame(matrix(NA, nrow = numguess, ncol = 6))
  
  colnames(fits)= c('p1','p2','p3','p4','convergence','sse')
  
  for (g in 1:numguess) {
    p1= runif(1, -1, log(1e12))
    p2= runif(1, -1, 1)
    p3= runif(1, -1, 4*24)
    p4= runif(1, -5, 1)
    
    params= c(p1, p2, p3, p4)
    
    fit= optim(params, sse.viralLoad, data=data)
    
    fits[g, ]= c(fit$par, fit$convergence, fit$value)
  }
  
  return(fits)
}

plot.optimization <- function(data) {
  fits= generate.possible.params(data)
  t=data$x * 24 - 19
  vl=log10(data$y)
  
  sortedfits= fits[order(fits$sse),]
  
  colors= c()
  for(i in c(800, 500, 200, 150, 100, 50, 1)) {
    fit=sortedfits[i, ]
    params= c(fit$p1, fit$p2, fit$p3, fit$p4)
    predictedvl= generate.estimates(params, t)
    plot(t, vl, xlim=c(0, 260), ylim=c(0, 12), xlab="Hours Since Inoculation", 
         ylab="Viral Shedding (log10)")
    lines(t, predictedvl, col="navyblue")
    legend(200, 2, legend=c("Fit"), col="navyblue", lty=1)
  }
}

plot.viralload <- function(data, pars, sympdata, max.modeltime, ymax, cex,
                           aggregated=FALSE, labels=TRUE) {
  
  if(data$dose[1] == '4800') {
    col="#E69F00"
  } else if (data$dose[1] == '48') {
    col="#871a6e" 
  } else if (data$dose[1] == '4.8') {
    col="#009E73"
  }
  
  t=data$x * 24 - 19
    
  vl= data$y
  
  if(aggregated) {
    points(t, log10(vl), col=col)
  } else {
    plot(t, log10(vl), col=col, xlim=c(0, 15 * 24 - 19), ylim=c(7, ymax),
       xlab="", ylab="", xaxt="n", yaxt="n", cex.lab=cex, cex.axis=cex, 
       cex.main=cex, font.lab=2, cex=cex*2, lwd=cex-2)
    
    linet =seq(min(data$x)*24-19, data[max.modeltime,]$x*24-19, by=1)
    predictedvl= generate.estimates(pars, linet)
    
    symp.onset= mean(sympdata$symp.onset) * 24 - 19
    symp.end= symp.onset + mean(sympdata$symp.dur) * 24
    plot.bestfit(linet, predictedvl, col, symp.onset, symp.end, cex)
    
    axis(1, cex.axis=cex, at=seq(0, 300, by=50), labels=c(expression(0),
                                                          expression(50),
                                                          expression(100),
                                                          expression(150),
                                                          expression(200),
                                                          expression(250),
                                                          expression(300)),
         lwd=0, line=6)
    axis(2, cex.axis=cex, at=seq(6, 14, by=1), 
         labels=c(expression(10^6), 
                  expression(10^7),
                  expression(10^8),
                  expression(10^9),
                  expression(10^10),
                  expression(10^11),
                  expression(10^12),
                  expression(10^13),
                  expression(10^14)),
         lwd=0, line=3, las=2)
    
    if (time.peak(linet, predictedvl) < symp.onset) {
      prepost="Pre-symptomatic"
    } else {
      prepost="Post-symptomatic"
    }
    
    mtext(text=paste("Participant ", data$id[1], " (",
                     dose, " RT-PCR Units): ", prepost,
                     sep=""), side=3, line=6, cex=cex - 2, font=2)
    
    mtext(text="Viral Load", side=2, line=20, cex=cex-2, font=2)
    
    if (labels) {
    mtext(text="Hours After Inoculation", side=1, line=17, cex=cex-2, font=2)
    }
    plot.symptoms(sympdata, dose, cex)
  }
  
}

plot.bestfit <- function(linet, predictedvl, col, symp.onset, symp.end, cex) {
  
  lines(linet, predictedvl, col=col, lwd=cex)
  
  lines(linet[linet <= symp.end & linet >= symp.onset], 
        predictedvl[linet <= symp.end & linet >= symp.onset], 
        col=col, lwd=cex*3)
  
  if (time.peak(linet, predictedvl) < symp.onset) {
    redpointstyle <- 8
  } else {
    redpointstyle <- 19
  }
  
  points(time.peak(linet, predictedvl), max(predictedvl), pch=redpointstyle, 
         col="red", cex=cex*3, lwd=cex*2)
}

plot.symptoms <- function(data, dose, cex) {
  if(dose == 4800) {
    col="#E69F00"
  } else if (dose == 48) {
    col="#871a6e" 
  } else if (dose == 4.8) {
    col="#009E73"
  }
  
  #Study day 1 begins 5-6 hours postinoculation
  start= data$symp.onset * 24 - 19
  length= data$symp.dur * 24
  
  abline(v=c(start, start + length), lty=2, col=col, lwd=cex-2)
}

time.peak <- function(time, vl) {
  return(time[which(vl==max(vl))])
}

delay.symponset.peakvl <- function(time, predictedvl, symp.onset) {
  peak.vl.time= time.peak(time, predictedvl)
  
  delay= peak.vl.time - (symp.onset * 24 - 19)
  
  return(delay)
}

plot.parameter.panel <- function(id, param.value, param.name, cex, lab) {
  linecolor <- "orchid"
  slopex <- 3
  labx <- 1
  laby <- 17
  cex <- 2
  lwd <- 3
  bold <- 2
  xlab <- expression(paste(log[10], " Infectious Dose"))
  
  plot(id, param.value, xlab="", ylab="", 
       cex.lab=cex, cex.axis=cex, cex=cex, font.lab=bold)
  abline(lm(param.value ~ id), col=linecolor, lwd=lwd)
  mtext(text=xlab, side=1, line=4, cex=cex)
  mtext(text=param.name, side=2, line=3, cex=cex)
}

first.min <- function(data) {
  
  t.peak <- which(data$y == max(data$y))
  
  points.after.peak <- data$y[-1:-t.peak]
  
  next.points <- append(data$y[-1:-(t.peak+1)], NA)
  
  return.val <- which(points.after.peak <= next.points)[1] + t.peak
  
  if (is.na(return.val)) {
    
    return(length(data$y))
    
  } else {
    
    return(return.val)
    
  }
}

make.param.df <- function(id, data, maxt) {
  
  patientparams <- matrix(nrow=length(id), ncol=6)
  
  for (i in 1:length(id)) {
    params <- get.params(data[data$id==id[i],]
                         [1:(maxt[maxt$id==id[i],]$max.point),])
    
    patientparams[i,] <- c(id[i], params[1:4],
                           data[data$id==id[i],]$dose[1])
  }
  
  patientparams <- as.data.frame(patientparams)
  colnames(patientparams) <- c("patient", "p1", "p2",
                               "p3", "p4", "dose")
  
  return(patientparams)
  
}

atmardata <- read.csv('atmardata.csv')
atmarasymp <- read.csv("atmardataasymp.csv")

# Adding the max time index that the viral trajectory model should be fit to

id <- unique(atmardata$id)

# Method 1: Stopping at the first local minimum
maxt.firstmin <- data.frame(id=id, max.point=rep(NA, length(id)))

for (p in id) {
  maxt.firstmin[maxt.firstmin$id==p,]$max.point <-  
    first.min(atmardata[atmardata$id==p,])
}

# Method 2: Stopping when viral load gets back down to the initial viral load
maxt.backtoinitial <- data.frame(id=id, max.point=rep(NA, length(id)))

for (p in id) {
  participantdat <- atmardata[atmardata$id==p,]
  
  t.peak <- which(participantdat$y == max(participantdat$y))
  
  points.after.peak <- participantdat$y[-1:-t.peak] 
  
  maxt.backtoinitial[maxt.backtoinitial$id==p,]$max.point <-  
    which(points.after.peak <= participantdat$y[1])[1] + t.peak
}

# Producing parameters with both methods of cutoff

patientparams.localmin <- make.param.df(id, atmardata, maxt.firstmin)

#write.csv(patientparams.localmin, 
#          "Atmar Patient Parameters local minimum.csv", row.names=FALSE)

patientparams.backtoinitial <- make.param.df(id, atmardata, maxt.backtoinitial)

write.csv(patientparams.backtoinitial, 
          "Atmar Patient Parameters back to initial.csv", row.names=FALSE)
patientparams.backtoinitial <- read.csv("Atmar Patient Parameters back to initial.csv")
# Importing symptom data and parameters

atmarsymptoms <- read.csv('atmarsymptoms.csv')
patientparams <- read.csv("Atmar Patient Parameters local minimum.csv")

#sse= data.frame("id"=id, "avg.sse"=rep(0, length(id)))

delay <- data.frame("id"=id, "delay"=rep(0, length(id)), 
                    "dose"=rep(0, length(id)))


#patients <- data.frame(matrix(NA, nrow = length(id), ncol = 7))
#colnames(patients)= c('patient','infectious dose','p1','p2','p3','p4', 'sse')
row <- 1

for (i in 1:length(id)) {
  
  patientid <- id[i]
  vldata <- atmardata[atmardata$id == patientid, ]
  sympdata <- atmarsymptoms[atmarsymptoms$id == patientid, ]
  dose <- vldata$dose[1]
  
  params <- unlist(patientparams[patientparams$patient==patientid, 2:5])

  #patients[row, ] = c(patientid, vldata[1, ]$dose, params)
  row <- row + 1
  #minsse= params[5]

  t=seq(min(vldata$x)*24-19, 
        vldata[maxt.firstmin[maxt.firstmin$id==patientid,]$max.point,]$x*24-19, by=1)
  #sse[sse$id==patientid, ]$avg.sse= minsse / length(t)
  delay[delay$id==patientid, ]$delay <-
     delay.symponset.peakvl(t, generate.estimates(params, t), sympdata$symp.onset)
   
  delay[delay$id==patientid, ]$dose <- dose
  
  #Plotting data points, Li and Handel/Holder model fit, and symptom onset/duration
  #png("Atmar703.png", width=800, height=800)
  #par(mar=c(5, 5, 4, 4))
  #plot.viralload(vldata, unlist(params), sympdata, 
  #maxt[maxt$id==patientid,]$max.point, 2, 14, labels=TRUE)
  #text(280, 11, labels="C", cex=2, font=2)
  #dev.off()
}

write.csv(delay, "Delays by ID Norovirus.csv")

legend(220, 11.3, 
       legend=c("4800 RT-PCR Units", "48 RT-PCR Units", "4.8 RT-PCR Units"), 
       lty= 1, col=c("#E69F00", "#871a6e", "#009E73"), lwd=3, cex=1.5)

dev.off()

infectious.dose <- rev(unique(delay$dose))
average.delay <- c()
prop.presymp <- c()

for (i in 1:length(infectious.dose)) {
  dose.info <- delay[delay$dose==infectious.dose[i],]
  average.delay[i] <- mean(dose.info$delay)
  prop.presymp[i] <- sum(dose.info$delay < 0) / nrow(dose.info)
}

write.csv(data.frame(infectious.dose = infectious.dose, 
                     average.delay = average.delay,
                     prop.presymp = prop.presymp), 
          "Average Delays Norovirus.csv")
#Plotting average viral shedding trajectory over time for each infectious dose
png("Atmar Modeling Averages.png", width=900, height=800)

par(mar=c(5.1, 6.1, 4.1, 2.1))
plot(1, type="n", xlab='', ylab='', xlim=c(0, 336), ylim=c(6,12), cex.axis=2)

mtext(text="Hours After Inoculation", side=1, line=4, cex=2, font=2)
mtext(text=expression(bold(paste("Viral Shedding (", log[10], ")"))), side=2, line=3, cex=2)

params4800 <- read.csv("Params 4800.csv")$x[1:4]
params48 <- read.csv("Params 48.csv")$x[1:4]
params4.8 <- read.csv("Params 4.8.csv")$x[1:4]

vldata4800 <- atmardata[atmardata$dose=='4800',]
vldata48 <- atmardata[atmardata$dose=='48',]
vldata4.8 <- atmardata[atmardata$dose=='4.8',]

sympdata4800 <- atmarsymptoms[atmarsymptoms$dose=='4800',]
sympdata48 <- atmarsymptoms[atmarsymptoms$dose=='48',]
sympdata4.8 <- atmarsymptoms[atmarsymptoms$dose=='4.8',]

plot.viralload(vldata4800, params4800, sympdata4800)

plot.viralload(vldata48, params48, sympdata48)

plot.viralload(vldata4.8, params4.8, sympdata4.8)

legend(220, 11.3, 
       legend=c("4800 RT-PCR Units", "48 RT-PCR Units", "4.8 RT-PCR Units"), 
       lty= 1, col=c("#E69F00", "#871a6e", "#009E73"), lwd=3, cex=1.5, text.font=2)

dev.off()
#Plotting the parameters against infectious dose
id <- log10(patientparams$infectious.dose)
p1 <- (patientparams$p1)
p2 <- (patientparams$p2)
p3 <- (patientparams$p3)
p4 <- (patientparams$p4)

png("Parameters vs. Infectious Dose.png", width = 900, height = 800)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 1, 2.1))

plot.parameter.panel(id, p1, "P1", cex, "A)")
plot.parameter.panel(id, p2, "P2", cex, "B)")
plot.parameter.panel(id, p3, "P3", cex, "C)")
plot.parameter.panel(id, p4, "P4", cex, "D)")

dev.off()

#Plotting patients 703, 722, and 724
patientexamples <- c('703', '722', '724')
label <- c("A)", "B)", "C)")

postscript(file="Example Atmar.eps", paper="special", 
           height=75, width=30, horizontal=FALSE)
par(mfrow=c(3, 1), mar=c(5.1, 8.1, 4.1, 3.1), oma=c(10, 10, 7, 5),
    mai=c(1.1, 2.6, 3.1, 1.1))
for(i in 1:length(patientexamples)) {

  p=patientexamples[i]
  vldata= atmardata[atmardata$id == p, ]
  sympdata= atmarsymptoms[atmarsymptoms$id == p, ]
  dose=vldata$dose[1]
  params=unlist(patientparams[patientparams$patient==p,][2:5])
  
  if(i == 3) {
    plot.viralload(vldata, params, sympdata, maxt.firstmin[maxt.firstmin$id==p,]$max.point, 12, 
                   cex=8, aggregated=FALSE)}
  else {
    plot.viralload(vldata, params, sympdata, maxt.firstmin[maxt.firstmin$id==p,]$max.point, 12,
                   cex=8, aggregated=FALSE, 
                   labels=FALSE)
  }
  
  text(320, 11.6, labels=label[i], font=2, cex=11)

}
dev.off()

plot.patient.forloop <- function(data, patientid, params, maxt) {
  vldata <- data[data$id==patientid,]
  sympdata <- atmarsymptoms[atmarsymptoms$id == patientid, ]
  dose <- vldata$dose[1]
  pars <- unlist(params[params$patient==patientid,2:5])
  
  plot.viralload(vldata, pars, sympdata, 
                 maxt[maxt$id==patientid,]$max.point,
                 14, cex=8, aggregated=FALSE, labels=TRUE)
}

# Making all the plots for each Atmar participant
for(patientid in id) {
  postscript(file=paste(as.character(patientid), ".eps"), paper="special", 
             height=30, width=32, horizontal=FALSE)
  par(mar=c(5.1, 8.1, 1.1, 1.1), oma=c(15, 15, 1, 5),
      mai=c(1.1, 2.6, 3.1, 1.1))
  plot.patient.forloop(atmardata, patientid, patientparams, maxt.firstmin)
    dev.off()
}

# Grouping plots by infectious dose
for (dose in 48) {
  patients <- unique(atmardata[atmardata$dose == dose,]$id)
  nrows <- ceiling(length(patients) / 3)
  
  postscript(file=paste(as.character(dose), ".eps"), paper="special", 
             height=28 * nrows, width=28 * 3, horizontal=FALSE)
  par(mfrow=c(nrows, 3), mar=c(5.1, 8.1, 1.1, 1.1), oma=c(15, 15, 1, 5),
      mai=c(2.1, 2.6, 3.1, 1.1))
  
  for (patientid in patients){
    plot.patient.forloop(atmardata, patientid, patientparams.backtoinitial, 
                         maxt.backtoinitial)
  }
  
  dev.off()
}
#Finding out if peak viral load and other parameters significantly differ between doses

peakvl <- data.frame(matrix(NA, nrow = length(unique(atmardata$id)), ncol = 4))
colnames(peakvl)= c('patient','dose', 'peakvl', 'peaktime')
row <- 1

for (patient in unique(atmardata$id)) {
  data= atmardata[atmardata$id==patient,]
  linet =seq(min(data$x)*24-19, 
             (data[maxt[maxt$id==patient,]$max.point,]$x)*24-19, by=1)
  params=unlist(patientparams[patientparams$patient==patient,][2:5])
  predictedvl= generate.estimates(params, linet)
  peak.vl = max(predictedvl)
  
  peakvl[row,]$patient= patient
  peakvl[row,]$peakvl= peak.vl
  peakvl[row,]$dose= data$dose[1]
  peakvl[row,]$peaktime= linet[which(predictedvl==peak.vl)]
  row = row + 1
}

#peakvl <- write.csv(peakvl, "peakvl.csv")

peakvl <- read.csv("peakvl.csv")

#looking at peak vl predicted by the model
peakvl.aov <- summary(aov(peakvl$peakvl ~ peakvl$dose))

#looking at peak vl time predicted by the model
peaktime.aov <- summary(aov(peakvl$peaktime ~ peakvl$dose))

#looking at p1 parameter
p1.aov <- summary(aov(patientparams$p1 ~ patientparams$dose))

#looking at p2 parameter
p2.aov <- summary(aov(patientparams$p2 ~ patientparams$dose))

#looking at p3 parameter
p3.aov <- summary(aov(patientparams$p3 ~ patientparams$dose))

#looking at p4 parameter
p4.aov <- summary(aov(patientparams$p4 ~ patientparams$dose))

#Variance of peak vl between infectious doses
lowdose.peakvl.var <- var(peakvl[peakvl$dose=="4.8",]$peakvl)
meddose.peakvl.var <- var(peakvl[peakvl$dose=="48",]$peakvl)
highdose.peakvl.var <- var(peakvl[peakvl$dose=="4800",]$peakvl)

#Variance of peak vl time between infectious doses
lowdose.peakvltime.var <- var(peakvl[peakvl$dose=="4.8",]$peaktime)
meddose.peakvltime.var <- var(peakvl[peakvl$dose=="48",]$peaktime)
highdose.peakvltime.var <- var(peakvl[peakvl$dose=="4800",]$peaktime)

#Saving variance to csv file
variance <- data.frame(dose=c(4.8, 48, 4800), 
                       peakvl.var=c(lowdose.peakvl.var, 
                                    meddose.peakvl.var, 
                                    highdose.peakvl.var),
                       peaktime.var=c(lowdose.peakvltime.var, 
                                      meddose.peakvltime.var, 
                                      highdose.peakvltime.var))

write.csv(variance, "Variance in peakvl atmar.csv")

#Looking at significance of variance of peakvl and peakvltime between infectious doses
variance <- read.csv("Variance in peakvl atmar.csv")

peakvl.var.aov <- summary(aov(variance$peakvl.var ~ variance$dose))
peakvltime.var.aov <- summary(aov(variance$peaktime.var ~ variance$dose))

#Making Swarm plot of symptom onset

ggplot(atmarsymptoms, aes(x=as.factor(dose), y=symp.onset, group=as.factor(dose))) +
  geom_boxplot() + xlab("Infectious Dose") + ylab("Day of Symptom Onset")

ggsave("Symptom Onset Box Plot.png")
  
# Seeing if pathogen load correlates to time of symptom onset

# Eliminating patients with no viral load projected at symptom onset

atmarsymptoms <- atmarsymptoms[atmarsymptoms$id != 731,]
atmarsympmtomaticdata <- atmardata[atmardata$id != 731,]
patientparamsymptomaticdata <- patientparams[patientparams$patient != 731,]

vl.atsymponset <- data.frame(dose= atmarsymptoms$dose, 
                             symp.onset=rep(NA, length(atmarsymptoms$id)),
                             vl=rep(NA, length(atmarsymptoms$id)),
                             presymp= rep(NA, length(atmarsymptoms$id)))

for (i in 1:length(atmarsymptoms$id)) {
  patientid <- atmarsymptoms$id[i]
  patientdat <- atmarsympmtomaticdata[atmarsympmtomaticdata$id==patientid,]
  patientpars <- unlist(patientparamsymptomaticdata[
    patientparamsymptomaticdata$patient==patientid, 2:5])
  patienttime <- seq(min(patientdat$x)*24-19, 
                     patientdat[maxt.firstmin[maxt.firstmin$id==patientid,]$max.point,]$x*24-19, 
                     by=1)
  vltrajectory <- generate.estimates(patientpars, patienttime)
  
  vl.atsymponset[i,]$symp.onset <- 
    atmarsymptoms[atmarsymptoms$id==patientid,]$symp.onset
  
  if (time.peak(patienttime, vltrajectory) < vl.atsymponset[i,]$symp.onset*24-19) {
    vl.atsymponset[i,]$presymp <- TRUE
  } else {
    vl.atsymponset[i,]$presymp <- FALSE
  }
  
  symp.onset.point <- 
    which.min(abs(patienttime - 
                    (atmarsymptoms[atmarsymptoms$id==patientid,]$
                       symp.onset*24 - 19)))
  vl.atsymponset[i,]$vl <- vltrajectory[vl.atsymponset[i,]$symp.onset]
}

write.csv(vl.atsymponset, "Viral Load at Symptom Onset.csv", row.names=FALSE)

vl.atsymponset <- read.csv("Viral Load at Symptom Onset.csv")

